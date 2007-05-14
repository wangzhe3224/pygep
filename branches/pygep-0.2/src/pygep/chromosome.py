# PyGEP: Gene Expression Programming for Python
# Copyright (C) 2007  Ryan J. O'Neil
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
from itertools import izip
from pygep.functions.linkers import default_linker
from pygep.gene import KarvaGene
from pygep.util import cache
import functools, random



def symbol(s):
    '''
    Decorator that assigns a symbol to a function for chromosome display.
    The symbol is stored in the function.symbol attribute.

        @symbol('/')
        def divide(x, y):
            return x / y
    '''
    def decorator(func):
        func.symbol = s
        return func

    return decorator


class MetaChromosome(type):
    '''
    Metaclass for computing various information about a chromosomal
    type.  Sets the following attributes on a chromosome class:
        - arity:   maximum functional arity
        - symbols: symbols that can reside in the head
    Also turns caching of fitness values on for all chromosomes.
    '''
    def __new__(typ, *args, **kwds):
        t = type.__new__(typ, *args, **kwds)
        t.symbols = t.functions + t.terminals

        # Find the max arity
        try:
            t.arity = max([f.func_code.co_argcount for f in t.functions])
        except ValueError:
            t.arity = 0

        # Cache fitness values
        t._fitness = cache(t._fitness)
        return t


class Chromosome(object):
    '''
    A Chromosome must provide these attributes:
        - functions: tuple of nonterminals
        - terminals: tuple of terminal symbols

    And override these functions:
        - _fitness: fitness of a given individual
        - _solved:  True if the problem is optimally solved (optional)

    An example Chromosome that evolves simple arithmetic expressions
    on data objects providing attributes 'a' and 'b' and the constants
    1 and 2:

        from pygep.functions.arithmetic import *
        from pygep import Chromosome

        class Calculation(Chromosome):
            functions = multiply, add, subtract, divide
            terminals = 'a', 1, 2

            def _fitness(self):
                # Evaluate chromosome fitness here.
                # This often involves calling self.evaluate(something)

            def _solved(self):
                # Not required, but useful if the problem can
                # be optimally solved.  Usually this just means
                # checking self.fitness.
    '''
    __metaclass__ = MetaChromosome
    __next_id = 1
    gene_type = KarvaGene


    functions = ()
    terminals = ()
    head = tail = length = 0


    @classmethod
    def generate(cls, head, genes=1, linker=default_linker):
        '''
        Returns a generator of random GEP chromosomes
        @param head:   head length (min=0)
        @param genes:  number of genes (min=1)
        @param linker: linking function
        '''
        tail = head * (cls.arity - 1) + 1

        while True:
            g = [None] * genes
            for i in xrange(genes):
                g[i] = cls.gene_type(
                    [random.choice(cls.symbols)   for _ in xrange(head)] + \
                    [random.choice(cls.terminals) for _ in xrange(tail)], head
                )

            yield cls(g, head, linker)


    def __init__(self, genes, head, linker):
        '''
        Instantiates a chromsome instance and analyzes it for evaluation.
        Sets the self.coding tuple to the last genes in the coding regions
        and various other internal data for the chromosome.  Note that it
        is generally unwise to instantiate chromosomes manually.  It is
        much more common to create them via calls to the static method
        Chromosome.generate(...).

        @param genes:  number of genes in the chromosome (min=1)
        @param head:   length (not index) of the gene heads (min=0)
        @param linker: linker function for gene evaluation
        '''
        # Must have at least one gene and a head length of 0
        if head < 0:
            raise ValueError('Head length must be at least 0')
        if not genes:
            raise ValueError('Must have at least 1 gene')

        self.genes  = genes
        self.head   = head
        self.linker = linker

        # Unique number of the organism
        self.__id = type(self).__next_id
        type(self).__next_id += 1


    def __len__(self):
        return sum(len(g) for g in self.genes)


    @cache
    def __repr__(self):
        return ''.join(repr(g) for g in self.genes)


    def _child(self, genes):
        '''Returns a child chromosome of self'''
        return type(self)(genes, self.head, self.linker)


    # Unique ID of the organism
    id = property(lambda self: self.__id, doc='Organism #')


    def __call__(self, obj):
        '''
        Evaluates a given GEP chromosome against some instance.  The
        terminals in the chromosome are assumed to be attributes on
        the object instance provided (unless they are numeric constants).

        @param obj: an object instance with terminal attributes set
        @return:    result of evaluating the chromosome
        '''
        return self.linker(*(g(obj) for g in self.genes))


    def _fitness(self):
        '''Abstract base method for fitness determination'''
        raise NotImplementedError('Must override Chromosome._fitness')

    def _solved(self):
        '''Defaults to False.  Override to terminate early.'''
        return False

    fitness = property(lambda self: self._fitness(), doc='Fitness value')
    solved  = property(lambda self: self._solved(),  doc='Problem solved')


    def mutate(self, rate):
        '''
        Produces a new chromosome via potential point mutation on each
        locus.  If nothing changes, the original chromosome is returned.

        @param rate: mutation rate per locus
        @return: child chromosome (or self)
        '''
        genes = list(self.genes)
        
        # Traverse the chromosome gene by gene
        for gene_idx, gene in enumerate(self.genes):
            # Then locus by locus
            replacements = []
            for i, allele in enumerate(gene):
                # Do we mutate this locus?
                if random.random() < rate:
                    # Mutation within the tail can only use terminals
                    if i >= self.head:
                        new_allele = random.choice(self.terminals)
                    else:
                        new_allele = random.choice(self.symbols)
                    
                    # Only use this if the mutation actually did something
                    if new_allele != allele:
                        replacements.append((i, [new_allele]))


            # If we have actual replacements to make, do them
            if replacements:
                genes[gene_idx] = gene.derive(replacements)
            
        # Create a child of this chromosome
        if genes != self.genes:
            return self._child(genes)
        
        return self


    def invert(self):
        '''Produces a new chromosome via head inversion'''
        if self.head < 2: # Head inversion does nothing in this case
            return self

        chromosome = list(self.chromosome)

        # Choose a random gene and two points within the head
        gene = random.choice(self._gene_starts)
        start, stop = random.sample(xrange(self.head), 2)

        # Order the indexes correctly
        if start > stop:
            start, stop = stop, start

        # Create the new chromosome
        chromosome[start:stop] = reversed(chromosome[start:stop])
        if chromosome != self.chromosome:
            return self._child(chromosome)
        else:
            return self


    def transpose_is(self, length):
        '''Produces a new chromosome via IS transposition'''
        # Since IS does not transpose to the root, it has no purpose
        # if the head length is less than 2.
        if self.head < 2:
            return self

        chromosome = list(self.chromosome)

        # Pick source and target genes
        source = random.choice(self._gene_starts)
        target = random.choice(self._gene_starts)

        # Extract a transposition sequence. Truncate if required.
        start = random.choice(xrange(self._gene_length))
        end   = start + length
        end   = self._gene_length if end > self._gene_length else end

        # Offset into target gene: in the head but not the root
        offset = random.choice(xrange(1, self.head))

        # Insert into the target gene's head
        chromosome[target:target+self.head] = (
            chromosome[target:target+offset]    +
            chromosome[source+start:source+end] +
            chromosome[target+offset:target+self.head]
        )[:self.head]

        if chromosome != self.chromosome:
            return self._child(chromosome)
        else:
            return self


    def transpose_ris(self, length):
        '''Produces a new chromosome via RIS transposition'''
        chromosome = list(self.chromosome)

        # Pick source and target genes
        source = random.choice(self._gene_starts)
        target = random.choice(self._gene_starts)

        # Extract a transposition sequence. Truncate if required.
        # For RIS the sequence must begin with a function.
        try:
            start = random.choice([
                i for i in xrange(source, source+self._gene_length)
                if callable(chromosome[i])
            ])
        except IndexError: # no functions!
            return self

        end, trunc = start+length, start+self._gene_length
        end = trunc if end > trunc else end

        # Insert into the target gene's head
        chromosome[target:target+self.head] = (
            chromosome[start:end] + chromosome[target:target+self.head]
        )[:self.head]

        if chromosome != self.chromosome:
            return self._child(chromosome)
        else:
            return self


    def transpose_gene(self):
        '''Produces a new chromosome via gene transposition'''
        s1 = random.choice(self._gene_starts)
        t1 = random.choice(self._gene_starts)
        if s1 == t1:
            return self

        # Switch these genes
        s2 = s1 + self._gene_length
        t2 = t1 + self._gene_length
        chromosome = list(self.chromosome)

        saved = chromosome[t1:t2]
        chromosome[t1:t2] = chromosome[s1:s2]
        chromosome[s1:s2] = saved

        if chromosome != self.chromosome:
            return self._child(chromosome)
        else:
            return self


    def crossover_one_point(self, other):
        '''
        Produces two children via one-point crossover
        @param other: second parent
        @return: child 1, child 2
        '''
        index = random.randint(0, len(self)-1)
        child1 = self.chromosome[:index] + other.chromosome[index:]
        child2 = other.chromosome[:index] + self.chromosome[index:]
        return self._child(child1), self._child(child2)


    def crossover_two_point(self, other):
        '''
        Produces two children via two-point crossover
        @param other: second parent
        @return: child 1, child 2
        '''
        if len(self) < 2:
            return self, other

        i1, i2 = random.sample(xrange(len(self)), 2)
        if i1 > i2:
            i1, i2 = i2, i1
        p1, p2 = self.chromosome, other.chromosome

        child1 = p1[:i1] + p2[i1:i2] + p1[i2:]
        child2 = p2[:i1] + p1[i1:i2] + p2[i2:]
        return self._child(child1), self._child(child2)


    def crossover_gene(self, other):
        '''
        Produces two children via full gene crossover
        @param other: second parent
        @return: child 1, child 2
        '''
        # Choose a random gene
        i1 = random.choice(self._gene_starts)
        i2 = i1 + self._gene_length
        p1, p2 = self.chromosome, other.chromosome

        child1 = p1[:i1] + p2[i1:i2] + p1[i2:]
        child2 = p2[:i1] + p1[i1:i2] + p2[i2:]
        return self._child(child1), self._child(child2)

