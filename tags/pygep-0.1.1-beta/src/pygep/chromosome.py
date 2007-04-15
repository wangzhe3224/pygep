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

    functions = ()
    terminals = ()
    head = tail = length = 0


    @classmethod
    def generate(cls, head, genes=1, linker=lambda x: x):
        '''
        Returns a generator of random GEP chromosomes
        @param genes:  number of genes
        @param head:   head length
        @param linker: linking function
        '''
        tail = head * (cls.arity - 1) + 1

        while True:
            chromosome = []
            for _ in xrange(genes):
                chromosome.extend(
                    [random.choice(cls.symbols)   for _ in xrange(head)] + \
                    [random.choice(cls.terminals) for _ in xrange(tail)]
                )
            yield cls(chromosome, head, genes, linker)


    def __init__(self, chromosome, head, genes, linker):
        '''
        Instantiates a chromsome instance and analyzes it for evaluation.
        Sets the self.coding tuple to the last genes in the coding regions
        and various other internal data for the chromosome.  Note that it
        is generally unwise to instantiate chromosomes manually.  It is
        much more common to create them via calls to the static method
        Chromosome.generate(...).

        @param chromosome: combined list of all genes
        @param head:       length (not index) of the gene heads
        @param genes:      number of genes in the chromosome
        @param linker:     linker function for gene evaluation
        '''
        self.chromosome = chromosome
        self.head       = head
        self.genes      = genes
        self.linker     = linker

        # Useful information: the length of each gene
        self._gene_length = len(self) / self.genes

        # Unique number of the organism
        self.__id = type(self).__next_id
        type(self).__next_id += 1

        # Determine coding indexes for each gene:
        self.coding = ()
        for start in self._gene_starts:
            # How to find the length of a single coding region:
            #
            # Start at the first gene and determine how many args it
            # requires. Then move forward that many args and sum their
            # required args. Continue this until there are no more
            # required args. The resulting index will be one gene past
            # the coding region for the current gene
            index, args = start, 1
            while args:
                next_args = 0
                for _ in xrange(args):
                    if callable(chromosome[index]):
                        next_args += chromosome[index].func_code.co_argcount
                    index += 1

                args = next_args

            self.coding += index - 1,

        # This attribute is a prepopulated list to make evaluating
        # the chromosome against data instances more efficient.
        self._eval = list(chromosome)


    def __len__(self):
        return len(self.chromosome)


    @cache
    def __repr__(self):
        s = ''
        for gene in self.chromosome:
            # Differentiate between functions and terminals
            try:
                name = gene.symbol
            except AttributeError:
                try:
                    name = gene.__name__
                except AttributeError:
                    name = str(gene)

            # Truncate any name longer than 1 char
            s += name[:1] or '?'

        return s


    _gene_starts = property(
        lambda self: xrange(0, len(self), self._gene_length),
        doc='Indexes where each individual gene commences'
    )


    def _child(self, chromosome):
        '''Returns a child chromosome of self'''
        return type(self)(chromosome, self.head, self.genes, self.linker)


    # Unique ID of the organism
    id = property(lambda self: self.__id, doc='Organism #')


    def evaluate(self, obj):
        '''
        Evaluates a given GEP chromosome against some instance.  The
        terminals in the chromosome are assumed to be attributes on
        the object instance provided (unless they are numeric constants).

        @param obj: an object instance with terminal attributes set
        @return:    result of evaluating the chromosome
        '''
        # Start by clearing out our eval list
        for i, gene in enumerate(self.chromosome):
            self._eval[i] = None
            if callable(gene):
                self._eval[i] = None
            elif not isinstance(gene, str):
                # could be a number
                self._eval[i] = gene
            else:
                self._eval[i] = getattr(obj, gene)


        # Evaluate each chromosome gene against obj in reverse
        for coding, start in izip(self.coding, self._gene_starts):
            index = coding + 1
            for i in reversed(xrange(start, index)):
                gene = self.chromosome[i]

                if callable(gene):
                    num  = gene.func_code.co_argcount
                    args = self._eval[index-num:index]

                    # Replace the operation in self._eval with its return val
                    self._eval[i] = gene(*args)
                    index -= num

        return self.linker(*self._eval[0:len(self):self._gene_length])


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
        chromosome = list(self.chromosome)

        # Traverse the chromosome gene by gene
        for start in self._gene_starts:
            # Then locus by locus
            for i in xrange(self._gene_length):
                # Do we mutate this locus?
                if random.random() < rate:
                    # Mutation within the tail can only use terminals
                    if i >= self.head:
                        gene = random.choice(self.terminals)
                    else:
                        gene = random.choice(self.symbols)

                    chromosome[start+i] = gene

        # Only create a new chromosome if it is actually different
        if chromosome != self.chromosome:
            return self._child(chromosome)
        else:
            return self


    def invert(self):
        '''Produces a new chromosome via head inversion'''
        chromosome = list(self.chromosome)

        # Choose a random gene and two points within the head
        gene = random.choice(self._gene_starts)
        indexes = random.sample(xrange(self.head), 2)
        start, stop = min(indexes), max(indexes)

        # Create the new chromosome
        chromosome[start:stop] = reversed(chromosome[start:stop])
        if chromosome != self.chromosome:
            return self._child(chromosome)
        else:
            return self


    def transpose_is(self, length):
        '''Produces a new chromosome via IS transposition'''
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
            chromosome[target+offset:]
        )[:self.head]

        if chromosome != self.chromosome:
            return self._child(chromosome)
        else:
            return self


    def transpose_is(self, length):
        '''Produces a new chromosome via IS transposition'''
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
        indexes = random.sample(xrange(len(self)), 2)
        i1, i2 = min(indexes), max(indexes)
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

