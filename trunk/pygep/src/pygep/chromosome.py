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
from pygep.util import cache
from random import choice


class MetaChromosome(type):
    '''
    Metaclass for computing various information about a chromosomal
    type.  Sets the following attributes on a chromosome class:
        - arity:   maximum functional arity
        - symbols: symbols that can reside in the head
    '''
    def __new__(typ, *args, **kwds):
        t = type.__new__(typ, *args, **kwds)
        t.symbols = t.functions + t.terminals

        # Find the max arity
        try:
            t.arity = max([f.func_code.co_argcount for f in t.functions])
        except ValueError:
            t.arity = 0

        return t


class Chromosome(object):
    '''
    A Chromosome must provide these attributes:
        - functions: sequence of nonterminals symbols
        - terminals: sequence of terminal symbols
    '''
    __metaclass__ = MetaChromosome

    functions = ()
    terminals = ()
    head = tail = length = 0


    @classmethod
    def generate(cls, n, head):
        '''
        Returns a generator of n random GEP chromosomes
        @param n:    number of chromosomes
        @param head: head length
        '''
        tail = head * (cls.arity - 1) + 1

        for _ in xrange(n):
            yield cls([choice(cls.symbols)   for _ in xrange(head)] + \
                      [choice(cls.terminals) for _ in xrange(tail)])


    def __init__(self, genes):
        '''
        Instantiates a chromsome instance and analyzes it for evaluation.
        Sets the self.coding index to the last gene in the coding region.
        '''
        self.genes = genes

        # Starts at the first gene and determines how many args it requires.
        # Then moves forward that many args and sums their required args.
        # Continues this process until there are no more required args.  The
        # resulting index will be one gene past the coding region.
        index, args = 0, 1
        while args:
            next_args = 0
            for _ in xrange(args):
                if callable(genes[index]):
                    next_args += genes[index].func_code.co_argcount
                index += 1

            args = next_args

        self.coding = index - 1

        # This attribute is a prepopulated list to make evaluating
        # the chromosome against data instances more efficient.
        self._eval = list(genes)


    def __len__(self):
        return len(self.genes)


    @cache
    def __repr__(self):
        s = ''
        for gene in self.genes:
            # Differentiate between functions and terminals
            name = gene.__name__ if callable(gene) else gene

            # It is discouraged, but names may be longer than 1
            if len(name) > 1:
                s += '{%s}' % name
            else:
                s += name

        return s


    def evaluate(self, obj):
        '''
        Evaluates a given GEP chromosome against some instance.  The
        terminals in the chromosome are assumed to be attributes on
        the object instance provided.
        @param obj: an object instance with terminal attributes set
        @return: result of evaluating the chromosome
        '''
        # Start by clearing out our eval list
        for i, gene in enumerate(self.genes):
            self._eval[i] = None if callable(gene) else getattr(obj, gene)

        # Evaluate the chromosome against obj in reverse
        index = self.coding + 1
        for i in reversed(xrange(index)):
            gene = self.genes[i]

            if callable(gene):
                num  = gene.func_code.co_argcount
                args = self._eval[index-num:index]

                # Replace the operation in self._eval with its return value
                self._eval[i] = gene(*args)
                index -= num

        return self._eval[i]


    def fitness(self):
        '''Abstract base method for fitness determination'''
        raise NotImplementedError('Must override Chromosome.fitness')

