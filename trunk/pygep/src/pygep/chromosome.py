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


class ChromosomeType(type):
    '''
    Metaclass for computing various information about a chromosomal
    type.  Sets the following attributes on a chromosome class:
        - arity:  maximum functional arity
        - tail:   length of the terminals tail
        - length: length of the head + length of the tail
        - head_symbols: symbols that can reside in the head
    '''
    def __new__(typ, *args, **kwds):
        t = type.__new__(typ, *args, **kwds)
        t.head_symbols = t.functions + t.terminals

        # Find the max arity
        try:
            t.arity = max([f.func_code.co_argcount for f in t.functions])
        except ValueError:
            t.arity = 0

        # And the tail length
        t.tail   = t.head * (t.arity - 1) + 1
        t.length = t.head + t.tail
        return t


class Chromosome(object):
    '''
    A Chromosome must provide these attributes:
        - functions: sequence of nonterminals symbols
        - terminals: sequence of terminal symbols
        - head:      length of the head
    '''
    __metaclass__ = ChromosomeType

    functions = ()
    terminals = ()
    head = 0


    @classmethod
    def generate(cls):
        '''Generates a random chromosome'''
        for i in xrange(cls.length):
            genes = [choice(cls.head_symbols) for _ in xrange(cls.head)] + \
                    [choice(cls.terminals) for _ in xrange(cls.tail)]
            return cls(genes)
        

    def __init__(self, genes):
        self.genes = genes


    def __len__(self):
        return self.length 


    @cache
    def __repr__(self):
        s = ''
        for gene in self.genes:
            # Differentiate between functions and terminals
            name = gene
            if callable(gene):
                name = gene.__name__

            # It is discouraged, but names may be longer than 1
            if len(name) > 1:
                s += '{%s}' % name
            else:
                s += name

        return s

    
    def fitness(self):
        '''Abstract base method for fitness determination'''
        raise NotImplementedError('Must override Chromosome.fitness')

