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
from operator import attrgetter
from pygep.util import stats
import random, string


class Population(object):
    '''
    A Population instance has the following default configuration:
        - selection_pressure: sigma-scaled pressure (1.2 to 2)
        - mutation_rate: probability that each organism mutates
        - crossover_one_point_rate: probability of 1-point crossover
        - crossover_two_point_rate: probability of 2-point crossover
        - crossover_gene_rate: probability of full gene crossover
    '''
    selection_pressure       = 1.2
    mutation_rate            = 0.005
    crossover_one_point_rate = 0.3
    crossover_two_point_rate = 0.3
    crossover_gene_rate      = 0.1


    def __init__(self, cls, size, head, genes=1, linker=lambda x: x):
        '''
        Generates a population of some chromsome class
        @param cls:    Chromosome type
        @param size:   population size
        @param head:   chromosome head length
        @param genes:  number of genes
        @param linker: multigenic results linker function
        '''
        self.size   = size
        self.head   = head
        self.genes  = genes
        self.linker = linker

        self.__age = 1

        # Start an initial population
        self.population = [i for _, i in izip(xrange(size),
                           cls.generate(head, genes, linker))]
        self._next_pop = [None] * size # placeholder for next generation
        self._scaled   = [None] * size # fitness scaling

        # Header for display purposes
        try:
            l = len(self.population[0])
            self.header = string.digits * (l / len(string.digits)) + \
                          string.digits[:(l % len(string.digits))]
            self.header += '\n' + '-' * len(self.header)
        except IndexError:
            raise ValueError('Empty populations are meaningless!')


    def __repr__(self):
        return '\n'.join([str(i) for i in [self.header] + self.population])


    def __len__(self):
        return self.size


    def __iter__(self):
        return iter(self.population)


    age  = property(lambda self: self.__age, doc='Generation number')
    best = property(
        lambda self: max(self.population, key=attrgetter('fitness')),
        doc='The best Chromosome of the current generation'
    )


    def solve(self, generations):
        '''Cycles a number of generations. Stops if chrom.solved()'''
        pass


    def cycle(self):
        '''Selects and recombines the next generation'''
        # Copy the best individual via simple elitism
        self._next_pop[0] = self.best

        # Fill in the rest through fitness scaling. First compute the
        # sigma-scaled fitness proportionate value for each chromosome.
        mean, stdev = stats.fitness_stats(self)
        for i, c in enumerate(self.population):
            scaled = c.fitness - (mean - (self.selection_pressure * mean))
            self._scaled[i] = max(scaled, 0.0)
        scaling = sum(self._scaled)

        # Then generate n-1 spins of the roulette wheel
        select = [random.random() * scaling for _ in xrange(self.size-1)]
        select.sort()

        # Moved through the current population, calculating a window
        # for each scaled fitness.  However many of the sorted roulette
        # spins are within that window yields the number of copies.
        window = self._scaled[0]
        source, target = 0, 1
        for scaled in select:
            # Move to the element that covers our range
            while window < scaled:
                source += 1
                try:
                    window += self._scaled[source]
                except IndexError: # possible floating-point errors
                    source = self.size - 1
                    break

            # Copy this element to the next generation
            self._next_pop[target] = self.population[source]
            target += 1

        # Recombination section:
        # First try and mutate each individual
        if self.mutation_rate:
            for i, c in enumerate(self._next_pop):
                self._next_pop[i] = c.mutate(self.mutation_rate)

        # Then try one|two-point and gene crossover
        if random.random() < self.crossover_one_point_rate:
           i1, i2 = random.sample(xrange(self.size), 2)
           p1, p2 = self._next_pop[i1], self._next_pop[i2]
           self._next_pop[i1], self._next_pop[i2] = p1.crossover_one_point(p2)

        if random.random() < self.crossover_two_point_rate:
           i1, i2 = random.sample(xrange(self.size), 2)
           p1, p2 = self._next_pop[i1], self._next_pop[i2]
           self._next_pop[i1], self._next_pop[i2] = p1.crossover_two_point(p2)

        if random.random() < self.crossover_gene_rate:
           i1, i2 = random.sample(xrange(self.size), 2)
           p1, p2 = self._next_pop[i1], self._next_pop[i2]
           self._next_pop[i1], self._next_pop[i2] = p1.crossover_gene(p2)

        # Switch to the next generation and increment age
        self._next_pop, self.population = self.population, self._next_pop
        self.__age += 1

