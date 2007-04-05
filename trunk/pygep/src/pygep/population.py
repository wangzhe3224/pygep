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
from itertools import chain
from random import choice
import string


class Population(object):
    def __init__(self, cls, size):
        '''
        Generates a population of some chromsome class
        @param cls:  Chromosome type
        @param size: population size
        '''
        self.size      = size

        # Start an initial population
        population = [None] * self.size
        for i in xrange(self.size):
            population[i] = cls.generate()
        
        self.population = population

        # Header for display purposes
        self.header = string.digits * (cls.length / len(string.digits)) + \
                      string.digits[:(cls.length % len(string.digits))]


    def __repr__(self):
        return '\n'.join([str(i) for i in [self.header] + self.population])
        
