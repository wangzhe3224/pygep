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
from pygep.chromosome import symbol
import math


'''
Provides basic arithmetic non-terminals for use in Chromosomes.
Any semantic exceptions (ZeroDivisionError, etc.) are passed up
the call chain to the user.  Typically one should catch any
exceptions when calling chromosome.evaluate() and set the fitness
of nonviable organisms to 0.
'''


__all__ = 'add', 'subtract', 'multiply', 'divide', 'power', 'root'


@symbol('*')
def multiply(x, y):
    '''Returns x * y'''
    return x * y


@symbol('+')
def add(x, y):
    '''Returns x + y'''
    return x + y


@symbol('-')
def subtract(x, y):
    '''Returns x - y'''
    return x - y


@symbol('/')
def divide(x, y):
    '''
    Returns x / y
    @raises ZeroDivisionError
    '''
    return x / y


@symbol('^')
def power(x, y):
    '''Returns x ** y'''
    return x ** y


@symbol('Q')
def root(x):
    '''
    Returns square root of x
    @raises ValueError
    '''
    return math.sqrt(x)
