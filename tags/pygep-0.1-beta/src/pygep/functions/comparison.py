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


'''
Provides basic comparison non-terminals.  If each is true, it returns
the first value given.  If false, it returns the second.
'''


__all__ = 'equal', 'unequal', 'less', 'greater', \
          'less_or_equal', 'greater_or_equal'


@symbol('=')
def equal(x, y):
    '''Returns x if x == y else y'''
    return x if x == y else y


@symbol('U')
def unequal(x, y):
    '''Returns x if x != y else y'''
    return x if x != y else y


@symbol('<')
def less(x, y):
    '''Returns x if x < y else y'''
    return x if x < y else y


@symbol('>')
def greater(x, y):
    '''Returns x if x > y else y'''
    return x if x > y else y


@symbol('L')
def less_or_equal(x, y):
    '''Returns x if x <= y else y'''
    return x if x <= y else y


@symbol('G')
def greater_or_equal(x, y):
    '''Returns x if x >= y else y'''
    return x if x >= y else y
