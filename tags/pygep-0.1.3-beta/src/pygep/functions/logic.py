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
Provides basic logical operators: and, or, not, if.  The former 
three will return boolean values, whereas if returns one of the
values passed in, so be careful mixing these with other operators.

Common logic non-terminal functions:
    - (&) and_op
    - (|) or_op
    - (!) not_op
    - (I) if_op
'''


__all__ = 'and_op', 'or_op', 'not_op', 'if_op'


@symbol('&')
def and_op(x, y):
    '''Returns x and y'''
    return x and y


@symbol('|')
def or_op(x, y):
    '''Returns x or y'''
    return x or y


@symbol('!')
def not_op(x):
    '''Returns not x'''
    return not x


@symbol('I')
def if_op(x, y, z):
    '''Returns y if x else z'''
    return y if x else z

