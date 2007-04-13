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


'''
Provides linkers for combining multigenic chromosomes.  In general
PyGEP linkers should accept and process any number of arguments,
thus the use of *args.
'''


def sum_linker(*args):
    '''Returns the sum of all sub-ETs'''
    return sum(args)


def or_linker(*args):
    '''Returns the OR of all given args'''
    for x in args:
        if x:
            return True
        return False
