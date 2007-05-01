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
import functools


def cache(func):
    '''
    Decorator for caching the return value of an instance level method in self.
    Assumes that there are no arguments passed to the method (not memoization).
    The return value is cached on self._{method}_cache where {method} is the
    name of the method.
    
        @cache
        def _get_something(self):
            ...
            return 'something'
    '''
    cache_name = '_%s_cache' % func.func_name

    @functools.wraps(func)
    def wrapper(self):
        try:
            return getattr(self, cache_name)
            
        except AttributeError:
            setattr(self, cache_name, func(self))
            return getattr(self, cache_name)

    return wrapper

