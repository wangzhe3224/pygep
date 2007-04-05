from pygep.chromosome import Chromosome
from pygep.population import Population
import sys

def X(x, y, z):
    return (x + y + z) / 3.

def M(x, y):
    return x * y

def A(x, y):
    return x + y

def S(x, y):
    return x - y

def D(x, y):
    try:
        return x / y
    except ZeroDivisionError:
        return sys.maxint


class Test(Chromosome):
    functions = M, A, S, D, X
    terminals = 'x', 'y'
    head = 20 


p = Population(Test, 30)
print p

