from pygep.chromosome import Chromosome
from pygep.population import Population
import math, sys

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

def Q(x):
    try:
        return math.sqrt(x)
    except ValueError:
        return 0.0

class Test(Chromosome):
    functions = M, A, S, D, Q
    terminals = 'a', 'b'

class Data(object):
    a = 1.
    b = 2.

p = Population(Test, 30, 20)
d = Data()
for i in p:
    print i, i.coding, i.evaluate(d)

