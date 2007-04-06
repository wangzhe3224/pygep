from pygep.chromosome import Chromosome, symbol
from pygep.population import Population
from pygep.util import cache
import math, random, sys

# The functin we are trying to find
def target(x):
    return 4 * x

# Generate a random sample to test against
class Data(object):
    def __init__(self, a):
        self.a = a

# The functions we use in our chromosome
@symbol('*')
def multiply(x, y):
    return x * y

@symbol('+')
def add(x, y):
    return x + y

@symbol('-')
def subtract(x, y):
    return x - y

@symbol('/')
def divide(x, y):
    try:
        return x / y
    except ZeroDivisionError:
        return sys.maxint

def Q(x):
    try:
        return math.sqrt(x)
    except ValueError:
        return 0.0

# The chromsomes: fitness is accuracy over the sample
class Regression(Chromosome):
    functions = multiply, add, subtract, divide, Q
    terminals = 'a',
    sample = [Data(float(random.randint(1, 10))) for _ in xrange(25)]

    @cache
    def fitness(self):
        diff = 0
        for data in self.sample:
            diff += float(abs(target(data.a) - self.evaluate(data))) / data.a

        return -diff

p = Population(Regression, 30, 10)
for i in p:
    print i, i.fitness()
