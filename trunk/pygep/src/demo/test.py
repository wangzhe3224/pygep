from pygep.chromosome import Chromosome, symbol
from pygep.population import Population
from pygep.util import cache
import math, random, sys

# The functin we are trying to find
def target(x, y):
    return 4 * x + y

# Generate a random sample to test against
class Data(object):
    def __init__(self, a):
        self.a = float(a)
        self.b = 1.0

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

# The chromsomes: fitness is accuracy over the sample
class Regression(Chromosome):
    functions = multiply, add, subtract, divide
    terminals = 'a', 'b'
    sample = [Data(float(random.randint(1, 10))) for _ in xrange(25)]

    @cache
    def fitness(self):
        good = 0
        for data in self.sample:
            desired = float(target(data.a, 1))
            closeness = abs((self.evaluate(data)-desired) / desired)
            if closeness < .1:
                good += 3
            elif abs(closeness) < .25:
                good += 1

        return good

#p = Population(Regression, 30, 7)
p = Population(Regression, 30, 7, 5, lambda *args: sum(args))
d = Data(2)
for i in p:
    print i, i.evaluate(d)

