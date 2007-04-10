from pygep.chromosome import Chromosome, symbol
from pygep.population import Population
import math, random


# The function we are trying to find
def target(x):
    #return 3 * (x ** 2) + (2 * x)

    # 3a^2 + 2b^3 + 1
    return 3 * (x.a ** 2) + 2 * (x.b ** 3) + 1


# Generate a random sample to test against
class Data(object):
    def __init__(self, a, b):
        self.a = float(a)
        self.b = float(b)


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
    return x / y



# The chromsomes: fitness is accuracy over the sample
class Regression(Chromosome):
    functions = multiply, add, subtract, divide
    terminals = 'a', 'b', 1, 2

    range_low, range_high = -10.0, 10.0
    range_size = range_high - range_low
    sample = [Data(range_low + (random.random() * range_size),
                   range_low + (random.random() * range_size)) for _ in xrange(10)]

    def _fitness(self):
        total = 0
        for x in self.sample:
            try:
                diff = 1 - (self.evaluate(x) / float(target(x)))
                if diff < 0.001:
                    total += 20
                if diff < 0.01:
                    total += 10
                elif diff < 0.05:
                    total += 5
                elif diff < 0.10:
                    total += 1
            except: # Let nature sort out semantic errors
                pass

        return total

    def _solved(self):
        return self.fitness >= (20 * len(self.sample))


if __name__ == '__main__':
    # Search for a solution
    p = Population(Regression, 20, 6, 3, lambda *args: sum(args))
    print p

    for _ in xrange(1000):
        if p.best.solved:
            break
        p.cycle()
        print
        print p

    if p.best.solved:
        print
        print 'SOLVED:', p.best
