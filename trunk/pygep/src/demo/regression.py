from pygep.chromosome import Chromosome, symbol
from pygep.population import Population
import math, random


# The function we are trying to find
def target(x):
    # 3x^2 + 2x + 1
    return 3*(x.a**2) + (2*x.a) + 1


# Generate a random sample to test against
class Data(object):
    def __init__(self, a):
        self.a = float(a)


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
    terminals = 'a', 1, 2

    selection_range = 100.0
    range_low, range_high = -10.0, 10.0
    range_size = range_high - range_low
    sample = [Data(range_low + (random.random() * range_size)) for _ in xrange(10)]
    
    def _fitness(self):
        try:
            return int(max(sum(self.selection_range - abs(self.evaluate(x)-target(x))
                               for x in self.sample), 0.0))
        except:
            return 0

    def _solved(self):
        return self.fitness >= (self.selection_range * len(self.sample))


if __name__ == '__main__':
    # Search for a solution
    p = Population(Regression, 20, 6, 3, lambda *args: sum(args))
    print p

    for _ in xrange(50):
        if p.best.solved:
            break
        p.cycle()
        print
        print p

    if p.best.solved:
        print
        print 'SOLVED:', p.best
