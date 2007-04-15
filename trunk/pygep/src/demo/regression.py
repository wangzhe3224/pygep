#!/usr/bin/env python2.5
from pygep.functions.arithmetic import *
from pygep.functions.linkers import sum_linker
from pygep import *
import random


# The function we are trying to find
def target(x):
    # 3x^2 + 2x + 1
    return 3*(x.a**2) + (2*x.a) + 1


# Data points to test against
SELECTION_RANGE = 100.0
RANGE_LOW, RANGE_HIGH = -10.0, 10.0
RANGE_SIZE = RANGE_HIGH - RANGE_LOW
class Data(object):
    def __init__(self, a):
        self.a = float(a)

NUM = 10
SAMPLE = [Data(RANGE_LOW + (random.random() * RANGE_SIZE)) for _ in xrange(NUM)]


# The chromsomes: fitness is accuracy over the sample
class Regression(Chromosome):
    functions = multiply, add, subtract, divide
    terminals = 'a', 1, 2

    selection_range = 100.0

    def _fitness(self):
        try:
            total = sum(SELECTION_RANGE - abs(self.evaluate(x)-target(x)) 
                        for x in SAMPLE)
            return int(max(total, 0.0))
        except:
            return 0

    def _solved(self):
        return self.fitness >= (SELECTION_RANGE * NUM)


if __name__ == '__main__':
    # Search for a solution
    p = Population(Regression, 20, 6, 3, sum_linker)
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

