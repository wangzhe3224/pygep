#!/usr/bin/env python2.5
from pygep.functions.arithmetic import *
from pygep.functions.linkers import sum_linker
from pygep import *
import random


# Data points to test against and target function
SAMPLE = []
class DataPoint(object):
    RANGE_LOW, RANGE_HIGH = -10.0, 10.0
    RANGE_SIZE = RANGE_HIGH - RANGE_LOW

    def __init__(self, x):
        self.x = float(x)

        # The function we are trying to find
        # f(x) = 3x^2 + 2x + 1
        self.y = 3*(x**2) + (2*x) + 1



# The chromsomes: fitness is accuracy over the sample
SELECTION_RANGE = 100.0
class Regression(Chromosome):
    functions = multiply, add, subtract, divide
    terminals = 'x', 1, 2

    def _fitness(self):
        try:
            total = sum(SELECTION_RANGE - abs(self.evaluate(x)-x.y) 
                        for x in SAMPLE)
            return int(max(total, 0.0))
        except:
            return 0

    def _solved(self):
        return self.fitness >= (SELECTION_RANGE * NUM)


if __name__ == '__main__':
    # Create a random sample of data points
    NUM = 10
    for _ in xrange(NUM):
        x = DataPoint.RANGE_LOW + (random.random() * DataPoint.RANGE_SIZE)
        SAMPLE.append(DataPoint(x))

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

