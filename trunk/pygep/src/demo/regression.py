#!/usr/bin/env python2.5
from pygep.functions.arithmetic import *
from pygep.functions.linkers import sum_linker
from pygep import *
import random


# Data points to test against and target function
class DataPoint(object):
    SAMPLE = []
    SAMPLE_SIZE = 10
    RANGE_LOW, RANGE_HIGH = -10.0, 10.0
    RANGE_SIZE = RANGE_HIGH - RANGE_LOW

    def __init__(self, x):
        self.x = float(x)

        # The function we are trying to find
        # f(x) = 4*(x**4) + 3x^2 + 2x + 1
        self.y = 4*(x**3) + 3*(x**2) + (2*x) + 1

    @staticmethod
    def populate():
        # Creates a random sample of data points
        DataPoint.SAMPLE = []
        for _ in xrange(DataPoint.SAMPLE_SIZE):
            x = DataPoint.RANGE_LOW + (random.random() * DataPoint.RANGE_SIZE)
            DataPoint.SAMPLE.append(DataPoint(x))


# The chromsomes: fitness is accuracy over the sample
class Regression(Chromosome):
    SELECTION_RANGE = 1000.0
    functions = multiply, add, subtract, divide
    terminals = 'x', 1, 2

    def _fitness(self):
        try:
            total = sum(self.SELECTION_RANGE - abs(self.evaluate(x)-x.y)
                        for x in DataPoint.SAMPLE)
            return int(max(total, 0.0))
        except:
            return 0

    def _solved(self):
        return self.fitness >= (self.SELECTION_RANGE * len(DataPoint.SAMPLE))


if __name__ == '__main__':
    DataPoint.populate()

    # Search for a solution
    p = Population(Regression, 20, 8, 4, sum_linker)
    print p

    for _ in xrange(100):
        if p.best.solved:
            break
        p.cycle()
        print
        print p

    if p.best.solved:
        print
        print 'SOLVED:', p.best

