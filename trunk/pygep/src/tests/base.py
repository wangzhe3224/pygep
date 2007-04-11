from pygep.chromosome import Chromosome, symbol


'''Base types for use in other tests'''


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

FUNCTIONS = multiply, add, subtract, divide


class Computation(Chromosome):
    ARITY = 2
    functions = FUNCTIONS
    terminals = 'a', 1, 2

