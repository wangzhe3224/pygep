from tests.base import Computation
import unittest


def z(x, y):
    pass


class MutationTest(unittest.TestCase):
    '''Verifies effects of mutation per allele'''
    def setUp(self):
        self.chromosome = Computation([z,z,'y',3,4], 2, 1, lambda x: x)


    def testAllMutate(self):
        newchr = self.chromosome.mutate(1.1)
        for allele in newchr.chromosome:
            self.assertTrue(allele not in [z,z,'y',3,4])
        self.assertNotEqual(self.chromosome.id, newchr.id)


    def testNoMutate(self):
        newchr = self.chromosome.mutate(0)
        self.assertEqual(newchr.chromosome, self.chromosome.chromosome)
        self.assertEqual(self.chromosome.id, newchr.id)


if __name__ == '__main__':
    unittest.main()

