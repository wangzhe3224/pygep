from tests.base import Computation
import unittest


class ChromosomeTest(unittest.TestCase):
    '''Tests basic components and functionality of a chromosomes'''
    head  = 5
    genes = 3
    symbols = Computation.functions + Computation.terminals
    tail_length = head * (Computation.arity-1) + 1
    gene_length = head + tail_length


    def setUp(self):
        self.chromosome = Computation.generate(self.head, self.genes).next()
        

    def testArity(self):
        # Max arity of functions
        self.assertEqual(Computation.ARITY, Computation.arity)


    def testLength(self):
        # Head / tail / gene length
        self.assertEqual(self.gene_length * self.genes, len(self.chromosome))


    def testStructure(self):
        self.assertEqual(Computation.symbols, self.symbols)

        for gene in self.chromosome.genes:
            # Head can contain both terminals and non-terminals
            for allele in gene[:self.head]:
                self.assertTrue(allele in self.symbols)

            # Tail only contains terminals
            for allele in gene[self.head:]:
                self.assertTrue(allele in Computation.terminals)
        

if __name__ == '__main__':
    unittest.main()
