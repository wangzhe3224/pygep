from tests.base import Computation
import unittest


class InversionTest(unittest.TestCase):
    def testShortInversion(self):
        # 5: Makes sure inversion works with head length=1
        i = Computation.generate(1, 1).next()
        self.assertEqual(i, i.invert())


if __name__ == '__main__':
    unittest.main()
