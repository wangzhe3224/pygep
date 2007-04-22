from pygep.util import cache
import unittest


class Foo(object):
    '''A dummy class for testing the cache'''
    x = 0

    @cache
    def bar(self):
        Foo.x += 1
        return 'baz'


class CacheTest(unittest.TestCase):
    '''Verifies that caching works on an instance level'''
    def testCache(self):
        f = Foo()
        self.assertEqual('baz', f.bar())
        self.assertEqual(1, f.x)
        self.assertEqual('baz', f.bar())
        self.assertEqual(1, f.x)

        f = Foo()
        self.assertEqual('baz', f.bar())
        self.assertEqual(2, f.x)
        self.assertEqual('baz', f.bar())
        self.assertEqual(2, f.x)


if __name__ == '__main__':
    unittest.main()

