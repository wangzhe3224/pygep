from pkg_resources import require
from setuptools import setup, find_packages

require('python>=2.5')

setup(
    name        ='PyGEP',
    version     ='0.1',
    package_dir = {'': 'src'},
    packages    = find_packages('src', exclude=['tests', 'tests.*']),
    zip_safe    = True,
    test_suite  = 'tests'
)
