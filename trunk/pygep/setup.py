from setuptools import setup, find_packages

setup(
    name        ='PyGEP',
    version     ='0.1',
    package_dir = {'': 'src'},
    packages    = find_packages('src', exclude=['tests.*']),
    scripts     = [],
    zip_safe    = True
)
