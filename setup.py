import os
from setuptools import find_packages, setup

# determining the directory containing setup.py
setup_path = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(setup_path, 'README.md'), encoding='utf-8') as f:
    readme = f.read()

setup(
    # package information
    name = 'nanogate',
    packages = find_packages(),
    version = '2.0.0-dev',
    description = 'Tool to verify DNA sequences after assembly',
    long_description = readme,
    license = 'MIT',
    url='git@github.com:stracquadaniolab/nanogate.git',
    keywords='',

    # author information
    author = 'Angelo Gaeta Giovanni Stracquadanio',
    author_email = 's1899269@ed.ac.uk',

    # installation info and requirements
    install_requires=[],
    setup_requires=[],

    # test info and requirements
    test_suite='tests',
    tests_require=[],

    # package deployment info
    include_package_data=True,
    zip_safe=False,

    # all tools have cli interface
    entry_points={
        'console_scripts': [
            'nanogate=nanogate.cli:main',
        ],
    },
)
