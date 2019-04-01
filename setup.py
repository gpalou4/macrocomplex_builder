#!/usr/bin/env python

from setuptools import setup

setup(name='Macrocomplex Builder',
	version='1.0',
	description='This program is able to reconstruct biological macrocomplexes of protein-protein interactions as well as protein-DNA/RNA interactions given a set of binary interactions and the desired number of chains of the target complex.',
	author='Guillermo Palou Marquez and Javier Sanchez Utges',
	author_email='guillepalou4@gmail.com',
	long_description=open('README.md').read(),
	install_requires=['biopython >= 1.73.0','argparse >= 1.1.0'],
	license='LICENSE.txt',
	url='https://github.com/gpalou4/macrocomplex_builder',
	scripts=['macrocomplex_builder.py','macrocomplex_functions.py'])