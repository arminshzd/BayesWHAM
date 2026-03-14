#!/usr/bin/env python3
"""
Setup script for BayesWHAM package
"""

from setuptools import setup, find_packages
import os

# Read the README file for long description
def read_readme():
    readme_path = os.path.join(os.path.dirname(__file__), 'README.md')
    if os.path.exists(readme_path):
        with open(readme_path, 'r', encoding='utf-8') as f:
            return f.read()
    return ""

# Read version from package
def get_version():
    version = {}
    with open(os.path.join('bayeswham', '__init__.py'), 'r') as f:
        for line in f:
            if line.startswith('__version__'):
                exec(line, version)
                return version['__version__']
    return "1.1.0"

setup(
    name='bayeswham',
    version=get_version(),
    description='Bayesian Weighted Histogram Analysis Method for molecular free energy estimation',
    long_description=read_readme(),
    long_description_content_type='text/markdown',
    author='Andrew L Ferguson',
    author_email='',
    url='https://github.com/FergusonAJ/BayesWHAM',
    license='MIT',

    packages=find_packages(),

    # Include package data
    package_data={
        'bayeswham.config': ['*.yaml'],
    },

    # Dependencies
    install_requires=[
        'numpy>=1.19.0',
        'scipy>=1.5.0',
        'pyyaml>=5.3.0',
    ],

    # Optional dependencies
    extras_require={
        'plotting': ['matplotlib>=3.3.0'],
        'dev': [
            'pytest>=6.0.0',
            'pytest-cov>=2.10.0',
        ],
    },

    # Console scripts / command-line entry points
    entry_points={
        'console_scripts': [
            'bayeswham=bayeswham.core.bayeswham:main',
            'bayesreweight=bayeswham.core.bayesreweight:main',
            'bayeswham-plot=bayeswham.plotting.bayeswham_plotter:main',
            'bayesreweight-plot=bayeswham.plotting.bayesreweight_plotter:main',
        ],
    },

    # Classifiers
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
    ],

    # Python version requirement
    python_requires='>=3.7',

    # Keywords
    keywords='molecular dynamics, free energy, umbrella sampling, WHAM, Bayesian statistics',

    # Project URLs
    project_urls={
        'Documentation': 'https://github.com/FergusonAJ/BayesWHAM/blob/master/README.md',
        'Source': 'https://github.com/FergusonAJ/BayesWHAM',
        'Tracker': 'https://github.com/FergusonAJ/BayesWHAM/issues',
    },
)
