"""
BayesWHAM - Bayesian Weighted Histogram Analysis Method

A Bayesian implementation of WHAM for estimating multidimensional molecular
free energy surfaces from umbrella sampling simulations.

Copyright (c) 2016 Andrew L Ferguson
Updated for Python 3 compatibility - 2025
"""

__version__ = "1.1.0"
__author__ = "Andrew L Ferguson"
__license__ = "MIT"

from .core import BayesWHAM, BayesReweight

__all__ = ['BayesWHAM', 'BayesReweight', '__version__']
