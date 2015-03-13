#!/usr/bin/python
"""
Build the csolve module for seismic-2d.

Author: Nimar Arora (feel free to modify, use, or redistribute)
Copyright 2015.
"""
from distutils.core import setup, Extension

module1 = Extension('csolve', sources = ['csolve.c'])

setup (name = 'C-Solver',
       version = '1.0',
       description = 'C-based solver for seismic-2d',
       ext_modules = [module1])
