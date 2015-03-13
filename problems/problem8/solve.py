#!/usr/bin/python
"""
Solves the seismic-2d episodes using a very simple algorithm.

./solve.py data/physics.data data/test.blind data/test.solution

Note: this version of the solution reads the underlying physics, which
is, in fact, cheating. A future version will learn the physics directly
from training data.

Author: Nimar Arora (feel free to modify, use, or redistribute)
"""
from __future__ import print_function, division

import sys

from scipy.stats import laplace
from scipy.stats import gamma, norm, invgamma, poisson, uniform, expon,\
     bernoulli, laplace, cauchy
from numpy import array, sqrt, log, exp, pi, arcsin, degrees, linspace, seterr,\
     inf
seterr(all='raise')

from util import *

try:
  from csolve import *
except ImportError:
  print("Warning: csolve module not found. Please build it with:\n"
        "  python setup.py build_ext --inplace\n"
        "Reverting to python-based solver module, which is very slow!!")
  from pysolve import *

def main():
  if len(sys.argv) != 4:
    print("Usage: solve.py physics.data test.blind test.solution",
          file = sys.stderr)
    sys.exit(1)

  phys_file_name, blind_file_name, solution_file_name = sys.argv[1:]
  
  physics = read_namedtuple("Physics", phys_file_name)
  
  testdata = read_episodes(blind_file_name)

  # we will write out each episode as we solve it so that it can
  # be concurrently analyzed
  fp = open(solution_file_name, "w")
  
  print("Episodes:", file=fp, end="\n\n")
  
  for episode in testdata:
    
    events, assocs = solve_episode(physics, STATIONS, episode.detections)
    
    write_single_episode(Episode(events, [], assocs), fp)
    
    fp.flush()
    
  fp.close()

if __name__ == "__main__":
  try:
    main()
  except SystemExit:
    raise
  except:
    import pdb, traceback, sys
    traceback.print_exc(file=sys.stdout)
    pdb.post_mortem(sys.exc_traceback)
    raise
