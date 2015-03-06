#!/usr/bin/python
"""
Evaluates a seismic bulletin versus the ground truth.

For example the following command will evaluate the bulletin mysolution.data
versus the ground truth in test.data

./evaluate.py data/test.data mysolution.data

Author: Nimar Arora (feel free to modify, use, or redistribute)
"""
from __future__ import print_function, division

import sys

from numpy import mean, std

import mwmatching
# optimum checking doesn't work with non-integer weights, so disable it
mwmatching.CHECK_OPTIMUM = False

from util import *

def main():
  if len(sys.argv) != 3:
    print("Usage: evaluate.py gold.data guess.data",
          file = sys.stderr)
    sys.exit(1)

  gold_file_name, guess_file_name = sys.argv[1:]
  
  guessitr = iterate_episodes(guess_file_name)
  
  tot_gold, tot_mat, tot_guess, time_errs, dist_errs, mag_errs\
            = 0, 0, 0, [], [], []
  
  for gold in iterate_episodes(gold_file_name):

    try:
      
      guess = guessitr.next()

    except StopIteration:
      
      print("Guess data has fewer episodes than gold data!!",
            file = sys.stderr)
      
      break
    
    matchable_gold, indices = match_episodes(gold, guess)
    
    tot_gold += matchable_gold
    tot_guess += len(guess.events)
    
    for goldidx, guessidx in indices:
      
      tot_mat += 1
      
      goldev = gold.events[goldidx]
      guessev = guess.events[guessidx]

      time_errs.append(abs(goldev.time - guessev.time))
      dist_errs.append(compute_distance((goldev.lon, goldev.lat),
                                        (guessev.lon, guessev.lat)))
      mag_errs.append(abs(goldev.mag - guessev.mag))

  if tot_gold == 0:
    recall = 1
  else:
    recall = tot_mat / tot_gold

  if tot_guess == 0:
    prec = 1
  else:
    prec = tot_mat / tot_guess

  if prec == 0 or recall == 0:
    f1 = 0
  else:
    f1 = 2 * prec * recall / (prec + recall)

  print("{:d} matchable events, {:d} guess events, and {:d} matched"
        .format(tot_gold, tot_guess, tot_mat))
  print("Precision {:.1f} % , Recall {:.1f} % , F1 {:.1f}"
        .format(100 * prec, 100 * recall, 100 * f1))
  print("Time Errors mean {:.1f} std {:.1f}"
        .format(mean(time_errs), std(time_errs)))
  print("Dist Errors mean {:.1f} std {:.1f}"
        .format(mean(dist_errs), std(dist_errs)))
  print("Mag Errors mean {:.1f} std {:.1f}"
        .format(mean(mag_errs), std(mag_errs)))

def match_episodes(gold, guess, MINASSOC=2, DMAX=5, TMAX=50):
  """
  We want a max cardinality min cost matching for the provided events.

  Returns the number of matchable events and the matching indices.
  """
  ## we will first construct a list of edges between gold and guess events
  ## that are allowed to be matched.
  
  edges = []
  matchable_gold = 0
  
  for goldnum, goldev in enumerate(gold.events):
    
    # events with very few detections are unmatchable
    if len(gold.assocs[goldnum]) < MINASSOC:
      continue
    
    matchable_gold += 1
    
    for guessnum, guessev in enumerate(guess.events):
      # check the separation in space and time
      dist = compute_distance((goldev.lon, goldev.lat),
                              (guessev.lon, guessev.lat))
      timesep = abs(goldev.time - guessev.time)
      
      # note that the guess events are placed after the gold events hence
      # they are numbered as "len(gold.events) + guessnum"
      if dist <= DMAX and timesep <= TMAX:
        edges.append((goldnum, len(gold.events) + guessnum,
                      - dist / DMAX - timesep / TMAX))
  
  ## call the max-cardinality "max" weight matching, note that we used negative
  ## weights above so we effectively get max cardinality "min" weight matching
  
  mat = mwmatching.maxWeightMatching(edges, maxcardinality=True)
  
  ## convert the indices back to the original numbering
  
  indices = []
  for i in range(len(gold.events)):
    if i < len(mat) and mat[i] >= 0:
      assert(mat[i] >= len(gold.events))
      indices.append((i, mat[i] - len(gold.events)))

  return matchable_gold, indices

if __name__ == "__main__":
  main()
  
