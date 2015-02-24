"""
Sample code to compare 1-d seismic bulletin using max-cardinality min-weight
matching.

A matching is possible between two events that are separated by at most DMAX
in distance and TMAX in time. The weight of the matching is the
normalized, absolute difference in space and time between matching events
(normalized by DMAX and TMAX).

For example if an event at (x,t) is matched to an event at (y, r) the weight
of this match is abs(x-y)/DMAX + abs(t-r)/TMAX.

Author: Nimar Arora (feel free to modify, use, or redistribute)
"""
from __future__ import print_function, division

import mwmatching
# optimum checking doesn't work with non-integer weights, so disable it
mwmatching.CHECK_OPTIMUM = False        


def evaluate(gold_events, guess_events, DMAX=.01, TMAX=.01):
  """
  We want a max cardinality min cost matching for the provided events.
  
  Events are assumed to be a list of tuples, where each tuple
  is (x-coord, time-coord)
  
  Returns precision, recall, f1-score, and a matching.
  Where a matching is a list of pairs of gold, guess indices
  """
  
  ## we will first construct a list of edges between gold and guess events
  ## that are allowed to be matched.
  
  edges = []
  for goldnum, gold in enumerate(gold_events):
    for guessnum, guess in enumerate(guess_events):
      # check the separation in space and time
      xdist = abs(gold[0] - guess[0])
      tdist = abs(gold[1] - guess[1])
      # note that the guess events are placed after the gold events hence
      # they are numbered as "len(gold_events) + guessnum"
      if xdist <= DMAX and tdist <= TMAX:
        edges.append((goldnum, len(gold_events)+guessnum,
                      -xdist / DMAX - tdist / TMAX))

  ## call the max-cardinality "max" weight matching, note that we used negative
  ## weights above so we effectively get max cardinality "min" weight matching
  
  mat = mwmatching.maxWeightMatching(edges, maxcardinality=True)
  
  ## convert the indices back to the original numbering
  
  indices = []
  for i in range(len(gold_events)):
    if i < len(mat) and mat[i] >= 0:
      assert(mat[i] >= len(gold_events))
      indices.append((i, mat[i] - len(gold_events)))

  ## now compute the f1, precision, and recall
  
  if not len(gold_events):
    recall = 1.0
  else:
    recall = len(indices) / len(gold_events)

  if not len(guess_events):
    prec = 1.0
  else:
    prec = len(indices) / len(guess_events)

  if prec == 0.0 or recall == 0:
    f1 = 0.
  else:
    f1 = 2 * prec * recall / (prec + recall)
  
  return prec, recall, f1, indices

if __name__ == "__main__":
  ## construct some sample input, and check the evaluation
  
  gold = [(.034, .015), (.048, .030),
          (.100, .100),
          (.6, .6)]
  
  guess = [(.039, .021), (.025, .006),
           (.101, .109), (.102, .105),
           (.5, .5)]

  prec, recall, f1, indices = evaluate(gold, guess)

  print ("Matching", indices)
  print ("Precision", prec)
  print ("Recall", recall)
  print ("F1", f1)

  """
  Expected output:
  
  Matching [(0, 1), (1, 0), (2, 3)]
  Precision 0.6
  Recall 0.75
  F1 0.666666666667
  """
