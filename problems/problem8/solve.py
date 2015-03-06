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

VERBOSE = True                          # print some chatty output

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
    
    episode = solve_episode(physics, episode.detections)
    
    write_single_episode(episode, fp)
    
    fp.flush()
    
  fp.close()

def solve_episode(physics, detections):

  events, assocs = [], []

  skip_detnums = set()

  while True:
    
    solution = find_best_remaining_event(physics, detections, skip_detnums)

    if VERBOSE:
      print(solution)
    
    # we exit the function when we have found all events
    # for now put a threshold on 20 ... there are lots of approximations
    # here a full solution would integrate the event parameters by
    # using the prior density
    if solution is None or solution[0] < 20:
      # we will return an episode without detections since we don't need
      # to write these out
      return Episode(events, [], assocs)
    
    llratio, detnums, event = solution

    # we don't want these detections to be made available to other events
    for detnum in detnums:
      
      skip_detnums.add(detnum)
      
    events.append(event)
    assocs.append(detnums)

def find_best_remaining_event(physics, detections, skip_detnums):
  """
  Finds a single event that best explains detections with the highest
  overall loglikelihood ratio. We are ignoring prior information and
  mis-detection probability in this version.

  This returns an event object and a list of detection numbers associated to
  it.
  """
  # first go through each detection and get a candidate list

  candidates = []
  
  for cand_det in detections:

    for cand_ev in generate_det_candidates(physics, cand_det, 9):

      detnums, llratio = find_best_detections(physics, detections,
                                              skip_detnums, cand_ev)

      candidates.append((llratio, detnums, cand_ev))

  if not len(candidates):
    return None

  else:
    return max(candidates)

def find_best_detections(physics, detections, skip_detnums, event):
  """
  Finds the best set of detections, at most one per station, that may be
  associated to this event.

  Returns the list of detection numbers and a log likehood ratio.
  """

  sta_llratio = [None for _ in STATIONS]
  
  for detnum, det in enumerate(detections):
    
    if detnum in skip_detnums:
      continue

    llratio = compute_ev_det_loglike(physics, event, det) \
              - compute_noise_loglike(physics, det)
    
    # have we found a better detection at this station?
    
    if llratio > 0 and (sta_llratio[det.stanum] is None
                        or sta_llratio[det.stanum][1] < llratio):
      
      sta_llratio[det.stanum] = (detnum, llratio)
  
  # now return the detnums and total llratio
  
  return ([x[0] for x in sta_llratio if x is not None],
          sum(x[1] for x in sta_llratio if x is not None))

def generate_det_candidates(physics, det, numperturb):
  """
  Generate multiple candidates for an event from a detection by perturbing
  the detection's observed slowness and azimuth. The numperturb parameter
  controls the number of perturbations to each of them. Note, we will generate
  approximately numperturb^2 candidates.
  
  The size of each perturbation will be determined by the uncertainty
  in the azimuth and slowness.
  """
  # We will perturb the azimuth and slowness along evenly spaced percentile
  # values. However, we must always try the 50th percentile as well
  # because that is the unperturbed value.
  percentiles = list(np.linspace(.1, .9, numperturb))
  if .5 not in percentiles:
    percentiles.append(.5)

  # Using the percentiles compute the appropriate perturbation to
  # azimuth and slowness.
  delaz_values = [laplace.ppf(perc, loc = physics.mu_z[det.stanum],
                              scale = physics.theta_z[det.stanum])
                  for perc in percentiles]

  delslow_values = [laplace.ppf(perc, loc = physics.mu_s[det.stanum],
                                scale = physics.theta_s[det.stanum])
                    for perc in percentiles]

  candidates = []

  for delslow in delslow_values:

    for delaz in delaz_values:

      event = invert_detection(physics, det, delslow, delaz)

      candidates.append(event)
  
  return candidates


def invert_detection(physics, det, delslow, delaz):
  """
  Invert a single detection to produce an event that could have produced
  the detection.

  The idea is very simple. First invert the slowness to get a distance
  estimate, and then invert the distance estimate and the azimuth reading
  to get a potential location. Finally, estimate the magnitude and event time
  using the travel time.

  delslow and delaz are the perturbations to be applied to the slowness and
  azimiuth respectively to get the candidate location.
  """
  sta = STATIONS[det.stanum]
  
  dist = invert_slowness(det.slowness + delslow)
  
  lon, lat = invert_dist_azimuth((sta.lon, sta.lat), dist, det.azimuth + delaz)
  
  ttime = compute_travel_time(dist)
  
  evtime = det.time - ttime
  
  evmag = (log(det.amp) - physics.mu_a0[det.stanum]
           - physics.mu_a2[det.stanum] * dist) / physics.mu_a1[det.stanum]
  
  return Event(lon, lat, evmag, evtime)

def compute_ev_det_loglike(physics, event, det):
  """
  Compute the log-likelihood that the given event generated the given detection.

  The likelihood calculation includes the probability that the event
  was detected at the station and the probability of each an every
  detection attribute being generated by the event.
  """
  try:
    return _compute_ev_det_loglike(physics, event, det)

  # a floating point underflow or overflow indicates that the likelihood
  # is very low
  except FloatingPointError:
    return -inf
  
def _compute_ev_det_loglike(physics, event, det):

  # load the station tuple and compute basic event-station attributes like
  # distance, travel time, azimuth, and azimuth difference

  stanum = det.stanum
  station = STATIONS[stanum]
  
  dist = compute_distance((station.lon, station.lat),
                          (event.lon, event.lat))

  ttime = compute_travel_time(dist)
  
  sta_to_ev_az = compute_azimuth((station.lon, station.lat),
                                 (event.lon, event.lat))

  # the azimuth difference of observed to theoretical
  degdiff = compute_degdiff(sta_to_ev_az, det.azimuth)
  
  loglike = 0

  # detection probability
  
  detprob = logistic(physics.mu_d0[stanum]
                     + physics.mu_d1[stanum] * event.mag
                     + physics.mu_d2[stanum] * dist)

  loglike += log(detprob)

  # detection time

  loglike += laplace.logpdf(det.time,
                            event.time + ttime + physics.mu_t[stanum],
                            physics.theta_t[stanum])


  # detection azimuth
  
  loglike += laplace.logpdf(degdiff, physics.mu_z[stanum],
                            physics.theta_z[stanum])

  # detection slowness

  loglike += laplace.logpdf(det.slowness,
                            compute_slowness(dist) + physics.mu_s[stanum],
                            physics.theta_s[stanum])

  # detection amplitude

  loglike += norm.logpdf(log(det.amp),
                         physics.mu_a0[stanum]
                         + physics.mu_a1[stanum] * event.mag
                         + physics.mu_a2[stanum] * dist,
                         physics.sigma_a[stanum])
  
  return loglike

def compute_noise_loglike(physics, det):
  """
  Compute the log-likelihood that the detection was generated by noise.
  """
  assert (det.time > 0 and det.time < physics.T and det.amp > 0)

  loglike = 0

  loglike += -log(physics.T)          # detection time

  loglike += -log(360)                # detection azimuth

  loglike += -log(compute_slowness(0) - compute_slowness(180)) # slowness

  # detection amplitude
  loglike += cauchy.logpdf(log(det.amp),
                         physics.mu_f[det.stanum], physics.theta_f[det.stanum])
  
  return loglike

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
