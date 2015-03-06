#!/usr/bin/python
"""
Utility functions.

Author: Nimar Arora (feel free to modify, use, or redistribute)

Two of the geophysical functions, compute_distance and invert_dist_azimuth
were modified from a version copied from geopy:

http://pydoc.net/Python/geopy/0.94.1/geopy.distance/

The license terms for geopy are:

Copyright (c) 2006-2010 Brian Beck
Copyright (c) 2010-2015 GeoPy Project and individual contributors

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
from __future__ import print_function, division

## IMPORTS
import numpy as np
from scipy.stats import norm
import re
from collections import namedtuple

## OBJECTS

# the Physics object represents the underlying physics of the problem which
# is generated once for all episodes
Physics = namedtuple("Physics", [
  "T", "R", "lambda_e", "mu_m", "theta_m", "mu_d0", "mu_d1", "mu_d2",
  "mu_t", "theta_t", "mu_z", "theta_z", "mu_s", "theta_s",
  "mu_a0", "mu_a1", "mu_a2", "sigma_a",
  "lambda_f", "mu_f", "theta_f",
  ])

Station = namedtuple("Station", ["name", "lon", "lat"])

# the Episode object is a collection of events, detections, and associations
# (of events to detections)
Episode = namedtuple("Episode", ["events", "detections", "assocs"])

Event = namedtuple("Event", ["lon", "lat", "mag", "time"])

Detection = namedtuple("Detection", ["stanum", "time", "azimuth", "slowness",
                                     "amp"])

# Note: the assocs object in the episode is simply a list of the same
# length as the events. Each item in the list is a list of detection numbers.


## CONSTANTS

# create a list of Station objects from the list of tuples below
STATIONS = map(lambda x: Station(*x),
               [("ASAR", 133.9, -23.7),
                ("CMAR", 98.9, 18.5),
                ("FINES", 26.1, 61.4),
                ("ILAR", -146.9, 64.8),
                ("MKAR", 82.3, 46.8),
                ("SONM", 106.4, 47.8),
                ("STKA", 141.6, -31.9),
                ("TORD", 1.7, 13.1),
                ("WRA", 134.3, -19.9),
                ("ZALV", 84.8, 53.9),
                ])

def logistic(a):
  return 1/(1 + np.exp(-a))

def compute_travel_time(dist):
  """
  The travel time between two points given the distance in degrees.

  Returns time in seconds.

  This is the function I_T(\delta) in the description.
  """
  return -0.023 * dist ** 2 + 10.7 * dist + 5

def compute_slowness(dist):
  """
  The derivative of the travel time. Slowness is returned in seconds per
  degree given the distance in degrees.
  """
  return -0.046 * dist + 10.7

def invert_slowness(slow):
  """
  Returns a slowness value corresponding to the given distance.
  """
  return (slow - 10.7) / -0.046

def write_namedtuple(obj, filename):
  """
  write out a named tuple to file
  """
  fp = open(filename, "w")
  
  for field in obj._fields:
    print("{0} = {1}".format(field, getattr(obj, field)), file = fp)

  fp.close()

def read_namedtuple(typename, filename):
  """
  returns a named tuple that has been written using write_namedtuple
  the format of the file is simply
  
  <attribute name> = <python expression>
  
  """
  all_varnames, all_values = [], []
  
  for line in file(filename):
    result = re.search("(.*) = (.*)", line.strip())

    if result is None or len(result.groups()) != 2:
      raise ValueError("named tuple filename is corrupted")

    varname, value = result.groups()

    all_varnames.append(varname)
    all_values.append(eval(value))

  typeobj = namedtuple(typename, all_varnames)

  return typeobj(*all_values)

def write_episodes(episodes, filename):
  """
  Write out a list of episodes to a file
  """
  
  fp = open(filename, "w")
  
  print("Episodes:", file=fp, end="\n\n")

  for episode in episodes:
    
    write_single_episode(episode, fp)

  fp.close()

def write_single_episode(episode, fp):
    print("Events:", file=fp)

    for event in episode.events:

      print(" ".join(map(str, event)), file=fp)

    print("Detections:",  file=fp)

    for detection in episode.detections:
      
      print(" ".join(map(str, detection)), file=fp)

    print("Assocs:", file=fp)
    
    for evnum, assoc in enumerate(episode.assocs):

      for detnum in assoc:
        print("{0} {1}".format(evnum, detnum), file=fp)

    # blank line at the end of each episode
    print(file=fp)

def read_episodes(filename):
  """
  Read episodes from a file.
  """
  return [epi for epi in iterate_episodes(filename)]

def iterate_episodes(filename):
  """
  Iterating over the episodes is more memory efficient if you plan to
  discard them as soon as you have read them.
  """
  fp = open(filename)
  
  line1 = fp.readline().strip()
  line2 = fp.readline().strip()
  line3 = fp.readline().strip()
  
  if line1 != "Episodes:" or line2 != "" or line3 != "Events:":
    raise ValueError("Corrupted episodes files")
  
  while True:
    
    # read events
    
    events = []
    
    while True:
      line = fp.readline().strip()
      if line == "Detections:":
        break
      
      events.append(Event(*map(float, line.split())))
    
    # read detections
    
    detections = []
    
    while True:
      line = fp.readline().strip()
      if line == "Assocs:":
        break
      
      # convert the fields of the detection to floating point, but the
      # first field, which is the station num needs to be converted to integer
      det = map(float, line.split())
      det[0] = int(det[0])
      
      detections.append(Detection(*det))
      
    # read associations
    assocs = [[] for _ in events]
      
    while True:
      line = fp.readline().strip()
      if line == "":
        break
      
      evnum, detnum = map(int, line.split())
      
      assocs[evnum].append(detnum)
      
    # now we have a complete episode
    yield Episode(events, detections, assocs)
    
    # check if we have more episodes to read
    line = fp.readline().strip()
    
    if line == "Events:":
      continue
    
    else:
      break
  
def compute_distance(loc1, loc2):
  """
  Compute the great circle distance between two point on the earth's surface
  in degrees.
  loc1 and loc2 are pairs of longitude and latitude
  >>> compute_distance((10,0), (20, 0))
  10.0
  >>> compute_distance((10,0), (10, 45))
  45.0
  >>> int(compute_distance((-78, -12), (-10.25, 52)))
  86
  >>> compute_distance((132.86521, -0.45606493), (132.86521, -0.45606493)) == 0
  True
  """
  lng1, lat1 = np.radians(loc1[0]), np.radians(loc1[1])
  lng2, lat2 = np.radians(loc2[0]), np.radians(loc2[1])
  
  sin_lat1, cos_lat1 = np.sin(lat1), np.cos(lat1)
  sin_lat2, cos_lat2 = np.sin(lat2), np.cos(lat2)
  
  delta_lng = lng2 - lng1
  cos_delta_lng, sin_delta_lng = np.cos(delta_lng), np.sin(delta_lng)
  
  # From http://en.wikipedia.org/wiki/Great_circle_distance:
  #   Historically, the use of this formula was simplified by the
  #   availability of tables for the haversine function. Although this
  #   formula is accurate for most distances, it too suffers from
  #   rounding errors for the special (and somewhat unusual) case of
  #   antipodal points (on opposite ends of the sphere). A more
  #   complicated formula that is accurate for all distances is: (below)
  
  angle = np.arctan2(np.sqrt((cos_lat2 * sin_delta_lng) ** 2 +
                             (cos_lat1 * sin_lat2 -
                              sin_lat1 * cos_lat2 * cos_delta_lng) ** 2),
                     sin_lat1 * sin_lat2 + cos_lat1 * cos_lat2 * cos_delta_lng)
  
  return np.degrees(angle)

def compute_azimuth(loc1, loc2):
  """
  Angle in degrees measured clockwise from north starting at
  loc1 towards loc2. loc1 and loc2 are (longitude, latitude) in degrees.

  See https://en.wikipedia.org/wiki/Great-circle_navigation.

  However, we want north = 0 degrees,
                   east = 90 degrees,
                   south = 180 degrees, and
                   west = 270 degrees.

  Return an angle in 0 to 360 degrees.
  
  >>> int(compute_azimuth((10,0), (20, 0)))
  90
  >>> int(compute_azimuth((20,0), (10, 0)))
  270
  >>> int(compute_azimuth((10,0), (10, 45)))
  0
  >>> compute_azimuth((10, 45), (10, 0))
  180.0
  >>> int(compute_azimuth((133.9, -23.665), (132.6, -.83)))
  356
  """
  lng1, lat1 = np.radians(loc1[0]), np.radians(loc1[1])
  lng2, lat2 = np.radians(loc2[0]), np.radians(loc2[1])
  
  delta_lon = lng2 - lng1;
  
  y = np.sin(delta_lon)
  
  x = np.cos(lat1) * np.tan(lat2) - np.sin(lat1)*np.cos(delta_lon);
  
  azi = np.degrees(np.arctan2(y, x))

  # azi is now in range -180/180
  # the following trick brings it in the 0 - 360 range

  return (azi + 360) % 360
  
def invert_dist_azimuth(loc1, dist, azi):
  """
  Return the location loc2=(longitude, latitude) such that
  compute_distance(loc1, loc2) == dist  and  compute_azimuth(loc1, loc2) = azi.

  Or, in other words if we move along the great-circle from 'loc1' along the
  azimuth 'azi' for a distance of 'dist' degrees then we will reach loc2.

  >>> map(int, invert_dist_azimuth((10, 0), 10.049369393181079,
  ... 95.740074136412659))
  [20, 0]

  """
  # convert everything to radians
  lng1, lat1 = np.radians(loc1[0]), np.radians(loc1[1])
  dist, azi = np.radians(dist), np.radians(azi)

  lat2 = np.arcsin(np.sin(lat1) * np.cos(dist)
                   + np.cos(lat1) * np.sin(dist) * np.cos(azi) )

  lng2 = lng1 + np.arctan2(np.sin(azi) * np.sin(dist) * np.cos(lat1),
                           np.cos(dist) - np.sin(lat1) * np.sin(lat2))

  return np.degrees(lng2), np.degrees(lat2)

def compute_degdiff(angle1, angle2):
  """
  The difference of two angles given in degrees. The answer is an angle from
  -180 to 180. Positive angles imply angle2 is clockwise from angle1 and -ve
  angles imply counter-clockwise.
  
  >>> int(compute_degdiff(40, 30))
  -10
  >>> int(compute_degdiff(30, 40))
  10
  >>> int(compute_degdiff(361, 40))
  39
  >>> int(compute_degdiff(40, 361))
  -39
  >>> compute_degdiff(40,250)
  -150
  >>> compute_degdiff(40,200)
  160
  >>> compute_degdiff(40, 219)
  179
  >>> compute_degdiff(40, 220)
  180
  >>> compute_degdiff(40, 221)
  -179
  """
  # bring the angle into the 0 to 360 range
  delta = ((angle2 - angle1) + 360) % 360
  # angles above 180 need to be shifted down by 360 degrees so that 181 is -179
  # 200 is -160 etc.
  return delta - (delta > 180) * 360

def mvar_norm_fit(data):
  """
  Fits a multivariate normal using row vectors as data.

  Returns a mean vector (also a row vector) and a covariance matrix.
  """
  assert(len(data) > 1)
  
  mu = sum(data) / len(data)

  # data should be a row vectors
  assert(len(mu.shape)==1)

  # ML estimator of covariance

  Sigma = sum( np.outer(x - mu, x-mu) for x in data ) / len(data)

  return mu, Sigma

def mvar_norm_sample(mu, Sigma):
  """
  Generates a random variable from a multivariate normal distribution with
  given mean vector and covariance matrix.
  
  Returns a vector of the same shape as the mean vector
  
  see http://en.wikipedia.org/wiki/Multivariate_normal_distribution  
  """
  # mu better be a row vector
  assert(len(mu.shape)==1)

  # the square root (roughly speaking) of Sigma
  # i.e. root * root.T = Sigma
  root = np.linalg.cholesky(Sigma)

  z = np.array([norm.rvs() for _ in mu])
  
  return mu + np.dot(root, z)
  
# if this file is run stand alone it tests the doc strings

if __name__ == "__main__":
  import doctest
  doctest.testmod()
