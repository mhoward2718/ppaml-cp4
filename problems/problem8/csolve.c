/*
  Solves the seismic-2d episodes.

  Author: Nimar Arora (feel free to modify, use, or redistribute)
  Copyright 2015.
*/
#include <Python.h>
#include <assert.h>
#include <gsl/gsl_randist.h>

#define VERBOSE
// #define VERYVERBOSE

/*
 * number of perturbations to apply to the slowness and azimuth of a detection
 * in order to get a list of potential events
 */
#define NUMPERTURB 9

/* arbitrary threshold on the log-likelihood ratio. This is to be replaced
 * by a more sophisticated calculation based on integrating over event priors.
 */
#define LLRATIO_THRESHOLD 20

/* Data types copied from util.py */

typedef struct Physics_t
{
  double T;
  double R;
  double lambda_e;
  double mu_m;
  double theta_m;
  double gamma_m;
  double * mu_d0;
  double * mu_d1;
  double * mu_d2;
  double * mu_t;
  double * theta_t;
  double * mu_z;
  double * theta_z;
  double * mu_s;
  double * theta_s;
  double * mu_a0;
  double * mu_a1;
  double * mu_a2;
  double * sigma_a;
  double * lambda_f;
  double * mu_f;
  double * theta_f;

} Physics_t;

typedef struct Station_t
{
  /* the memory for the name is a borrowed reference it must not be freed */
  char * name;
  double lon;
  double lat;
} Station_t;

typedef struct Detection_t
{
  int staidx;
  double time;
  double azimuth;
  double slowness;
  double amp;

  /* cached values */
  double logamp;
  double noise_loglike;

} Detection_t;

typedef struct Event_t
{
  double lon;
  double lat;
  double mag;
  double time;

  /* The size of sta_det is the number of stations. Each entry is either -1 
     or the detidx at that station */
  int * sta_det;         
  double sum_llratio;
  
} Event_t;

/* Data conversion routines to read data in and out of Python */

static double read_float_attr(PyObject * record_obj, char * attr_name)
/* read a single floating point attribute from a Python record and return
 * a C double */
{
  /* new reference */
  PyObject * item_obj = PyObject_GetAttrString(record_obj, attr_name);
  
  double retval = PyFloat_AsDouble(item_obj);

  Py_DECREF(item_obj);

#ifdef VERYVERBOSE
  printf(" %s = %lf\n", attr_name, retval);
#endif

  return retval;
}

static double * read_list_float_attr(PyObject * record_obj, char * attr_name)
/* Read a list of floating point values from a Python record given the
 * attribute name and return a C list of doubles. The C list has to be
 * deallocated with free.
 */
{
  /* new reference */
  PyObject * item_list_obj = PyObject_GetAttrString(record_obj, attr_name);
  
  int num_item = (int) PySequence_Length(item_list_obj);
  
  double * retval = (double *) malloc(num_item * sizeof(double));

  int idx;
  
  for (idx = 0; idx < num_item; idx++)
  {
    /* new reference */
    PyObject * item_obj = PySequence_GetItem(item_list_obj, (Py_ssize_t) idx);

    retval[idx] = PyFloat_AsDouble(item_obj);

    Py_DECREF(item_obj);
  }
  
  Py_DECREF(item_list_obj);

#ifdef VERYVERBOSE
  printf(" %s = [", attr_name);
  for (idx = 0; idx < num_item; idx++)
    printf("%lf, ", retval[idx]);
  printf("]\n");
#endif

  return retval;
}

static Physics_t * read_physics(PyObject * physics_obj)
/*
 * Read a Python physics named tuple and return a C Physics_t object.
 * The caller is responsible to free the returned object by calling 
 * free_physics(...).
 */
{
  Physics_t * physics = (Physics_t *) malloc(sizeof(Physics_t));

  physics->T = read_float_attr(physics_obj, "T");
  physics->R = read_float_attr(physics_obj, "R");
  physics->lambda_e = read_float_attr(physics_obj, "lambda_e");
  physics->mu_m = read_float_attr(physics_obj, "mu_m");
  physics->theta_m = read_float_attr(physics_obj, "theta_m");
  physics->gamma_m = read_float_attr(physics_obj, "gamma_m");

  physics->mu_d0 = read_list_float_attr(physics_obj, "mu_d0");
  physics->mu_d1 = read_list_float_attr(physics_obj, "mu_d1");
  physics->mu_d2 = read_list_float_attr(physics_obj, "mu_d2");

  physics->mu_t = read_list_float_attr(physics_obj, "mu_t");
  physics->theta_t = read_list_float_attr(physics_obj, "theta_t");

  physics->mu_z = read_list_float_attr(physics_obj, "mu_z");
  physics->theta_z = read_list_float_attr(physics_obj, "theta_z");

  physics->mu_s = read_list_float_attr(physics_obj, "mu_s");
  physics->theta_s = read_list_float_attr(physics_obj, "theta_s");

  physics->mu_a0 = read_list_float_attr(physics_obj, "mu_a0");
  physics->mu_a1 = read_list_float_attr(physics_obj, "mu_a1");
  physics->mu_a2 = read_list_float_attr(physics_obj, "mu_a2");
  physics->sigma_a = read_list_float_attr(physics_obj, "sigma_a");
  
  physics->lambda_f = read_list_float_attr(physics_obj, "lambda_f");

  physics->mu_f = read_list_float_attr(physics_obj, "mu_f");
  physics->theta_f = read_list_float_attr(physics_obj, "theta_f");
  
  return physics;
}

static void free_physics(Physics_t * physics)
/* free the object previously allocated by read_physics(...) */
{
  free(physics->mu_d0);
  free(physics->mu_d1);
  free(physics->mu_d2);
  free(physics->mu_t);
  free(physics->theta_t);
  free(physics->mu_z);
  free(physics->theta_z);
  free(physics->mu_s);
  free(physics->theta_s);
  free(physics->mu_a0);
  free(physics->mu_a1);
  free(physics->mu_a2);
  free(physics->sigma_a);
  free(physics->lambda_f);
  free(physics->mu_f);
  free(physics->theta_f);
  
  free(physics);
}

static Station_t * read_sta_list(PyObject * sta_list_obj, int * out_num_sta)
/* Read a Python list of Station tuples (name, lon, lat) and return a C
 * list of Station_t objects that has to be freed with free_sta_list */
{
  int num_sta = (int) PySequence_Length(sta_list_obj);
  
  Station_t * sta_list = (Station_t *)malloc(num_sta * sizeof(Station_t));

  int staidx;
  
  for (staidx = 0; staidx < num_sta; staidx ++)
  {
    /* new reference */
    PyObject * sta_obj = PySequence_GetItem(sta_list_obj, (Py_ssize_t) staidx);

    PyArg_ParseTuple(sta_obj, "sdd", &sta_list[staidx].name,
                     &sta_list[staidx].lon, &sta_list[staidx].lat);

#ifdef VERYVERBOSE
    printf("name %s lon %lf lat %lf\n", sta_list[staidx].name,
           sta_list[staidx].lon, sta_list[staidx].lat);
#endif

    Py_DECREF(sta_obj);
  }

  *out_num_sta = num_sta;
  return sta_list;
}

static void free_sta_list(Station_t * sta_list)
{
  free(sta_list);
}

static Detection_t * read_det_list(PyObject * det_list_obj, int * out_num_det)
/* Read a Python list of Detection tuples and return a C
 * list of Detection_t objects that has to be freed with free_det_list */
{
  int num_det = (int) PySequence_Length(det_list_obj);
  
  Detection_t * det_list = (Detection_t *)malloc(num_det * sizeof(Detection_t));

  int detidx;
  
  for (detidx = 0; detidx < num_det; detidx ++)
  {
    /* new reference */
    PyObject * det_obj = PySequence_GetItem(det_list_obj, (Py_ssize_t) detidx);

    PyArg_ParseTuple(det_obj, "idddd", &det_list[detidx].staidx,
                     &det_list[detidx].time, &det_list[detidx].azimuth,
                     &det_list[detidx].slowness, &det_list[detidx].amp);


#ifdef VERYVERBOSE
    printf("staidx %d time %lf azimuth %lf slowness %lf amp %lf\n", 
           det_list[detidx].staidx,
           det_list[detidx].time, det_list[detidx].azimuth, 
           det_list[detidx].slowness, det_list[detidx].amp);
#endif

    Py_DECREF(det_obj);
  }

  *out_num_det = num_det;
  return det_list;
}

static void free_det_list(Detection_t * det_list)
{
  free(det_list);
}

static PyObject * make_event_obj(const Event_t * event)
/* returns a Python tuple consisting of the event fields */
{
  return Py_BuildValue("dddd", event->lon, event->lat, event->mag, event->time);
}

static PyObject * make_assoc_obj(const Event_t * event, int num_sta)
/* returns a Python list of detnums associated to the event */
{
  PyObject * assoc_obj = PyList_New(0);
 
  int staidx;
  for (staidx=0; staidx<num_sta; staidx++)
  {
    int detnum = event->sta_det[staidx];
    
    if (detnum >= 0)
    {
      PyObject * item = Py_BuildValue("i", detnum); /* new reference */
      
      PyList_Append(assoc_obj, item);
      
      Py_DECREF(item);
    }
  }

  return assoc_obj;
}


static void free_event(Event_t * event)
{
  free(event->sta_det);
  free(event);
}

#ifdef VERBOSE

static void print_event(const Event_t * event, int num_sta)
{
  printf("lon %lf lat %lf mag %lf time %lf [", event->lon, event->lat,
         event->mag, event->time);

  int staidx;
  
  for (staidx=0; staidx<num_sta; staidx++)
  {
    int detidx = event->sta_det[staidx];
    
    if (detidx >= 0)
      printf("%d ", detidx);
  }
  
  printf("] sum_llratio %lf", event->sum_llratio);
}

#endif

/* end of data conversion routines */

/* geometry and geography routines */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static inline double to_degrees(double radians)
{
  return (radians / M_PI) * 180;
}

static inline double to_radians(double degrees)
{
  return (degrees / 180) * M_PI ;
}

static inline double square(double val)
{
  return val * val;
}

static inline double compute_travel_time(double dist)
{
  return (-0.023 * square(dist) + 10.7 * dist + 5);
}

static inline double compute_slowness(double dist)
{
  return (-0.046 * dist + 10.7);
}

static inline double invert_slowness(double slow)
{
  return ((slow - 10.7) / -0.046);
}
 
static double compute_distance(double lon1, double lat1, double lon2,
                               double lat2)
/*
  Compute the great circle distance between two point on the earth's surface
  in degrees. Both the input and the output values are in degrees.

  Returns a number between 0 and 180
*/
{
  lon1 = to_radians(lon1);
  lat1 = to_radians(lat1);

  lon2 = to_radians(lon2);
  lat2 = to_radians(lat2);
  
  double sin_lat1 = sin(lat1);
  double cos_lat1 = cos(lat1);

  double sin_lat2 = sin(lat2);
  double cos_lat2 = cos(lat2);
  
  double delta_lon = lon2 - lon1;
  
  double cos_delta_lon = cos(delta_lon);
  double sin_delta_lon = sin(delta_lon);
  
  double angle = atan2(sqrt(square(cos_lat2 * sin_delta_lon) +
                            square(cos_lat1 * sin_lat2 -
                                   sin_lat1 * cos_lat2 * cos_delta_lon)),
                       sin_lat1 * sin_lat2 
                       + cos_lat1 * cos_lat2 * cos_delta_lon);
  
  return to_degrees(angle);
}

static double compute_azimuth(double lon1, double lat1, double lon2,
                              double lat2)
/*
  Angle in degrees measured clockwise from north starting at
  loc1 towards loc2. loc1 and loc2 are (longitude, latitude) in degrees.

  See https://en.wikipedia.org/wiki/Great-circle_navigation.

  However, we want north = 0 degrees,
                   east = 90 degrees,
                   south = 180 degrees, and
                   west = 270 degrees.

  Return an angle in 0 to 360 degrees.
  
*/
{
  lon1 = to_radians(lon1);
  lat1 = to_radians(lat1);

  lon2 = to_radians(lon2);
  lat2 = to_radians(lat2);
  
  double delta_lon = lon2 - lon1;
  
  double y = sin(delta_lon);
  
  double x = cos(lat1) * tan(lat2) - sin(lat1) * cos(delta_lon);
  
  double azi = to_degrees(atan2(y, x));

  /*
    azi is now in range -180/180
    the following trick brings it in the 0 - 360 range
  */

  return fmod(azi + 360 , 360);
}

static void invert_dist_azimuth(double lon1, double lat1, double dist,
                                double azi, double * out_lon2,
                                double * out_lat2)
/*
  Return the location loc2=(longitude, latitude) such that
  compute_distance(loc1, loc2) == dist  and  compute_azimuth(loc1, loc2) = azi.

  Or, in other words if we move along the great-circle from 'loc1' along the
  azimuth 'azi' for a distance of 'dist' degrees then we will reach loc2.

  Example:
     invert_dist_azimuth(10, 0, 10.049369393181079, 95.740074136412659, ..)
       --> 20, 0
*/
{
  /* convert everything to radians */
  lon1 = to_radians(lon1);
  lat1 = to_radians(lat1);
  dist = to_radians(dist);
  azi = to_radians(azi);

  double lat2 = asin(sin(lat1) * cos(dist) + cos(lat1) * sin(dist) * cos(azi));

  double lon2 = lon1 + atan2(sin(azi) * sin(dist) * cos(lat1),
                             cos(dist) - sin(lat1) * sin(lat2));

  *out_lon2 = to_degrees(lon2);
  *out_lat2 = to_degrees(lat2);
}

static inline double compute_degdiff(double angle1, double angle2)
/*
  The difference of two angles given in degrees. The answer is an angle from
  -180 to 180. Positive angles imply angle2 is clockwise from angle1 and -ve
  angles imply counter-clockwise.
  
  compute_degdiff(40, 30) -> -10

  compute_degdiff(30, 40) -> 10

  compute_degdiff(361, 40) -> 39

  compute_degdiff(40, 361) -> -39

  compute_degdiff(40,250) -> -150

  compute_degdiff(40,200) -> 160

  compute_degdiff(40, 219) -> 179

  compute_degdiff(40, 220) -> 180

  compute_degdiff(40, 221) -> -179
*/
{
  /* bring the angle difference into the 0 to 360 range */
  double delta = fmod(angle2 - angle1 + 360, 360);
  
  /*
    angles above 180 need to be shifted down by 360 degrees so that 181 is -179
    200 is -160 etc.
  */
  if (delta > 180)
    delta -= 360;
  
  return delta;
}

/* end geographical functions */

/* begin statistical functions */

static inline double logistic(double a)
{
  return 1.0 / ( 1.0 + exp(-a) );
}

static inline double laplace_ppf(double perc, double loc, double scale)
/*
 * Given the location and scale of a laplace compute the value that
 * correponds to a required percentile. This is simply the inverse of the
 * cumulative function.
 */
{
  assert((perc > 0) && (perc < 1.0));
  
  if (perc > .5)
    return loc - scale * log(2 * (1 - perc));

  else if (perc < .5)
    return loc + scale * log(2 * perc);

  else
    return loc;
}

static inline double laplace_logpdf(double x, double loc, double scale)
{
  return -log(2 * scale) - fabs(x - loc) / scale;
}

static inline double norm_logpdf(double x, double mu, double sigma)
{
  return - log(sigma) - 0.5 * log(2 * M_PI) - 0.5 * square( (x - mu) / sigma);
}

static inline double cauchy_logpdf(double x, double loc, double scale)
{
  return -log(M_PI) - log(scale) - log(1 + square((x-loc)/scale));
}

/* end statistical functions */

static Event_t * invert_detection(const Physics_t * physics, int num_sta,
                                  const Station_t * sta_list,
                                  const Detection_t * det,
                                  double delslow, double delaz)
/*
  Invert a single detection to produce an event that could have produced
  the detection.

  The idea is very simple. First invert the slowness to get a distance
  estimate, and then invert the distance estimate and the azimuth reading
  to get a potential location. Finally, estimate the magnitude and event time
  using the travel time.

  delslow and delaz are the perturbations to be applied to the slowness and
  azimiuth respectively to get the candidate location.
*/
{
  /* allocate an event and set its sta_det array to -1. Note: memset
   * can be used to set an int array to -1 because setting each byte
   * of an integer to -1 or 0xff is equivalent to setting the whole
   * integer to -1 */
  Event_t * event = (Event_t *) malloc(sizeof(Event_t));
  event->sta_det = (int *)malloc(num_sta * sizeof(int));
  memset(event->sta_det, -1, num_sta * sizeof(int));

  event->sum_llratio = 0; /* no events; so log-likelihood ratio is zero */

  int staidx = det->staidx;
  
  const Station_t * sta = sta_list + staidx;
  
  double dist = invert_slowness(det->slowness + delslow);
  
  invert_dist_azimuth(sta->lon, sta->lat, dist, det->azimuth + delaz,
                      &event->lon, &event->lat);
  
  event->time = det->time - compute_travel_time(dist);
  
  event->mag = (det->logamp - physics->mu_a0[staidx]
                - physics->mu_a2[staidx] * dist) / physics->mu_a1[staidx];

  /* don't exceed the maximum magnitude */
  if (event->mag > physics->gamma_m)
    event->mag = physics->gamma_m;
  
  return event;
}

static double compute_noiseloglike(const Physics_t * physics,
                                   int num_sta, const Station_t * sta_list,
                                   const Detection_t * det)
/* compute the log-likelihood that the detection was generated by a
 * noise process at its station */
{
  assert((det->time > 0) && (det->time < physics->T) && (det->amp > 0));

  double loglike = 0;

  loglike += -log(physics->T);          /* detection time */

  loglike += -log(360);                /* detection azimuth */

  loglike += -log(compute_slowness(0) - compute_slowness(180)); /* slowness */

  /* detection amplitude */
  loglike += cauchy_logpdf(det->logamp, physics->mu_f[det->staidx],
                           physics->theta_f[det->staidx]);
  
  return loglike;
}

static void init_det_list(const Physics_t * physics,
                          int num_sta, const Station_t * sta_list,
                          int num_det, Detection_t * det_list)
/* initialize some cached values inside the Detection objects */
{
  int detidx;
  for (detidx=0; detidx<num_det; detidx++)
  {
    Detection_t * det = det_list + detidx;

    /* initialize logamp to avoid unnecessary repeated computation */
    det->logamp = log(det->amp);

    /* also compute and store the noise loglikelihood as this is
     * referenced a lot */
    det->noise_loglike = compute_noiseloglike(physics, num_sta, sta_list, det);
  }
}

static double compute_ev_det_loglike(const Physics_t * physics,
                                     int num_sta, const Station_t * sta_list,
                                     const Event_t * event,
                                     const Detection_t * det)
/* returns the log likelihood that the event generated that given detection */
{
  /* first compute basic event-station attributes like distance,
   * travel time, azimuth, and azimuth difference
   */
  int staidx = det->staidx;
  
  const Station_t * sta = sta_list + staidx;
  
  double dist = compute_distance(sta->lon, sta->lat, event->lon, event->lat);

  double ttime = compute_travel_time(dist);
  
  double sta_to_ev_az = compute_azimuth(sta->lon, sta->lat, event->lon,
                                        event->lat);
  
  /* the azimuth difference of observed vs. theoretical */
  double degdiff = compute_degdiff(sta_to_ev_az, det->azimuth);
  
  double loglike = 0;

  /* detection probability */
  
  double detprob = logistic(physics->mu_d0[staidx]
                            + physics->mu_d1[staidx] * event->mag
                            + physics->mu_d2[staidx] * dist);
  

  loglike += log(detprob);

  /* detection time */

  loglike += laplace_logpdf(det->time,
                            event->time + ttime + physics->mu_t[staidx],
                            physics->theta_t[staidx]);
  


  /* detection azimuth */
  
  loglike += laplace_logpdf(degdiff, physics->mu_z[staidx],
                            physics->theta_z[staidx]);
  
  
  /* detection slowness */

  loglike += laplace_logpdf(det->slowness,
                            compute_slowness(dist) + physics->mu_s[staidx],
                            physics->theta_s[staidx]);

  /* detection amplitude */

  loglike += norm_logpdf(det->logamp,
                         physics->mu_a0[staidx]
                         + physics->mu_a1[staidx] * event->mag
                         + physics->mu_a2[staidx] * dist,
                         physics->sigma_a[staidx]);
  
  return loglike;
}

static void find_best_detections(const Physics_t * physics,
                                 int num_sta, const Station_t * sta_list,
                                 int num_det, const Detection_t * det_list,
                                 const int * skip_detnums, Event_t * event)
/*
  Finds the best set of detections, at most one per station, that may be
  associated to this event. We are not allowed to use any detnum whose flag
  is set to True in skip_detnums.

  note: event->sta_det is assumed to be all -1 and
  event->sum_llratio is assumed to be zero when this function is called
  
  Modifies the event->sta_det field and updates event->sum_llratio
*/
{
  double * sta_best_llratio = (double *) malloc(num_sta * sizeof(double));

  int detidx;
  for(detidx = 0; detidx < num_det; detidx ++ )
  {
    if (skip_detnums[detidx])
      continue;
    
    const Detection_t * det = det_list + detidx;
    
    double llratio = compute_ev_det_loglike(physics, num_sta, sta_list,
                                            event, det) - det->noise_loglike;
    
    /* have we found a better detection at this station */
    int staidx = det->staidx;
    
    if ((llratio > 0) && ((event->sta_det[staidx] < 0) 
                          || (sta_best_llratio[staidx] < llratio)))
    {
      event->sta_det[staidx] = detidx;
      sta_best_llratio[staidx] = llratio;
    }
  }
  
  /* finally add up the llratio of the best event at each station */
  int staidx;
  for (staidx=0; staidx<num_sta; staidx++)
  {
    if (event->sta_det[staidx] >= 0)
      event->sum_llratio += sta_best_llratio[staidx];
  }
 
  free(sta_best_llratio);
  
}



static Event_t * find_best_remaining_event(const Physics_t * physics,
                                           int num_sta,
                                           const Station_t * sta_list,
                                           int num_det,
                                           const Detection_t * det_list,
                                           int * skip_detnums)
/*
  Finds a single event that best explains detections with the highest
  overall loglikelihood ratio. We are ignoring prior information and
  mis-detection probability in this version.

  This returns an event object with a list of detection numbers associated to
  it.

  The boolean array skip_detnums is updated to reflect the detections
  associated with the returned event object.
*/
{
  Event_t * best_event = NULL;
  
  /* The perturbation that we will apply to the slowness and azimuth
   * depend on the percentiles of their relevant distributions. Here
   * we first compute those percentiles between 10% and 90% */
  double PERCENTILES[NUMPERTURB];
  int percidx;
  for (percidx=0; percidx<NUMPERTURB; percidx++)
  {
    PERCENTILES[percidx] = .1 + percidx * (.9 - .1) / (NUMPERTURB-1);
  }
  
  /* go through each detection and invert it to get a suitable candidate */
  int inv_detidx;
  for (inv_detidx=0; inv_detidx < num_det; inv_detidx++)
  {
    const Detection_t * inv_det = det_list + inv_detidx;
    int staidx = inv_det->staidx;

#ifdef VERYVERBOSE
    printf("detidx %d noise log-likelihood %lf\n", inv_detidx,
           inv_det->noise_loglike);
#endif

    int sloidx;
    for (sloidx=0; sloidx<NUMPERTURB; sloidx++)
    {
      double delslow = laplace_ppf(PERCENTILES[sloidx], physics->mu_s[staidx],
                                   physics->theta_s[staidx]);
        
      int azidx;
      for (azidx=0; azidx<NUMPERTURB; azidx++)
      {
        double delaz = laplace_ppf(PERCENTILES[azidx], physics->mu_z[staidx],
                                   physics->theta_z[staidx]);
        
        Event_t * event = invert_detection(physics, num_sta, sta_list,
                                           inv_det, delslow, delaz);

        /* find the best set of detections for this event */
        find_best_detections(physics, num_sta, sta_list, num_det, det_list,
                             skip_detnums, event);

#ifdef VERYVERBOSE        
        printf("Inverted Event: ");
        print_event(event, num_sta);
        printf("\n");
#endif

        if ((event->sum_llratio > LLRATIO_THRESHOLD)
            && ((NULL == best_event) 
                || (event->sum_llratio > best_event->sum_llratio)))
        {
          if (NULL != best_event)
            free_event(best_event);
          
          best_event = event;
        }
        
        else
          free_event(event);
      }
    }
    
  }

  /* update skip_detnums with the detections of this best event */
  if (NULL != best_event)
  {
    int staidx;
    for (staidx=0; staidx < num_sta; staidx++)
    {
      if (best_event->sta_det[staidx] >= 0)
      {
        int detidx = best_event->sta_det[staidx];

        assert(detidx < num_det);
        
        skip_detnums[detidx] = 1;
      }
    }
  }
  
  return best_event;
}


static PyObject * csolve_solve_episode(PyObject * self, PyObject * args)
{
  PyObject * physics_obj;
  PyObject * sta_list_obj;
  PyObject * det_list_obj;

  /* note: the physics object and the in_episode object are borrowed references
   * for the lifetime of this function call
   */
  if (!PyArg_ParseTuple(args, "OOO", &physics_obj, &sta_list_obj,
                        &det_list_obj))
    return NULL;

  Physics_t * physics = read_physics(physics_obj);

  int num_sta;
  Station_t * sta_list = read_sta_list(sta_list_obj, &num_sta);

  int num_det;
  Detection_t * det_list = read_det_list(det_list_obj, &num_det);

  /* initialize some cached values in the detections */
  init_det_list(physics, num_sta, sta_list, num_det, det_list);
  
  /* allocate a boolean list indicating if a detection number is to be skipped
   * initially this list is all False (because calloc clears memory to 0) */
  int * skip_detnums = (int *)calloc(num_det,  sizeof(int));

  PyObject * events_list_obj = PyList_New(0);
  PyObject * assocs_list_obj = PyList_New(0);

  while (1)
  {
    Event_t * event = find_best_remaining_event(physics, num_sta, sta_list,
                                                num_det, det_list,
                                                skip_detnums);

    if (NULL == event)
      break;

#ifdef VERBOSE
    print_event(event, num_sta);
    printf("\n");
#endif

    PyObject * event_obj = make_event_obj(event);

    PyList_Append(events_list_obj, event_obj);
    Py_DECREF(event_obj);

    PyObject * assoc_obj = make_assoc_obj(event, num_sta);

    PyList_Append(assocs_list_obj, assoc_obj);
    Py_DECREF(assoc_obj);
    
    free_event(event);
  }

  free(skip_detnums);
  
  free_physics(physics);

  free_sta_list(sta_list);

  free_det_list(det_list);

  /* we are stealing references to the events_list and assocs_list_obj
   * below by using "N" hence we don't need to DECREF them */
  return Py_BuildValue("NN", events_list_obj, assocs_list_obj);
}


/* this method definition and the init-module function are placed at
 * the end of the file since they are simply Python interfacing
 * fluff */

static PyMethodDef csolve_methods[] = {

  {"solve_episode", csolve_solve_episode, METH_VARARGS,
   "Solve a seismic-2d episode."
  },

  /* Sentinel, marking the end of the list */
  {NULL, NULL, 0, NULL
  }
};

PyMODINIT_FUNC initcsolve(void)
{
  PyObject *m;
  
  m = Py_InitModule("csolve", csolve_methods);
  
  if (NULL == m)
    return;
  
  /* do something more with m here if needed, for example add error objects */
}

