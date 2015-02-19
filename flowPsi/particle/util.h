// this file contains various utility functions used
// in the particle tracking code
#ifndef PARTICLE_UTIL_H
#define PARTICLE_UTIL_H

#include "particle_config.h"

#include <Loci.h>
#include <list>
#include <string>
#include <sstream>

#ifdef NO_CSTDLIB
#include <stdlib.h>
#include <math.h>
#else
#include <cstdlib>
#include <cmath>
#endif

#include <sys/time.h>

#ifdef USE_PAPI
#include <papi.h>
#endif

namespace lagrangianP {

  // first we define a stop watch for timing purpose
  // below is a basic stop watch with accumulate elapsed time count
  // this one is not thread safe

  // operator defined for timing variables
  inline double
  operator-(const timeval& t1, const timeval& t2) {
    double dt1 = t1.tv_sec + t1.tv_usec*1e-6 ;
    double dt2 = t2.tv_sec + t2.tv_usec*1e-6 ;
    return dt1-dt2 ;
  }

  class basic_timer {
    bool ticking ;
    double accumu_time ;
    double step_time ;
    timeval start_t ;
    timeval stop_t ;
    basic_timer(const basic_timer& t) ;
    basic_timer& operator=(const basic_timer& t) ;
  public:
    basic_timer():ticking(false),accumu_time(0),step_time(0) {}
    void reset() {
      ticking = false ;
      accumu_time = 0 ;
      step_time = 0 ;
    }
    void start() {
      if(!ticking) {
        ticking = true ;
        step_time = 0 ;
        gettimeofday(&start_t,NULL) ;
      }
    }
    void stop() {
      if(ticking) {
        gettimeofday(&stop_t,NULL) ;
        ticking = false ;
        step_time = stop_t - start_t ;
        accumu_time += step_time ;
      }
    }
    // glance the clock without stoping it
    // returns a value since the start of the timer
    double glance() const {
      if(ticking) {
        timeval cur_t ;
        gettimeofday(&cur_t,NULL) ;
        return cur_t - start_t ;
      } else {
        return 0 ;
      }
    }
    bool is_ticking() const {
      return ticking ;
    }
    double read() const {
      return step_time ;
    }
    double read_accumu() const {
      return accumu_time ;
    }
  } ;

#ifdef USE_PAPI
  // this is a PAPI timer, it should use the most accurate timer
  // available on the platform according to PAPI documentation
  class papi_timer {
    bool ticking ;
    // in order to be compatible with the basic_timer,
    // the returned timing value are also in seconds.
    double accumu_time ;
    double step_time ;
    long_long start_t ;
    long_long stop_t ;
    papi_timer(const papi_timer& t) ;
    papi_timer& operator=(const papi_timer& t) ;
  public:
    papi_timer():ticking(false),accumu_time(0),step_time(0) {}
    void reset() {
      ticking = false ;
      accumu_time = 0 ;
      step_time = 0 ;
    }
    void start() {
      if(!ticking) {
        ticking = true ;
        step_time = 0 ;
        start_t = PAPI_get_real_usec() ;
      }
    }
    void stop() {
      if(ticking) {
        stop_t = PAPI_get_real_usec() ;
        ticking = false ;
        step_time = (stop_t - start_t) * 1e-6 ;
        accumu_time += step_time ;
      }
    }
    // glance the clock without stoping it
    // returns a value since the start of the timer
    double glance() const {
      if(ticking) {
        long_long cur_t = PAPI_get_real_usec() ;
        return (cur_t - start_t) * 1e-6 ;
      } else {
        return 0 ;
      }
    }
    bool is_ticking() const {
      return ticking ;
    }
    double read() const {
      return step_time ;
    }
    double read_accumu() const {
      return accumu_time ;
    }
  } ;  
#endif
  
  // normalization of 3d vector
  inline void
  normalize(Loci::vector3d<double>& v) {
    if( (fabs(v.x) <= 0.0) &&
        (fabs(v.y) <= 0.0) &&
        (fabs(v.z) <= 0.0)) {
      //cerr << "WARNING: normalizing zero vector, nothing done!" << endl ;
      return ;
    }
    double t = sqrt(v.x*v.x + v.y*v.y + v.z*v.z) ;
    v.x /= t ;
    v.y /= t ;
    v.z /= t ;
  }

  // a small utility function to test if a position is
  // within an axis-aligned box
  inline bool
  position_in_box(const Loci::vector3d<double>& pos,
                  const Loci::Array<double,6>& box) {
    if(pos.x < box[0])
      return false ;
    if(pos.x > box[1])
      return false ;
    if(pos.y < box[2])
      return false ;
    if(pos.y > box[3])
      return false ;
    if(pos.z < box[4])
      return false ;
    if(pos.z > box[5])
      return false ;

    return true ;
  }

  // this function remaps the image of a Map
  // but it does not remap the domain of a Map
  void
  remap_Map_image(Map& m, const dMap& remap) ;
  // this function remaps the image of a multiMap
  void
  remap_multiMap_image(Loci::multiMap& m, const dMap& remap) ;
  // these are for the storeRepP version
  void
  remap_Map_image(Loci::storeRepP m, const dMap& remap) ;
  void
  remap_multiMap_image(Loci::storeRepP m, const dMap& remap) ;

  // a little utility that converts a number to a string
  template<typename T> inline std::string
  num2str(T n) {
    std::stringstream ss ;
    ss << n ;
    return ss.str() ;
  }
  
  // a char* -> int
  inline int
  chars2int(char* s) {
    std::stringstream ss ;
    ss << s ;
    int ret ;
    ss >> ret ;
    
    return ret ;
  }
  
  // a char* -> long
  inline long
  chars2long(char* s) {
    std::stringstream ss ;
    ss << s ;
    long ret ;
    ss >> ret ;
    
    return ret ;
  }

  // here is a utility function that rotates a point
  // around an arbitrary axis for an angle (in degree)
  const double PI =
  3.141592653589793238462643383279502884197169399375105820974944592308 ;

  Loci::vector3d<double>
  rotate_point3d(const Loci::vector3d<double>& p,
                 const Loci::vector3d<double>& u, double angle) ;

  // given a vector in 3d space, computes two orthogonal of its vectors
  // adapted from the metrics.cc (orthogonal_coords class) in Chem/src
  void
  orthogonal_coords(const Loci::vector3d<double>& n,
                    Loci::vector3d<double>& u, Loci::vector3d<double>& v) ;

  // here is a small utility function that generates random
  // integers and doubles in the range [0,N)
  template <typename T> T
  get_random_number(T N) {
    double d = Loci::random() ;
    return static_cast<T>(d*N) ;
  }

  double inv_erf(double t) ;

  // Gaussian Distribution centered at zero
  inline double GSD(double d = 0, double sigma = 0.25) {
    double r = Loci::random() ; 
    double z = inv_erf(2.*r-1.0) ;
    return (sigma*sqrt(2.)*z+d) ;
  }
  // Log normal distribution centered about d
  inline double LND(double d, double sigma = 0.25) {
    double r = Loci::random() ;
    double z = inv_erf(2.*r-1.0) ;
    return exp(z*sigma*sqrt(2.)+log(d)) ;
  }

} // end of namespace lagrangianP

#endif
