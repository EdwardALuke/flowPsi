// implementation file for util.h
#include "util.h"
using std::list ;
using std::cerr ;
using std::cout ;
using std::endl ;
namespace lagrangianP {
  
  entitySet
  remap_entitySet(const entitySet& es, const dMap& remap) {
    entitySet ret ;
    for(entitySet::const_iterator ei=es.begin();ei!=es.end();++ei)
      ret += remap[*ei] ;
    
    return ret ;
  }
  // an overloaded version for Map
  entitySet
  remap_entitySet(const entitySet& es, const Map& remap) {
    entitySet ret ;
    for(entitySet::const_iterator ei=es.begin();ei!=es.end();++ei)
      ret += remap[*ei] ;
    
    return ret ;
  }
  // remap of a sequence
  sequence
  remap_sequence(const sequence& seq, const dMap& remap) {
    sequence ret ;
    for(sequence::const_iterator si=seq.begin();si!=seq.end();++si)
      ret += remap[*si] ;

    return ret ;
  }
  sequence
  remap_sequence(const sequence& seq, const Map& remap) {
    sequence ret ;
    for(sequence::const_iterator si=seq.begin();si!=seq.end();++si)
      ret += remap[*si] ;

    return ret ;
  }

  void
  remap_Map_image(Map& m, const dMap& remap) {
    entitySet domain = m.domain() ;
    for(entitySet::const_iterator ei=domain.begin();
        ei!=domain.end();++ei)
      m[*ei] = remap[m[*ei]] ;
  }
  
  void
  remap_multiMap_image(multiMap& m, const dMap& remap) {
    entitySet domain = m.domain() ;
    for(entitySet::const_iterator ei=domain.begin();
        ei!=domain.end();++ei) {
      int sz = m[*ei].size() ;
      for(int i=0;i<sz;++i)
        m[*ei][i] = remap[m[*ei][i]] ;
    }
  }
  
  void
  remap_Map_image(Loci::storeRepP m, const dMap& remap) {
    Map local_m(m) ;
    remap_Map_image(local_m, remap) ;
  }

  void
  remap_multiMap_image(Loci::storeRepP m, const dMap& remap) {
    Loci::multiMap local_m(m) ;
    remap_multiMap_image(local_m, remap) ;
  }
#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif
  Loci::vector3d<double>
  rotate_point3d(const Loci::vector3d<double>& p,
                 const Loci::vector3d<double>& u, double angle) {
    // first convert angle to radian
    angle = (angle * M_PI) / 180 ;

    Loci::vector3d<double> np ;
    double c = cos(angle) ;
    double s = sin(angle) ;

    np.x = (c + (1-c)*u.x*u.x)*p.x +
      ( (1-c)*u.y*u.x - s*u.z)*p.y + ( (1-c)*u.z*u.x + s*u.y)*p.z ;
  
    np.y = ( (1-c)*u.x*u.y + s*u.z)*p.x +
      (c + (1-c)*u.y*u.y)*p.y + ( (1-c)*u.z*u.y - s*u.x)*p.z ;
  
    np.z = ( (1-c)*u.x*u.z - s*u.y)*p.x +
      ( (1-c)*u.y*u.z + s*u.x)*p.y + (c + (1-c)*u.z*u.z)*p.z ;

    return np ;
  }
  // given a vector in 3d space, computes two orthogonal of its vectors
  // adapted from the metrics.cc (orthogonal_coords class) in Chem/src
  void
  orthogonal_coords(const Loci::vector3d<double>& n,
                    Loci::vector3d<double>& u, Loci::vector3d<double>& v) {
    // Attempt to minimize cancellation error when finding orthogonal vector.
    // Find the coordinate direction which is largest and base orthognality
    // on that direction. (The three following cases are for x largest,
    // then y largest, and finally z largest.
    if(abs(n.x)>abs(n.y) && abs(n.x)>abs(n.z)) {
      if(abs(n.y)>abs(n.z)) {
        u.y = n.y ;
        u.z = -n.x ;
      } else {
        u.y = -n.x ;
        u.z = n.z ;
      }
      u.x = -(u.y * n.y + u.z * n.z) / n.x ;
    } else if(abs(n.y)>abs(n.x) && abs(n.y)>abs(n.z)) {
      if(abs(n.x)>abs(n.z)) {
        u.x = n.x ;
        u.z = -n.y ;
      } else {
        u.x = -n.y ;
        u.z = n.z ;
      }
      u.y = -(u.x * n.x + u.z * n.z) / n.y ;
    } else {
      if(abs(n.x)>abs(n.y)) {
        u.x = n.x ;
        u.y = -n.z ;
      } else {
        u.x = -n.z ;
        u.y = n.y ;
      }
      u.z = -(u.x * n.x + u.y * n.y) / n.z ;
    }
  
    const double usr = 1./sqrt(dot(u,u)) ;
    u *= usr ;  //normalize the vector
    v = cross(u,n) ;  // compute v at last
  }  
  

  double inv_erf(double t) {
    if(t < 0)
      return -inv_erf(-t) ;
    if(t == 0)
      return 0 ;
    if(t >=1.0)
      return 1e300;
    double x = t ;
    const double x2 = x*x ;
    const double A1 = 1.0 ;
    const double A2 = 1./12. ;
    const double A3 = 7./480. ;
    const double A4 = 127./40320. ;
    const double A5 = 4369./5806080. ;
    const double A6 = 34807./182476800. ;
    const double pi = M_PI ;
    const double a = .5*x*(A1+pi*x2*(A2+pi*x2*(A3+pi*x2*(A4+pi*x2*(A5+pi*x2*A6))))) ;
    x = sqrt(pi)*a ;
    // x is an initial guess to the inv_erf, now use newton method to converge
    const int ITMAX=100 ;
    int i ;
    for(i=0;i<ITMAX;++i) {
      double efx = erf(x) ;
      double f = t-efx ;
      if(f < 1e-9*efx)
        break ;
      double fp =  -1.12837916709 * exp(-x*x) ;
      x -= f/fp ;
    }
    if(i == ITMAX)
      cerr << "inv_erf failed convergence, t = " << t << endl ;
    return x ;
  }
  
} // end of namespace lagrangianP
