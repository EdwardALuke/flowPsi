#ifndef GRID_COMPONENT_H
#define GRID_COMPONENT_H
#include "Loci.h"
#include "flowTypes.h"
#include <vector>
namespace flowPsi {
  struct componentXform {
    vect3d cg, new_cg ;
    Quaternion q ;
    componentXform() {
      cg = vect3d(0,0,0) ;
      new_cg = vect3d(0,0,0) ;
      q.x = 0 ;
      q.y = 0 ;
      q.z = 0 ;
      q.w = 1 ;
    }
    vect3d applyXform(vect3d pos) const {
      return (q*(pos-cg)) + new_cg ;
    }
    vect3d applyRotation(vect3d vec) const {
      return q*vec ;
    }
  } ;

  inline std::ostream &operator<<(std::ostream &s, const componentXform x) {
    Loci::Abort() ;
    return s ;
  }
  inline std::istream &operator>>(std::istream &s, const componentXform x) {
    Loci::Abort() ;
    return s ;
  }

  class geometry_type : public Loci::CPTR_type {
  public:
    virtual geometry_type *applyXform(componentXform xform) const = 0 ;
    virtual bool inGeometry(vect3d pt) const = 0 ;
    virtual real distToSurface(vect3d pt) const = 0 ;
  } ;

  class cylinder_type : public geometry_type {
    vect3d pt1,pt2 ;
    double radius ;
    double radius2 ; // radius squared
    double rlen2 ; // reciprocal of len(pt1-pt2)^2
  public:
    cylinder_type() {
      pt1=vect3d(0,0,0) ;
      pt2=vect3d(0,0,1) ;
      radius = 0 ;
      radius2 = 0 ;
      rlen2 = 1;
    }
    cylinder_type(vect3d p1, vect3d p2, double r) : pt1(p1),pt2(p2),radius(r)
    { rlen2 = 1./dot(pt2-pt1,pt2-pt1) ; radius2 = r*r; }
    geometry_type *applyXform(componentXform xform) const ;
    bool inGeometry(vect3d pt) const ;
    real distToSurface(vect3d pt) const ;
  } ;

  class planelist_type : public geometry_type {
    std::vector<vect3d> pts ;
    std::vector<vect3d> normals ;
  public:
    planelist_type() {
    }
    planelist_type(std::vector<vect3d> &ps, std::vector<vect3d> &ns) :
      pts(ps),normals(ns) {} 
    geometry_type *applyXform(componentXform xform) const ;
    bool inGeometry(vect3d pt) const ;
    real distToSurface(vect3d pt) const ;
  } ;

  class sphere_type : public geometry_type {
    vect3d center ;
    double radius ;
    double radius2 ; // radius squared
  public:
    sphere_type() {
      center=vect3d(0,0,0) ;
      radius = 0 ;
      radius2 = 0 ;
    }
    sphere_type(vect3d c,double r) : center(c),radius(r)
    { radius2 = r*r; }
    geometry_type *applyXform(componentXform xform) const ;
    bool inGeometry(vect3d pt) const ;
    real distToSurface(vect3d pt) const ;
  } ;

  class revolution_type : public geometry_type {
    vect3d p1,p2 ;
    double rlen2 ; // reciprocal of len(pt1-pt2)^2
    std::vector<pair<double,double> > radius_pairs ;
  public:
    revolution_type() {
      p1=vect3d(0,0,0) ;
      p2=vect3d(0,0,0) ;
    }
    revolution_type(vect3d p1i, vect3d p2i,
                    std::vector<pair<double,double> > rp ) : p1(p1i),p2(p2i),
                                                        radius_pairs(rp)
    { std::sort(radius_pairs.begin(),radius_pairs.end()) ; 
      rlen2 = 1./dot(p2-p1,p2-p1) ; }
    geometry_type *applyXform(componentXform xform) const ;
    bool inGeometry(vect3d pt) const ;
    real distToSurface(vect3d pt) const ;
  } ;
  // Find the time interval that we are splining
  int findt(const std::vector<double> &t, double tval) ;
  // Compute spline derivatives for hermite spline
  void splineD(std::vector<double> &xp, const std::vector<double> &x,
               const std::vector<double> &t) ;
  // cubic spline
  double spline(int ind,double tval,const std::vector<double> &t,
                const std::vector<double> &x,
                const std::vector<double> &xp) ;
  
  struct motionSplines {
    std::vector<double> t,x,y,z,q0,q1,q2,q3 ;
    std::vector<double> xp,yp,zp,q0p,q1p,q2p,q3p ;
    void initialize(std::vector<double> ti,std::vector<double> xi,
                    std::vector<double> yi,std::vector<double> zi,
                    std::vector<double> q0i,std::vector<double> q1i,
                    std::vector<double> q2i,std::vector<double> q3i) {
      t = ti ;
      x = xi ;
      splineD(xp,x,t) ;
      y = yi ;
      splineD(yp,y,t) ;
      z = zi ;
      splineD(zp,z,t) ;
      q0 = q0i ;
      splineD(q0p,q0,t) ;
      q1 = q1i ;
      splineD(q1p,q1,t) ;
      q2 = q2i ;
      splineD(q2p,q2,t) ;
      q3 = q3i ;
      splineD(q3p,q3,t) ;
    }
    void getMotion(vect3d &cg,Quaternion &q,double tv) const {
      int ind = findt(t,tv) ;
      cg.x = spline(ind,tv,t,x,xp) ;
      cg.y = spline(ind,tv,t,y,yp) ;
      cg.z = spline(ind,tv,t,z,zp) ;
      q.x = spline(ind,tv,t,q0,q0p) ;
      q.y = spline(ind,tv,t,q1,q1p) ;
      q.z = spline(ind,tv,t,q2,q2p) ;
      q.w = spline(ind,tv,t,q3,q3p) ;
      q.Normalize() ;
    }
  } ;

      
}
namespace Loci {
  template<> struct data_schema_traits<flowPsi::componentXform> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(flowPsi::componentXform()) ;
      LOCI_INSERT_TYPE(ct,flowPsi::componentXform,cg) ;
      LOCI_INSERT_TYPE(ct,flowPsi::componentXform,new_cg) ;
      LOCI_INSERT_TYPE(ct,flowPsi::componentXform,q) ;

      return DatatypeP(ct) ;
    }
  } ;
}
#endif
