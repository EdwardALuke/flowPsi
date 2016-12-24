#include <Loci.h>
#include <flowTypes.h>

namespace flowPsi {
  inline void inviscidFlux(Loci::Array<real,5> &iflux,
			   real pg, real T, vect3d U,
			   vect3d an, real area,
			   real pambient, real Rt, real gamma,
			   real Us_n) {
    real gm1 = gamma-1 ;
    real rho = (pg+pambient)/(Rt*T) ;
    
    real h0 = Rt*T*gamma/gm1 + 0.5*dot(U,U) ;
    const real e0 = h0-Rt*T ;
    const real uta = dot(U,an) ;
    const real ut = uta-Us_n ; 
    real mdot = area*rho*ut ;
    real p = pg+pambient ;

    iflux[0] = mdot ;
    iflux[1] = mdot*U.x+area*pg*an.x ;
    iflux[2] = mdot*U.y+area*pg*an.y ;
    iflux[3] = mdot*U.z+area*pg*an.z ;
    iflux[4] = mdot*e0+area*p*uta ;
  }
    


  void roeCT_flux(Loci::Array<real,5> &iflux,
		 real pgl, real Tl, vect3d Ul,
		 real pgr, real Tr, vect3d Ur,
		 vect3d an, real area,
		 real pambient, real Rt,real gamma, real Us_n,
		 real Beta) ;
  void roeCT_jacobian(Loci::Mat<real_fj> &Fj,
		      real pgl, real Tl, vect3d Ul,
		      real pgr, real Tr, vect3d Ur,
		      vect3d an, real area,
		      real pambient, real Rt,real gamma, real Us_n,
		      real Beta, int dir) ;

  void hllc_flux(Loci::Array<real,5> &iflux,
		 real pgl, real Tl, vect3d Ul,
		 real pgr, real Tr, vect3d Ur,
		 vect3d an, real area,
		 real pambient, real Rt,real gamma, real Us_n,
		 real Beta) ;
  void hllc_fjm(Mat<real_fj> &fj,
		real pgl, real Tl, vect3d Ul,
		real pgr, real Tr, vect3d Ur,
		vect3d an, real area,
		real pambient, real Rt,real gamma, real Us_n,
		real Beta) ;

  void hllc_fjp(Mat<real_fj> &fj,
		real pgl, real Tl, vect3d Ul,
		real pgr, real Tr, vect3d Ur,
		vect3d an, real area,
		real pambient, real Rt,real gamma, real Us_n,
		real Beta) ;

  inline void inviscidRiemannFlux(Loci::Array<real,5> &iflux,
				  real pgl, real Tl, vect3d Ul,
				  real pgr, real Tr, vect3d Ur,
				  vect3d an, real area,
				  real pambient, real Rt,real gamma, real Us_n,
				  real Beta) {
    hllc_flux(iflux,pgl,Tl,Ul,pgr,Tr,Ur,an,area,pambient,Rt,gamma,Us_n,Beta) ;
    //roeCT_flux(iflux,pgl,Tl,Ul,pgr,Tr,Ur,an,area,pambient,Rt,gamma,Us_n,Beta) ;
  }

  inline void inviscidRiemannFjp(Mat<real_fj> fj,
				 real pgl, real Tl, vect3d Ul,
				 real pgr, real Tr, vect3d Ur,
				 vect3d an, real area,
				 real pambient, real Rt,real gamma, real Us_n,
				 real Beta) {
    roeCT_jacobian(fj,pgl,Tl,Ul,pgr,Tr,Ur,an,area,pambient,Rt,gamma,Us_n,Beta,1) ;
    //hllc_fjp(fj,pgl,Tl,Ul,pgr,Tr,Ur,an,area,pambient,Rt,gamma,Us_n,Beta) ;
  }
  inline void inviscidRiemannFjm(Mat<real_fj> fj,
				 real pgl, real Tl, vect3d Ul,
				 real pgr, real Tr, vect3d Ur,
				 vect3d an, real area,
				 real pambient, real Rt,real gamma, real Us_n,
				 real Beta) {
    roeCT_jacobian(fj,pgl,Tl,Ul,pgr,Tr,Ur,an,area,pambient,Rt,gamma,Us_n,Beta,-1) ;
    //hllc_fjm(fj,pgl,Tl,Ul,pgr,Tr,Ur,an,area,pambient,Rt,gamma,Us_n,Beta) ;
  }
}
