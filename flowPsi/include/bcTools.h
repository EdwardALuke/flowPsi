#include <Loci.h>
#include "flowTypes.h"

namespace flowPsi {
  inline void cv_inflow(real &pgb, real &Tb, vect3d &ub,
		 real pgi, real Ti, vect3d ui,
		 real pg0, real T0, vect3d u0,
		 vect3d n, real us_n, 
		 real Pambient, real gamma, real Rtilde) {
    const real uit = dot(ui,n)-us_n ; //dot(uin-us,n) ;
    const real u0t = dot(u0,n)-us_n ; //dot(uref-us,n) ;

    const real a02 = gamma*Rtilde*T0 ;
    const real a0 = sqrt(a02) ;
    const real r0 = (pg0+Pambient)/(Rtilde*T0) ;
    const real MinPressure = 0.1*(pg0+Pambient)-Pambient ;

    pgb = max(0.5*(pg0 + pgi - r0*a0*(u0t-uit)),MinPressure) ;

    ub = u0 + ((pgb-pg0)/(r0*a0))*n ;

    const real rhob = r0 + (pgb-pg0)/a02 ;
    
    Tb = (pgb+Pambient)/(Rtilde*rhob) ;
  }
  inline void cv_outflow(real &pgb, real &Tb, vect3d &ub,
		  real pgi, real Ti, vect3d ui,
		  real pg0, real T0, vect3d u0,
		  vect3d n, real us_n, 
		  real Pambient, real gamma, real Rtilde) {
    const real uit = dot(ui,n)-us_n ; //dot(uin-us,n) ;
    const real u0t = dot(u0,n)-us_n ; //dot(uref-us,n) ;

    const real a02 = gamma*Rtilde*T0 ;
    const real a0 = sqrt(a02) ;
    const real r0 = (pg0+Pambient)/(Rtilde*T0) ;
    const real MinPressure = 0.1*(pgi+Pambient)-Pambient ;

    pgb = max(0.5*(pg0 + pgi - r0*a0*(u0t-uit)),MinPressure) ;
  
    ub = ui + ((pgi-pgb)/(r0*a0))*n ;
    const real rhob = r0 + (pgb-pgi)/a02 ;
    Tb = (pgb+Pambient)/(Rtilde*rhob) ;
  }      
}
