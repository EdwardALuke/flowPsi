#include <Loci.h>
#include "flowTypes.h"

namespace flowPsi {
  inline void cv_inflow(real &pgb, real &Tb, vect3d &ub,
			real pgi, real Ti, vect3d ui,
			real pg0, real T0, vect3d u0,
			vect3d n, real us_n, 
			real Pambient, real gamma, real Rtilde,
			real Beta) {
    const real uit = dot(ui,n)-us_n ; //dot(uin-us,n) ;
    const real u0t = dot(u0,n)-us_n ; //dot(uref-us,n) ;

    const real a02 = gamma*Rtilde*T0 ;
    const real a0 = sqrt(a02) ;
    const real r0 = (pg0+Pambient)/(Rtilde*T0) ;
    const real MinPressure = 0.1*(pg0+Pambient)-Pambient ;
    const real BetaM = 0.5*(1.0-Beta) ;
    const real uitBeta  = uit*(1.-Beta) ;
    const real sigma2 = uitBeta*uitBeta + 4.*Beta*a0*a0 ;
    const real sigma = sqrt(sigma2) ;

    pgb = 0.5*(pg0 + pgi - 0.5*r0*sigma*(u0t-uit)) 
      -((pgi-pg0) + r0*(uit-u0t)*uit*BetaM)*uit*BetaM/(2.*sigma) ;
    pgb = max(pgb,MinPressure) ;

    ub = u0 + 2.*n*(pgb-pg0)/(r0*(sigma+uit*BetaM)) ;
    real rhob = r0 + (pgb-pg0)*(sigma-uit*BetaM)/
      (Beta*a0*a0*(sigma+uit*BetaM)) ;

    Tb = (pgb+Pambient)/(Rtilde*rhob) ;
  }
  inline void cv_outflow(real &pgb, real &Tb, vect3d &ub,
			 real pgi, real Ti, vect3d ui,
			 real pg0, real T0, vect3d u0,
			 vect3d n, real us_n, 
			 real Pambient, real gamma, real Rtilde,
			 real Beta) {
    const real uit = dot(ui,n)-us_n ; //dot(uin-us,n) ;
    const real u0t = dot(u0,n)-us_n ; //dot(uref-us,n) ;

    const real a02 = gamma*Rtilde*T0 ;
    const real a0 = sqrt(a02) ;
    const real r0 = (pg0+Pambient)/(Rtilde*T0) ;
    const real MinPressure = 0.1*(pgi+Pambient)-Pambient ;
    const real BetaM = 0.5*(1.0-Beta) ;

    const real uitBeta  = uit*(1.-Beta) ;
    const real sigma2 = uitBeta*uitBeta + 4.*Beta*a0*a0 ;
    const real sigma = sqrt(sigma2) ;

    pgb = .5*(pgi+pg0-0.5*r0*sigma*(u0t-uit)) 
      +((pgi-pg0) + r0*(uit-u0t)*uit*BetaM)*uit*BetaM/(2.*sigma) ;
    pgb = max(pgb,MinPressure) ;
    ub = ui + 2.*n*(pgi-pgb)/(r0*(sigma+uit*BetaM)) ;

    const real rhob = r0 + (pgb-pgi)*(sigma-uit*BetaM)/
      (Beta*a0*a0*(sigma+uit*BetaM)) ;

    Tb = (pgb+Pambient)/(Rtilde*rhob) ;
  }      
}
