#ifndef TURB_PARAM_H
#define TURB_PARAM_H

#include <iostream>
#include "flowTypes.h"
#include <Tools/tools.h>
#include <Loci.h>

namespace flowPsi {
  struct rke_param {
    real a1,ce1,ce2,sigmak,sigmae,ctau,amu,ae,aet,cs ;
    rke_param() {
      a1 = 1.25 ;
      ce1 = 1.44 ;
      ce2 = 1.92 ;
      sigmak = 1.0 ;
      sigmae = 1.3 ;
      ctau = 1.414213562 ; //sqrt(2.0)
      amu = 0.01 ;
      ae = 0.3 ;
      aet = 0.15 ;
      cs = 0.05 ;
    }
  } ;

}
namespace Loci {
  template<> struct data_schema_traits<flowPsi::rke_param> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(flowPsi::rke_param()) ;
      LOCI_INSERT_TYPE(ct,flowPsi::rke_param,a1) ;
      LOCI_INSERT_TYPE(ct,flowPsi::rke_param,ce1) ;
      LOCI_INSERT_TYPE(ct,flowPsi::rke_param,ce2) ;
      LOCI_INSERT_TYPE(ct,flowPsi::rke_param,sigmak) ;
      LOCI_INSERT_TYPE(ct,flowPsi::rke_param,sigmae) ;
      LOCI_INSERT_TYPE(ct,flowPsi::rke_param,ctau) ;
      LOCI_INSERT_TYPE(ct,flowPsi::rke_param,amu) ;
      LOCI_INSERT_TYPE(ct,flowPsi::rke_param,ae) ;
      LOCI_INSERT_TYPE(ct,flowPsi::rke_param,aet) ;
      LOCI_INSERT_TYPE(ct,flowPsi::rke_param,cs) ;
      return DatatypeP(ct) ;
    }
  } ;
}

#endif
