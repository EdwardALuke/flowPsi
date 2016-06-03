#ifndef MTR_PARAM_H
#define MTR_PARAM_H

#include <iostream>
#include "flowTypes.h"
#include <Tools/tools.h>
#include <Loci.h>

namespace flowPsi {

  struct MTR1eq_param {
    real aplus,kappa,c1,c2,c3,sigma ; //coefficients in Menter 1 Eq model
    MTR1eq_param() {
      aplus = 13. ;
      kappa = 0.41 ;
      c1 = 0.144 ;
      c2 = 1.86 ;
      c3 = 7. ;
      sigma = 1. ;
    }
  } ;
}

namespace Loci {
  template<> struct data_schema_traits<flowPsi::MTR1eq_param> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(flowPsi::MTR1eq_param()) ;
      LOCI_INSERT_TYPE(ct,flowPsi::MTR1eq_param,aplus) ;
      LOCI_INSERT_TYPE(ct,flowPsi::MTR1eq_param,kappa) ;
      LOCI_INSERT_TYPE(ct,flowPsi::MTR1eq_param,c1) ;
      LOCI_INSERT_TYPE(ct,flowPsi::MTR1eq_param,c2) ;
      LOCI_INSERT_TYPE(ct,flowPsi::MTR1eq_param,c3) ;
      LOCI_INSERT_TYPE(ct,flowPsi::MTR1eq_param,sigma) ;
      return DatatypeP(ct) ;
    }
  } ;
}
#endif
