#ifndef WALL_FUNC_PARAM_H
#define WALL_FUNC_PARAM_H

#include <iostream>
#include "flowTypes.h"
#include <Tools/tools.h>
#include <Loci.h>

namespace flowPsi {
  struct wall_func_param {
    real kappa,B,E ; //coeffients in wall law function
    real Cmu ; // assumed turbulence coefficient Cmu
    wall_func_param() ;
    std::istream &Input(std::istream &s) ;
    std::ostream &Print(std::ostream &s) const ;
  } ;

  inline std::ostream & operator<<(std::ostream &s, 
				   const wall_func_param &wall)
    {return wall.Print(s) ; }
  inline std::istream & operator>>(std::istream &s, 
				   wall_func_param &wall) 
    {return wall.Input(s) ; }
}

namespace Loci {
  template<> struct data_schema_traits<flowPsi::wall_func_param> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(flowPsi::wall_func_param()) ;
      LOCI_INSERT_TYPE(ct,flowPsi::wall_func_param,kappa) ;
      LOCI_INSERT_TYPE(ct,flowPsi::wall_func_param,B) ;
      LOCI_INSERT_TYPE(ct,flowPsi::wall_func_param,E) ;
      return DatatypeP(ct) ;
    }
  } ;
}

#endif
