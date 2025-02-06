//#############################################################################
//#
//# Copyright 2015-2019, Mississippi State University
//#
//# This file is part of the flowPsi computational fluid dynamics solver.
//#
//# The flowPsi solver is free software: you can redistribute it and/or modify
//# it under the terms of the GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The flowPsi solver is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License
//# along with the flowPsi solver.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################

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
