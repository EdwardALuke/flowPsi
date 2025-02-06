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

#ifndef SA_PARAM_H
#define SA_PARAM_H

namespace flowPsi {

  struct Spa_All_param {
    real cb1,k_coeff,cv1,cb2,cw2,cw3,sigma ; //coefficients in Spalart_Allmaras
    //                                         //model

    Spa_All_param()   {
      cb1=0.1355; //production coeffient
      k_coeff=0.41; // coeffient for constructing vorticity \tilde{S}
      cv1=7.1;  //coeffient for constructing f_v1 which is the ratio of
      //        //eddy viscoisty to molecular viscosity
      cb2=0.622; // coeffient in diffusion term
      cw2=0.3; // coeffient in destruction term
      cw3=2.0; //coeffienct in destruction term
      sigma=2./3.; //coeffienct in diffusion term
    }
  } ;
}

namespace Loci {
  template<> struct data_schema_traits<flowPsi::Spa_All_param> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(flowPsi::Spa_All_param()) ;
      LOCI_INSERT_TYPE(ct,flowPsi::Spa_All_param,cb1) ;
      LOCI_INSERT_TYPE(ct,flowPsi::Spa_All_param,k_coeff) ;
      LOCI_INSERT_TYPE(ct,flowPsi::Spa_All_param,cv1) ;
      LOCI_INSERT_TYPE(ct,flowPsi::Spa_All_param,cb2) ;
      LOCI_INSERT_TYPE(ct,flowPsi::Spa_All_param,cw2) ;
      LOCI_INSERT_TYPE(ct,flowPsi::Spa_All_param,cw3) ;
      LOCI_INSERT_TYPE(ct,flowPsi::Spa_All_param,sigma) ;
      return DatatypeP(ct) ;
    }
  } ;
}

#endif
