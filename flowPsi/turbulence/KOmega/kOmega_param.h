//#############################################################################
//#
//# Copyright 2015, Mississippi State University
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

#ifndef TURB_PARAM_H
#define TURB_PARAM_H

namespace flowPsi {
  const real EPSILON = 1e-30 ;
  const real pi = M_PI ;

  struct sst_param {
    real sigmak, sigmae, beta, gama ;
  } ;

  struct sst1_param {
    real sigmak, sigmae, beta, betas, kappa, a1 ;
    sst1_param() {
      //sigmak = 0.5 ;
      sigmak = 0.85 ;
      sigmae = 0.5 ;
      beta = 0.075 ;
      betas = 0.09 ;
      kappa = 0.41 ;
      a1 = 0.31 ;
    }
  } ;

  struct sst2_param {
    real sigmak, sigmae, beta, betas, kappa ;
    sst2_param() {
      sigmak = 1.0 ;
      sigmae = 0.856 ;
      beta = 0.0828 ;
      betas = 0.09 ;
      kappa = 0.41 ;
    }
  } ;


}

namespace Loci {

  template<> struct data_schema_traits<flowPsi::sst_param> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(flowPsi::sst_param()) ;
      LOCI_INSERT_TYPE(ct,flowPsi::sst_param,sigmak) ;
      LOCI_INSERT_TYPE(ct,flowPsi::sst_param,sigmae) ;
      LOCI_INSERT_TYPE(ct,flowPsi::sst_param,beta) ;
      LOCI_INSERT_TYPE(ct,flowPsi::sst_param,gama) ;
      return DatatypeP(ct) ;
    }
  } ;


  template<> struct data_schema_traits<flowPsi::sst1_param> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(flowPsi::sst1_param()) ;
      LOCI_INSERT_TYPE(ct,flowPsi::sst1_param,sigmak) ;
      LOCI_INSERT_TYPE(ct,flowPsi::sst1_param,sigmae) ;
      LOCI_INSERT_TYPE(ct,flowPsi::sst1_param,beta) ;
      LOCI_INSERT_TYPE(ct,flowPsi::sst1_param,betas) ;
      LOCI_INSERT_TYPE(ct,flowPsi::sst1_param,kappa) ;
      LOCI_INSERT_TYPE(ct,flowPsi::sst1_param,a1) ;
      return DatatypeP(ct) ;
    }
  } ;

  template<> struct data_schema_traits<flowPsi::sst2_param> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(flowPsi::sst2_param()) ;
      LOCI_INSERT_TYPE(ct,flowPsi::sst2_param,sigmae) ;
      LOCI_INSERT_TYPE(ct,flowPsi::sst2_param,beta) ;
      LOCI_INSERT_TYPE(ct,flowPsi::sst2_param,betas) ;
      LOCI_INSERT_TYPE(ct,flowPsi::sst2_param,kappa) ;
      return DatatypeP(ct) ;
    }
  } ;

}
#endif
