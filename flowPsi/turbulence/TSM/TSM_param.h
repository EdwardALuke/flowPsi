//#############################################################################
//#
//# Copyright 2016, Mississippi State University
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
#ifndef TSM_PARAM_H
#define TSM_PARAM_H

#include <iostream>
#include <Tools/tools.h>
#include <Loci.h>
#include <flowTypes.h>

namespace flowPsi {

  struct TSM3eq_param {
    real ao,as,anu,abp,anat,ats,cbpcr,cnc,cnatcr ;       // Table 1, Column 1 coefficients in TSM 3 equation model
    real cint,ctscr,crnat,c11,c12,cr,cat,css,ctl ;       // Table 1, Column 2 " "
    real cw1,cw2,cw3,cwr,clam,cmustd,prt,sigmak,sigmaw;  // Table 1, Column 3 " "
    TSM3eq_param() {
      ao     = 4.04   ;
      as     = 2.12   ;
      anu    = 6.75   ;
      abp    = 0.6    ;
      anat   = 200    ;
      ats    = 200    ;
      cbpcr  = 5.0    ; //xw
      cnc    = 0.1    ;
      cnatcr = 1250   ;
      cint   = 0.75   ;
      ctscr  = 1000   ;
      crnat  = 0.02   ;
      c11    = 3.4e-6 ;
      c12    = 1.e-10 ;
      cr     = 0.12   ;
      cat    = 0.035  ;
      css    = 1.5    ;
      ctl    = 4360   ;
      cw1    = 0.44   ;
      cw2    = 0.92   ;
      cw3    = 0.3    ;
      cwr    = 1.5    ;
      clam   = 2.495  ;
      cmustd = 0.09   ;
      prt    = 0.85   ;
      sigmak = 1      ;
      sigmaw = 1.17   ; 
    }
  } ;
}
namespace Loci {
  template<> struct data_schema_traits<flowPsi::TSM3eq_param> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(flowPsi::TSM3eq_param()) ;
      LOCI_INSERT_TYPE(ct,flowPsi::TSM3eq_param,ao) ;
      LOCI_INSERT_TYPE(ct,flowPsi::TSM3eq_param,as) ;
      LOCI_INSERT_TYPE(ct,flowPsi::TSM3eq_param,anu) ;
      LOCI_INSERT_TYPE(ct,flowPsi::TSM3eq_param,abp) ;
      LOCI_INSERT_TYPE(ct,flowPsi::TSM3eq_param,anat) ;
      LOCI_INSERT_TYPE(ct,flowPsi::TSM3eq_param,ats) ;
      LOCI_INSERT_TYPE(ct,flowPsi::TSM3eq_param,cbpcr) ;
      LOCI_INSERT_TYPE(ct,flowPsi::TSM3eq_param,cnc) ;
      LOCI_INSERT_TYPE(ct,flowPsi::TSM3eq_param,cnatcr) ;
      LOCI_INSERT_TYPE(ct,flowPsi::TSM3eq_param,cint) ;
      LOCI_INSERT_TYPE(ct,flowPsi::TSM3eq_param,ctscr) ;
      LOCI_INSERT_TYPE(ct,flowPsi::TSM3eq_param,crnat) ;
      LOCI_INSERT_TYPE(ct,flowPsi::TSM3eq_param,c11) ;
      LOCI_INSERT_TYPE(ct,flowPsi::TSM3eq_param,c12) ;
      LOCI_INSERT_TYPE(ct,flowPsi::TSM3eq_param,cr) ;
      LOCI_INSERT_TYPE(ct,flowPsi::TSM3eq_param,cat) ;
      LOCI_INSERT_TYPE(ct,flowPsi::TSM3eq_param,css) ;
      LOCI_INSERT_TYPE(ct,flowPsi::TSM3eq_param,ctl) ;
      LOCI_INSERT_TYPE(ct,flowPsi::TSM3eq_param,cw1) ;
      LOCI_INSERT_TYPE(ct,flowPsi::TSM3eq_param,cw2) ;
      LOCI_INSERT_TYPE(ct,flowPsi::TSM3eq_param,cw3) ;
      LOCI_INSERT_TYPE(ct,flowPsi::TSM3eq_param,cwr) ;
      LOCI_INSERT_TYPE(ct,flowPsi::TSM3eq_param,clam) ;
      LOCI_INSERT_TYPE(ct,flowPsi::TSM3eq_param,cmustd) ;
      LOCI_INSERT_TYPE(ct,flowPsi::TSM3eq_param,prt) ;
      LOCI_INSERT_TYPE(ct,flowPsi::TSM3eq_param,sigmak) ;
      LOCI_INSERT_TYPE(ct,flowPsi::TSM3eq_param,sigmaw) ;
      return DatatypeP(ct) ;
    }
  } ;
}
#endif
