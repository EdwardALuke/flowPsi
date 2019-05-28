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
#ifndef INTERPOLATE_FILE_H
#define INTERPOLATE_FILE_H
#include <Loci.h>
#include "flowTypes.h"
#include <string>
#include <vector>
#include <interpolate.h>

namespace flowPsi {

  using Loci::get_stencil ;
  using Loci::stencil_weights ;
  extern void broadcast_storeRep(Loci::storeRepP rep) ;
  
  struct stencil_info {
    std::vector<Loci::Array<real,4> > weights ;
    std::vector<Loci::Array<int ,4> > stencils ;
    std::vector<int> send_info, req_sizes, snd_sizes ;
    Loci::storeRepP slookup ;
  } ;
  
  using Loci::sendStencilData ;

}
#endif


