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


