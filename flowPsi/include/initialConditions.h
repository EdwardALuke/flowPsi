#ifndef INITIALCONDITIONS_HEADER
#define INITIALCONDITIONS_HEADER

#include <Loci.h>
#include "flowTypes.h"
#include <vector>

namespace flowPsi {
  class geomTest : public Loci::CPTR_type {
  public:
    virtual bool inGeomPt(vect3d pt) const = 0 ;
  } ;
  Loci::CPTR<geomTest> geomTestFactory(string name, const options_list ol) ;

  struct ICstate_info {
    string name ;
    Loci::options_list state_info ;
  } ;
  struct ICgeom_info {
    Loci::CPTR<geomTest> geomTestFunc;
    int id ; // State for this geometry
  } ;
  
}



#endif
