//#############################################################################
//#
//# Copyright 2014-2019, Mississippi State University
//#
//# The GridMover is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The GridMover software is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the GridMover software. If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################
#ifndef NODEDATA_H
#define NODEDATA_H

#include "rotor.h"

namespace gridMotion {

  struct nodeDataBase {
    vector3d<realF> pos ;
    realF weight ;
    int order ;
  } ;

  struct NodeData {
    vector3d<real> pos;
    vector3d<realF> disp;
    Rotor            rot;
    realF           weight ;
  };

  inline std::ostream & operator<<(std::ostream & lhs, const NodeData & rhs) {
    lhs << rhs.pos       << std::endl;
    lhs << rhs.disp      << std::endl;
    lhs << rhs.rot.alpha << std::endl;
    lhs << rhs.rot.beta  << std::endl;
    lhs << rhs.weight    << std::endl;
    return lhs;
  }

  inline std::istream & operator>>(std::istream & lhs, NodeData & rhs) {
    lhs >> rhs.pos >> rhs.disp >> rhs.rot.alpha >> rhs.rot.beta >> rhs.weight;
    return lhs;
  }

  template<class T>
  struct Append {
    void operator()(T & lhs, const T & rhs) {
      lhs.insert(lhs.end(),rhs.begin(),rhs.end());
    }
  };
  
  
}

#endif // NODEDATA_H
