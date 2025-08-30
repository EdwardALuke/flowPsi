#ifndef PLANE_ROTATE_H
#define PLANE_ROTATE_H

#include <vector>

namespace LGSI {

  template <typename REAL> 
  inline void GeneratePlaneRotation(REAL dx, REAL dy, REAL &cs, REAL &sn) {
    if (dy == 0.0) {
      cs = 1.0;
      sn = 0.0;
    } else if (abs(dy) > abs(dx)) {
      REAL temp = dx / dy;
      sn = 1.0 / sqrt( 1.0 + temp*temp );
      cs = temp * sn;
    } else {
      REAL temp = dy / dx;
      cs = 1.0 / sqrt( 1.0 + temp*temp );
      sn = temp * cs;
    }
  }

  template <typename REAL> inline void
  backSolveH(std::vector<REAL> &y,
             const std::vector<std::vector<REAL> > &h,
             const std::vector<REAL> &s,
             int k) {
    y[k] = s[k]/h[k][k] ;
    for(int i=k-1; i>=0; i--) {
      double tmp = s[i] ;
      for(int j=k; j>=i+1; j--) {
	tmp -= h[j][i]*y[j];
      }
      y[i] = tmp/h[i][i];
    } 
  }

}
#endif
