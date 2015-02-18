#ifndef ROOT_H
#define ROOT_H

#include <iostream>
#include <Tools/tools.h>

namespace flowPsi {
  // Root finder based on Riddler's Method
  template <class T, class FUNC> T find_root(FUNC &func,T x1, T x2, T tol) {
    T fl = func(x1) ;
    T fh = func(x2) ;
    T xl = x1 ;
    T xh = x2 ;
    if((fl < 0 && fh < 0) || (fl > 0 && fh > 0)) {
      std::cerr << "root not bracketed in interval" << std::endl ;
      return 0.5*(x2+x1) ;
    }

    T sol = 0.5*(x2+x1) ;
    
    const int ITMAX = 100 ;
    for(int i=0;i<ITMAX;++i) {
      const T xm = 0.5*(xl+xh) ;
      const T fm = func(xm) ;
      if(fm == 0) return xm ;
      const T sr = sqrt(abs(fm*fm - fl*fh)) ;
      if(sr == 0) return xm ;
      sol = xm + (xm - xl)*(((fl-fh)>0?1:-1)*fm/sr) ;
      T fnew = func(sol) ;
      if(fnew == 0) return sol ;
      if((fnew < 0 && fm > 0) ||(fnew > 0 && fm < 0)) {
        xl = xm ;
        fl = fm ;
        xh = sol ;
        fh = fnew ;
      } else if((fnew < 0 && fl > 0) || (fnew > 0 && fl < 0)) {
        xh = sol ;
        fh = fnew ;
      } else {
        xl = sol ;
        fl = fnew ;
      }
      if(abs(xh-xl) <= tol)
        return sol ;
    }
    std::cerr << "find_root failed to converge" << std::endl ;
    return sol ;
  }
}
#endif
