#include <iostream>
using namespace std ;

int main() {
  int ni = 101 ;
  int nj = 11 ;
  int nk = 2 ;

  cout << 1 << endl ;
  cout << ni << ' ' << nj << ' ' << nk << endl ;
  double dx = 0.01 ;
  for(int k=0;k<nk;++k) 
    for(int j=0;j<nj;++j)
      for(int i=0;i<ni;++i) {
	cout << dx*double(i) << endl ;
      }
  for(int k=0;k<nk;++k) 
    for(int j=0;j<nj;++j)
      for(int i=0;i<ni;++i) {
	cout << dx*double(j) << endl ;
      }
  for(int k=0;k<nk;++k) 
    for(int j=0;j<nj;++j)
      for(int i=0;i<ni;++i) {
	cout << dx*double(k)-0.5*dx << endl ;
      }
}
