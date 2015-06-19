#include <Rcpp.h>
using namespace Rcpp;

unsigned int countBits(unsigned int x) {
  unsigned int c;
  x = x - ((x >> 1) & 0x55555555);                    // reuse input as temporary
  x = (x & 0x33333333) + ((x >> 2) & 0x33333333);     // temp
  c = ((x + (x >> 4) & 0xF0F0F0F) * 0x1010101) >> 24; // count
  return (c);
}

// [[Rcpp::export]]
NumericVector bDist(RawMatrix mat) {
  int nr(mat.nrow()), nc(mat.ncol());
  NumericVector res(nr * (nr - 1) / 2);
  long a(0);
  for (int i = 0; i < nr - 1; i++) {
    for (int j = i + 1; j < nr; j++) {
      unsigned int sx(0), so(0);
      for (int k = 0; k < nc; k++) {
        unsigned int o = mat(i, k) | mat(j, k);
        if (o) {
          so = so + countBits(o);
          unsigned int x = mat(i, k) ^ mat(j, k);
          if (x) {
            sx = sx + countBits(x);
          }
        }
      }
      res(a++) = (double)sx / so;
    }
  }
  return (res);
}