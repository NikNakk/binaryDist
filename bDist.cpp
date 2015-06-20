#include <Rcpp.h>
using namespace Rcpp;

//countBits function taken from https://en.wikipedia.org/wiki/Hamming_weight#Efficient_implementation

const uint64_t m1  = 0x5555555555555555; //binary: 0101...
const uint64_t m2  = 0x3333333333333333; //binary: 00110011..
const uint64_t m4  = 0x0f0f0f0f0f0f0f0f; //binary:  4 zeros,  4 ones ...
const uint64_t h01 = 0x0101010101010101; //the sum of 256 to the power of 0,1,2,3...

int countBits(uint64_t x) {
  x -= (x >> 1) & m1;             //put count of each 2 bits into those 2 bits
  x = (x & m2) + ((x >> 2) & m2); //put count of each 4 bits into those 4 bits 
  x = (x + (x >> 4)) & m4;        //put count of each 8 bits into those 8 bits 
  return (x * h01)>>56;  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... 
}

// [[Rcpp::export]]
int countBitsFromRaw(RawVector rv) {
  uint64_t* x = (uint64_t*)RAW(rv);
  return(countBits(*x));
}

// [[Rcpp::export]]
NumericVector bDist(RawMatrix mat) {
  int nr(mat.nrow()), nc(mat.ncol());
  int nw = nr / 8;
  NumericVector res(nc * (nc - 1) / 2);
  // Access the raw data as unsigned 64 bit integers
  uint64_t* data = (uint64_t*)RAW(mat);
  uint64_t a(0);
  // Work through each possible combination of columns (rows in the original integer matrix)
  for (int i = 0; i < nc - 1; i++) {
    for (int j = i + 1; j < nc; j++) {
      uint64_t sx = 0;
      uint64_t so = 0;
      // Work through each 64 bit integer and calculate the sum of (x ^ y) and (x | y)
      for (int k = 0; k < nw; k++) {
        uint64_t o = data[nw * i + k] | data[nw * j + k];
        // If (x | y == 0) then (x ^ y) will also be 0
        if (o) {
          // Use Hamming weight method to calculate number of set bits
          so = so + countBits(o);
          uint64_t x = data[nw * i + k] ^ data[nw * j + k];
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