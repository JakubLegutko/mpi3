#include"Helper.h"

double Helper::getDistanceSQ( double *x, double *y, double *z, int i, int j ) {
   double dx = x[ i ] - x[ j ];
   double dy = y[ i ] - y[ j ];
   double dz = z[ i ] - z[ j ];
   return dx*dx + dy*dy + dz*dz;
}
      
double Helper::getDistance( double *x, double *y, double *z, int i, int j ) {
   return sqrt( getDistanceSQ( x, y, z, i, j ) );
}

double Helper::middle( double *x, int i, int j ) {
  return ( x[ i ] + x[ j ] ) * 0.5;
}
