#ifndef HELPER_H_
#define HELPER_H_

#include<math.h>

class Helper {
   public:
      static double getDistanceSQ( double *x, double *y, double *z, int i, int j );
      
      static double getDistance( double *x, double *y, double *z, int i, int j );
      
      static double middle( double *x, int i, int j );

};

#endif

