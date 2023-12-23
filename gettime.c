#include  <sys/time.h>
#include <stdio.h>
#pragma weak gettime_=gettime
double gettime(void)
{
   struct timeval tv;
   double elapsed;
   gettimeofday(&tv, NULL);
   elapsed = ((double) tv.tv_sec) + 1.0e-6*((double) tv.tv_usec);
   return elapsed;
}

