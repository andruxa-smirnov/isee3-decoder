#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "timeformat.h"

void hms(int *hours, int *minutes, double *seconds,double t){
  assert(hours != NULL && minutes != NULL && seconds != NULL);

  *hours = t / 3600.;
  t -= *hours * 3600;

  *minutes = t / 60.;
  t -= *minutes * 60;

  *seconds = t;
}

// Format a time in seconds as hh:mm:ss.ss
// Note: not threadsafe (returns pointer to internal static)
char *format_hms(double t){
  static char output[1000];
  int hours,minutes,count,remain;
  double seconds;
  char *cp;

  hms(&hours,&minutes,&seconds,t);
  
  cp = output;
  remain = sizeof(output);
  if(hours > 0){
    count = snprintf(cp,remain,"%d:",hours);
    cp += count;
    remain -= count;
  }
  count = snprintf(cp,remain,"%02d:",minutes);
  cp += count;
  remain -= count;
  count = snprintf(cp,remain,"%02.3lf",seconds);
  cp += count;
  remain -= count;

  return output;

}
