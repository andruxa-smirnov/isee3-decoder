// Time conversion: convert floating point seconds into
// dd:hh:mm:ss.sss and pretty-print into a string

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "timeformat.h"

void hms(int *days,int *hours, int *minutes, double *seconds,double t){
  assert(days != NULL && hours != NULL && minutes != NULL && seconds != NULL);

  *days = t / 86400.;
  t -= *days * 86400;

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
  int days,hours,minutes,count,remain;
  double seconds;
  char *cp;

  hms(&days,&hours,&minutes,&seconds,t);
  
  cp = output;
  remain = sizeof(output);
  if(days > 0){
    count = snprintf(cp,remain,"%d:",days);
    cp += count;
    remain -= count;
  }

  if(days > 0 || hours > 0){
    count = snprintf(cp,remain,"%02d:",hours);
    cp += count;
    remain -= count;
  }
  count = snprintf(cp,remain,"%02d:",minutes);
  cp += count;
  remain -= count;
  if(seconds < 10.0){
    // print the leading zero; %f doesn't do this
    count = snprintf(cp,remain,"0");
    cp += count;
    remain -= count;
  }
  count = snprintf(cp,remain,"%.3lf",seconds);
  cp += count;
  remain -= count;

  return output;

}
