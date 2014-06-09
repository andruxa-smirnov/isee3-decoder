// Generate complex sinusoid for testing
#define _GNU_SOURCE 1
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>
#include <fcntl.h>
#include <errno.h>
#include <stdint.h>
#include <fftw3.h>
#include <complex.h>

// Complex pair of 16-bit, 2's complement receiver samples
struct sample {
  int16_t i;    // In-phase channel
  int16_t q;    // Quadrature channel
};


double Carrier = 2000.0; // Carrier freq, Hz
double Samprate = 32768;  // Sample rate, Hz
double Amplitude = 20000; // Digital peak amplitude, units
double Startphase = 0.;   // Initial carrier phase, radians
int Nsamp;                // Number of samples to emit

int main(int argc,char *argv[]){
  int i;
  double cphase;   // Starting carrier phase
  double cstep;           // Phase increment in radians/sample
  double complex v;
  struct sample s;

  cphase = Startphase;
  cstep = 2 * M_PI * Carrier / Samprate;
  Nsamp = 10 * Samprate;

  fprintf(stderr,"carrier %lf Hz, sample rate %lf Hz, amplitude %lf, phaseinc %lg rad/samp\n",
	  Carrier,Samprate,Amplitude,cstep);
  for(i=0;i<Nsamp;i++){
    v = Amplitude * cexp(I*cphase);
    s.i = creal(v);
    s.q = cimag(v);
    fwrite(&s,sizeof(s),1,stdout);

    // increment phase and reduce mod 2 pi
    cphase += cstep;
    if(cphase >= 2*M_PI)
      cphase -= 2*M_PI;
  }
  exit(0);
}
