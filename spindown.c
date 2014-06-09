// Complex frequency shifter
// Quick hack from icedemod.c
// 2 June 2014, Phil Karn, KA9Q
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
#include <complex.h> // Must be before fftw3.h so it understands native complex double
#include <fftw3.h>
#include <getopt.h>
#include <locale.h>
#include "code.h"

#undef I
#define Q _Complex_I  // Less likely to cause error than the default I

// Complex pair of 16-bit, 2's complement receiver samples
struct sample {
  int16_t i;    // In-phase channel
  int16_t q;    // Quadrature channel
};

double Samprate = 250000; // Sample rate, Hz
int Fftsize = 131072;              // Size of FFT for carrier search, depends on sample rate and bin size
int Flip_samples;         // Exchange I & Q samples when set (invert spectrum)
double Binsize = 1;       // bin size for FFT searching for carrier (Hz)
int Quiet = 0;
double Shift = 0;


int main(int argc,char *argv[]){
  char *filename,*locale;
  struct sample *samples = MAP_FAILED;
  off_t length;
  struct stat statbuf;
  double complex *buffer = NULL;
  double complex *spectrum = NULL;
  double complex cpstep;
  double carrier_freq,cstep;
  int i,fd,start,nsamples,exitcode;
  double complex carrier;
  
  exitcode = 0;
  fd = -1;
  if((locale = getenv("LANG")) != NULL)
    setlocale(LC_ALL,locale);
  else
    setlocale(LC_ALL,"en_US.utf8"); // The world revolves around the USA...

  while((i = getopt(argc,argv,"c:qr:f")) != EOF){
    switch(i){
    case 'c':
      Shift = atof(optarg);
      break;
    case 'q':
      Quiet = 1;
      break;
    case 'r': // Sample rate
      Samprate = atof(optarg);
      break;
    case 'f':
      Flip_samples++;
      break;
    default:
      fprintf(stderr,"%s: unknown option %c\n",argv[0],i);
      exit(1);
    }
  }

  filename = argv[optind];
  if(lstat(filename,&statbuf) == -1){
    fprintf(stderr,"lstat(%s) failed: %s\n",filename,strerror(errno));
    exitcode = 1;
    goto done;
  }
  if(!S_ISREG(statbuf.st_mode)){
    fprintf(stderr,"%s is not an ordinary file\n",filename);
    exitcode = 1;
    goto done;
  }
  length = statbuf.st_size;
  if((fd = open(filename,O_RDONLY)) == -1){
    fprintf(stderr,"open(%s,readonly) failed; %s\n",filename,strerror(errno));
    exitcode = 1;
    goto done;
  }
  if((samples = mmap(NULL,length,PROT_READ,MAP_SHARED,fd,0)) == MAP_FAILED){
    fprintf(stderr,"mmap(%s,%lld) failed: %s\n",filename,
	    (long long)length,strerror(errno));
    exitcode = 1;
    goto done;
  }
  nsamples = length / sizeof(struct sample);
  if(Flip_samples && !Quiet)
    fprintf(stderr,"I & Q samples swapped (spectrum inverted)\n");

  if((buffer = fftw_alloc_complex(Fftsize)) == NULL){
    fprintf(stderr,"fftw_alloc_complex(%d) failed\n",Fftsize);
    exitcode = 2;
    goto done;
  }
  if((spectrum = fftw_alloc_complex(Fftsize)) == NULL){
    fprintf(stderr,"fftw_alloc_complex(%d) failed\n",Fftsize);
    exitcode = 2;
    goto done;
  }
  if(!Quiet)
    fprintf(stderr,"demodulating %s: %'lld bytes, %'lld samples, %'.2lf sec @ %'.1lf Hz\n",
	    argv[optind],(long long)length,(long long)nsamples,nsamples/Samprate,Samprate);
  // For each block of input samples
  for(start=0; start < nsamples; start += Fftsize){

    if(!Flip_samples){
      for(i=0; i<Fftsize; i++)
	buffer[i] = samples[start+i].i + Q*samples[start+i].q;
    } else {
      for(i=0; i<Fftsize; i++)
	buffer[i] = samples[start+i].q + Q*samples[start+i].i;
    }
  
    carrier_freq = Shift;

    // Spin down
    cstep = 2*M_PI * carrier_freq/Samprate;  // carrier phase step per sample, radians
    cpstep = cos(cstep) - Q * sin(cstep); // Complex carrier step per sample
    carrier = 1;                          // Unity-amplitude carrier starts at 0 radians
    for(i=0;i<Fftsize;i++){
      buffer[i] *= carrier;         // Spin down to baseband
      carrier *= cpstep;                  // Increment carrier phase
    }
    for(i=0; i<Fftsize; i++){
      double s;

      s = creal(buffer[i]);
      fwrite(&s,sizeof(s),1,stdout);
      s = cimag(buffer[i]);
      fwrite(&s,sizeof(s),1,stdout);
    }
  }

 done:; // Clean up and exit
  if(samples != MAP_FAILED && munmap(samples,length) == -1){
    fprintf(stderr,"munmap(%p,%ld) failed: %s\n",samples,length,strerror(errno));
  }
  if(fd != -1)
    close(fd);

  if(buffer != NULL)
    fftw_free(buffer);

  exit(exitcode);
}
