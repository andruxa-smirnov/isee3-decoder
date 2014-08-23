// Autocorrelate baseband signal
// Phil Karn, KA9Q, June 2014

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#define __USE_GNU   1
#include <math.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>
#include <fcntl.h>
#include <errno.h>
#include <complex.h> // Must be before fftw3.h
#include <fftw3.h>
#include <assert.h>
#include <locale.h>
#include "code.h"
#include "viterbi224.h"

#undef I
#define Q _Complex_I  // Less likely to cause error than the default I

#define SPECTRUM "spectrum.plot"
#define AUTOSPECT "autospect.plot"
#define AUTOCORRELATION "autocorr.plot"

int Verbose;
double Samprate = 250000;    // Downconverted after PM demodulation by sox

int Pagesize;


int main(int argc,char *argv[]){
  // Set up FFT correlator
  char *filename;
  int fd;
  struct stat statbuf;
  size_t length;
  short *samples;
  long long i,k;
  char *locale;
  size_t nsamples;
  off_t offset;
  FILE *plot;
  long long corr_size = 0;
  // FFT correlator stuff
  fftw_plan corr_ffd,corr_ffr;
  double *corr_time;
  double complex *corr_freq;

  Pagesize = sysconf(_SC_PAGE_SIZE);

  offset = 0;
  while((i = getopt(argc,argv,"o:r:s:")) != EOF){
    switch(i){
    case 's':
      corr_size = atoll(optarg);
      break;
    case 'r':
      Samprate = atof(optarg);
      break;
    case 'o':
      offset = ~(Pagesize-1) & (atoll(optarg) * sizeof(*samples)); // Starting sample, default 0
      break;
    }
  }
  if((locale = getenv("LANG")) != NULL)
    setlocale(LC_ALL,locale);
  else
    setlocale(LC_ALL,"en_US.utf8");

  // Read baseband samples
  filename = argv[optind];
  if(lstat(filename,&statbuf) == -1){
    fprintf(stderr,"lstat(%s) failed: %s\n",filename,strerror(errno));
    exit(1);
  }
  if(!S_ISREG(statbuf.st_mode)){
    fprintf(stderr,"%s is not an ordinary file\n",filename);
    exit(1);
  }
  length = statbuf.st_size;
  if((fd = open(filename,O_RDONLY)) == -1){
    fprintf(stderr,"open(%s,readonly) failed; %s\n",filename,strerror(errno));
    exit(1);
  }
  if((samples = mmap(NULL,length,PROT_READ,MAP_SHARED,fd,offset)) == MAP_FAILED){
    fprintf(stderr,"mmap(%s,%lld) failed: %s\n",filename,
	    (long long)length,strerror(errno));
    close(fd);
    exit(1);
  }
  nsamples = length / sizeof(*samples);
  fprintf(stderr,"%s: %'lld samples, %'.3lf seconds @ %.1lf Hz\n",filename,
	  (long long)nsamples,nsamples/Samprate,Samprate);

  corr_size = ceil(log2(nsamples)); // Round up to next power of 2
  corr_size = 1LL << corr_size;
  fprintf(stderr,"Correlator size = %'lld\n",corr_size);

  fftw_import_system_wisdom();

  // Allocate buffer for time-domain data, real only
  corr_time = fftw_alloc_real(corr_size);
  assert(corr_time != NULL);

  // Allocate buffer for complex frequency domain data
  corr_freq = fftw_alloc_complex(corr_size);
  assert(corr_freq != NULL);
  
  corr_ffd = fftw_plan_dft_r2c_1d(corr_size,corr_time,corr_freq,FFTW_ESTIMATE);
  assert(corr_ffd != NULL);

  corr_ffr = fftw_plan_dft_c2r_1d(corr_size,corr_freq,corr_time,FFTW_ESTIMATE);
  assert(corr_ffr != NULL);

  for(i=0;i< nsamples;i++)
    corr_time[i] = samples[i];

  munmap(samples,length);
  close(fd);
  fprintf(stderr,"File read\n");

  for(; i< corr_size;i++)
    corr_time[i] = 0; // Pad with zeroes

  fftw_execute(corr_ffd); // Compute transform of input


  fprintf(stderr,"FFT executed\n");

  // Display spectrum
  plot = fopen(SPECTRUM,"w");
  fprintf(plot,"double double\ntitle\nSpectrum\nxlabel\nHz\n");
  for(i=0;i<corr_size/2;i++)
    fprintf(plot,"dot %lf %lf\n",(double)i*Samprate/corr_size,cabs(corr_freq[i]));
  fclose(plot);
  fprintf(stderr,"spectrum plot in %s\n",SPECTRUM);

  // Multiply transform by its complex conjugate
  for(k=0;k<corr_size;k++)
    corr_freq[k] *= conj(corr_freq[k]);

  plot = fopen(AUTOSPECT,"w");
  fprintf(plot,"double double\ntitle\nAutocorr spectrum\nxlabel\nHz\n");
  for(i=0;i<corr_size/2;i++)
    fprintf(plot,"dot %lf %lf\n",(double)i*Samprate/corr_size,cabs(corr_freq[i]));
  fclose(plot);
  fprintf(stderr,"autocorelation spectrum plot in %s\n",AUTOSPECT);

  // Inverse transform
  fftw_execute(corr_ffr);
  fprintf(stderr,"Autocorrelation executed\n");

  plot = fopen(AUTOCORRELATION,"w");
  fprintf(plot,"double double\ntitle\nAutocorrelation\nxlabel\nsec\n");
  // freqbinsize = Samprate / corr_size
  // Skip 0 sample as it will be huge
  for(i=1;i<corr_size/2;i++){
    fprintf(plot,"dot %lf %lf\n",(double)i/Samprate,corr_time[i]);
    //    printf("dot %lf %lf\n",Samprate/(double)(i,corr_time[i]);
  }
  fclose(plot);
  fprintf(stderr,"Autocorrelation plot in %s\n",AUTOCORRELATION);

  fftw_free(corr_freq);
  fftw_free(corr_time);  // Free input as soon as we don't need it

  exit(0);
}
