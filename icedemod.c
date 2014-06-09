// ICE PM demodulator
// Reads stereo signal with I&Q channels, finds residual carrier, spins it down to zero,
// outputs baseband data on stdout as 64-bit doubles in machine format
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
int Fftsize;              // Size of FFT for carrier search, depends on sample rate and bin size
int Flip_samples;         // Exchange I & Q samples when set (invert spectrum)
double Binsize = 1;       // bin size for FFT searching for carrier (Hz)
int Quiet = 0;

// aux function for Quinn's second estimator
double tau(double x){
  return 0.25*log(3*x*x + 6*x + 1)
    - sqrt(6.)/24 * log((x + 1 - sqrt(2/3.)) / (x + 1 + sqrt(2/3.)));
}

// Return square of amplitude of complex, i.e., convert complex amplitude to energy
// Same as cabs(x)*cabs(x) avoiding the sqrt()
double cenergy(double complex x){
  return creal(x)*creal(x) + cimag(x)*cimag(x);
}

int main(int argc,char *argv[]){
  char *filename,*locale;
  struct sample *samples = MAP_FAILED;
  off_t length;
  struct stat statbuf;
  fftw_plan ff = NULL;
  double complex *buffer = NULL;
  double complex *spectrum = NULL;
  double complex dc,cpstep;
  double energy,maxenergy,carrier_freq,ap,dp,am,dm,d,cstep;
  int peak,next,prev,i,fd,start,nsamples,exitcode,lfftsize;
  double carrier_sq;
  double carrier_sum,cn0,cn0_threshold;
  double complex carrier;
  
  exitcode = 0;
  fd = -1;
  cn0_threshold = 21; // Corresponds to roughly Eb/N0 = 0 dB for mod of 1.1 radians
  if((locale = getenv("LANG")) != NULL)
    setlocale(LC_ALL,locale);
  else
    setlocale(LC_ALL,"en_US.utf8"); // The world revolves around the USA...

  while((i = getopt(argc,argv,"F:r:fb:qt:")) != EOF){
    switch(i){
    case 't':
      cn0_threshold = atof(optarg);
      break;
    case 'q':
      Quiet = 1;
      break;
    case 'b':
      Binsize = atof(optarg);
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
  // FFTW3 can handle arbitrary buffer sizes, but round to a power of 2
  // just to keep things fast
  lfftsize = nearbyint(log2(Samprate / Binsize));
  Fftsize = 1 << lfftsize;
  if(!Quiet)
    fprintf(stderr,"Sample rate %'.1lf Hz, requested carrier search bin %'.4lf Hz, FFT size rounded to 2^%'d, bin size = %'.4lf Hz\n",
	    Samprate,Binsize,lfftsize,Samprate/Fftsize);
  Binsize = Samprate/Fftsize;    // Frequency width of each FFT bin, Hz

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
  // Set up complex FFT
  fftw_import_system_wisdom();
  if((ff = fftw_plan_dft_1d(Fftsize,buffer,spectrum,FFTW_FORWARD,FFTW_ESTIMATE)) == NULL){
    fprintf(stderr,"Can't set up FFT\n");
    exitcode = 2;
    goto done;
  }
  if(!Quiet)
    fprintf(stderr,"demodulating %s: %'lld bytes, %'lld samples, %'.2lf sec @ %'.1lf Hz\n",
	    argv[optind],(long long)length,(long long)nsamples,nsamples/Samprate,Samprate);
  // For each block of input samples
  for(start=0; start < nsamples; start += Fftsize){
  
    // FFT for carrier recovery
    // Samples per FFT depends on specified carrier loop BW
    if(!Flip_samples){
      for(i=0; i<Fftsize; i++)
	buffer[i] = samples[start+i].i + Q*samples[start+i].q;
    } else {
      for(i=0; i<Fftsize; i++)
	buffer[i] = samples[start+i].q + Q*samples[start+i].i;
    }
    fftw_execute(ff);
    // The most energetic frequency bin contains the carrier
    maxenergy = 0;
    peak = -1;
    for(i=0;i<Fftsize;i++){
      energy = cenergy(spectrum[i]);
      if(energy > maxenergy){
	// Save new winner
	maxenergy = energy;
	peak = i;
      }
    }
    // Interpolate with Quinn's second estimator
    // http://www.dspguru.com/dsp/howtos/how-to-interpolate-fft-peak
    next = (peak + 1) % Fftsize;
    prev = (Fftsize + peak - 1) % Fftsize;
#if 0
    fprintf(stderr,"prev %'lg @ fft bin %'d; I = %'lg, Q = %'lg angle %'lg rad\n",
	    cenergy(spectrum[prev]),prev,creal(spectrum[prev]),cimag(spectrum[prev]),carg(spectrum[prev]));
    fprintf(stderr,"peak %'lg @ fft bin %'d; I = %'lg, Q = %'lg angle %'lg rad\n",
	    maxenergy,peak,creal(spectrum[peak]),cimag(spectrum[peak]),carg(spectrum[peak]));
    fprintf(stderr,"next %'lg @ fft bin %'d; I = %'lg, Q = %'lg angle %'lg rad\n",
	    cenergy(spectrum[next]),next,creal(spectrum[next]),cimag(spectrum[next]),carg(spectrum[next]));
#endif    
    ap = (creal(spectrum[next]) * creal(spectrum[peak]) + cimag(spectrum[next]) * cimag(spectrum[peak])) / maxenergy;
    dp = -ap/(1-ap);
    am = (creal(spectrum[prev]) * creal(spectrum[peak]) + cimag(spectrum[prev]) * cimag(spectrum[peak])) / maxenergy;
    dm = am/(1-am);
    d = (dp + dm)/2 + tau(dp*dp) - tau(dm*dm);
    carrier_freq = Binsize * (peak + d);
    if(carrier_freq > Samprate/2)
      carrier_freq -= Samprate;            // Frequency is negative

    // Spin down to baseband with arbitrary phase
    cstep = 2*M_PI * carrier_freq/Samprate;  // carrier phase step per sample, radians
    cpstep = cos(cstep) - Q * sin(cstep); // Complex carrier step per sample
    carrier = 1;                          // Unity-amplitude carrier starts at 0 radians
    dc = 0;                               // Accumulate arbitrary phase carrier at DC as we go
    for(i=0;i<Fftsize;i++){
      dc += buffer[i] *= carrier;         // Spin down to baseband, accumulate carrier at DC
      carrier *= cpstep;                  // Increment carrier phase
    }
    dc /= Fftsize;                        // Find average carrier amplitude

    // Normalize carrier amplitude and take conjugate so we'll rotate it back onto the I axis
    dc = conj(dc) / cabs(dc);
    // Rotate phase to put carrier on I axis and data on Q axis, estimate C/N0
    carrier_sq = 0;
    carrier_sum = 0;
    for(i=0; i<Fftsize; i++){
      double complex cs;

      cs = buffer[i] *= dc;
      carrier_sq += creal(cs)*creal(cs); // Sum of carrier squares
      carrier_sum += creal(cs);          // Sum of carrier
    }
    carrier_sq /= Fftsize;  // mean of carrier squares
    carrier_sum /= Fftsize; // mean carrier
    carrier_sum *= carrier_sum; // square of mean carrier
    cn0 = 10*log10(0.5*Samprate*carrier_sum/(carrier_sq-carrier_sum)); // C/N0 estimate in dB
    if(!Quiet)
      fprintf(stderr,"sample %'d (%'.3lf sec) carrier %'.1lf Hz C/No = %'.2lf dB\n",
	      start,start/Samprate,carrier_freq,cn0);
    if(cn0 > cn0_threshold){
      for(i=0; i<Fftsize; i++){
	double s;
	s = cimag(buffer[i]);
	fwrite(&s,sizeof(s),1,stdout);
      }
    } else {
      // Receiver squelch closed, emit equal number of 0's
      for(i=0;i<Fftsize;i++){
	double s;

	s = 0;
	fwrite(&s,sizeof(s),1,stdout);
      }

    }
  }
 done:; // Clean up and exit
  if(samples != MAP_FAILED && munmap(samples,length) == -1){
    fprintf(stderr,"munmap(%p,%ld) failed: %s\n",samples,length,strerror(errno));
  }
  if(fd != -1)
    close(fd);

  if(ff != NULL)
    fftw_destroy_plan(ff);

  if(buffer != NULL)
    fftw_free(buffer);

  exit(exitcode);
}
