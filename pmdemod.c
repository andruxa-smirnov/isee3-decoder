// ICE PM demodulator
// Reads stereo signal with I&Q channels, finds residual carrier, spins it down to zero,
// outputs baseband data on stdout as 16-bit signed short ints in machine format
// 26 June 2014, Phil Karn, KA9Q
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#define _USE_GNU 1
#include <math.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>
#include <fcntl.h>
#include <errno.h>
#include <complex.h> // Must be before fftw3.h so it understands native complex double
#include <fftw3.h>
#include <locale.h>
#include <assert.h>
#include "timeformat.h"


#undef I
#define Q _Complex_I  // Less likely to cause error than the default I

// Complex pair of 16-bit, 2's complement receiver samples
struct sample {
  int16_t i;    // In-phase channel
  int16_t q;    // Quadrature channel
};

double Samprate;          // Sample rate, Hz
int Fftsize;              // Size of FFT for carrier search, depends on sample rate and bin size
int Flip_samples;         // Exchange I & Q samples when set (invert spectrum)
double Binsize;           // bin size for FFT searching for carrier (Hz)
int Quiet;
double Carrier_search_freq;     // Starting carrier frequency estimate, Hz
double Finish_carrier;    // Ending carrier frequency estimate, Hz
double Search_width;      // Carrier frequency search range around estimate, Hz
double Doppler_rate;      // Doppler rate, Hz/s

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
  char *locale;

  struct stat statbuf;
  fftw_plan ff = NULL;
  double complex *buffer = NULL;
  double complex *spectrum = NULL;
  double carrier_freq;
  int i,exitcode,lfftsize;
  double cn0_threshold,cn0 = -999;
  complex double loaccel;         // Doppler adjustment, rad/sample^2 (frequency rate)
  FILE *input = NULL;
  long long total_samples = 0;

  if((locale = getenv("LANG")) != NULL)
    setlocale(LC_ALL,locale);
  else
    setlocale(LC_ALL,"en_US.utf8"); // The world revolves around the USA...

  exitcode = 0;

  // Some explicit defaults
  Quiet = 0;
  Samprate = 250000;
  Carrier_search_freq = 0;
  Search_width = 0;
  Doppler_rate = 0;
  Binsize = 4;
  Samprate = 250000;
  cn0_threshold = 21; // Corresponds to roughly Eb/N0 = 0 dB for mod of 1.1 radians at 512 bps

  while((i = getopt(argc,argv,"S:W:D:r:fb:qt:")) != EOF){
    switch(i){
    case 'S': // Starting carrier frequency estimate, Hz
      Carrier_search_freq = atof(optarg);
      break;
    case 'W': // Frequency tolerance around estimate, Hz
      Search_width = atof(optarg);
      break;
    case 'D': // Doppler rate, Hz/s
      Doppler_rate = atof(optarg);
      break;
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
  if(fabs(Carrier_search_freq) > Samprate/2){
    fprintf(stderr,"%s: Carrier frequency of %'.1lf Hz is outside Nyquist bandwidth at %'.1lf Hz sample rate. Must be between +/- %'.1lf Hz\n",
	    argv[0],Carrier_search_freq,Samprate,Samprate/2);
    exitcode = 1;
    goto done;
  }
  if(Search_width < 0)
    Search_width = fabs(Search_width); // Be tolerant, you know what he means

  if(Search_width > Samprate/2){
    fprintf(stderr,"%s: Search width > 1/2 Nyquist rate; reduced to +/- %'.1lf Hz\n",argv[0],Samprate/2);
    Search_width = Samprate/2;
  }
  lfftsize = nearbyint(log2(Samprate / Binsize));
  Fftsize = 1 << lfftsize;
  Binsize = Samprate/Fftsize;    // Frequency width of each FFT bin, Hz

  // If carrier is not specified, we'll lock onto the first frequency that exceeds our C/No threshold
  // FFTW3 can handle arbitrary buffer sizes, but round to a power of 2
  // just to keep things fast
  if(!Quiet)
    fprintf(stderr,"%s: FFT bin size %'.4lf Hz; Start carrier %'.4lf Hz; Doppler %'.6lf Hz/s; Search range +/-%'.1lf Hz\n",
	    argv[0],Samprate/Fftsize,Carrier_search_freq,Doppler_rate,Search_width);

  // Convert Doppler values to internal units
  {
    double drate;
    drate = Doppler_rate * 2 * M_PI / (Samprate*Samprate); // Doppler in radians/sample^2
    loaccel = cos(drate) + Q*sin(drate); // Increment in LO frequency per sample
  }
  if(Flip_samples && !Quiet)
    fprintf(stderr,"%s: I & Q samples swapped (spectrum inverted)\n",argv[0]);

  if((buffer = fftw_alloc_complex(Fftsize)) == NULL){
    fprintf(stderr,"%s: fftw_alloc_complex(%d) failed\n",argv[0],Fftsize);
    exitcode = 2;
    goto done;
  }
  if((spectrum = fftw_alloc_complex(Fftsize)) == NULL){
    fprintf(stderr,"%s: fftw_alloc_complex(%d) failed\n",argv[0],Fftsize);
    exitcode = 2;
    goto done;
  }
  // Set up complex FFT. Note: non-destructive of input buffer, as we'll need it
  fftw_import_system_wisdom();
  if((ff = fftw_plan_dft_1d(Fftsize,buffer,spectrum,FFTW_FORWARD,FFTW_ESTIMATE)) == NULL){
    fprintf(stderr,"%s: Can't set up FFT\n",argv[0]);
    exitcode = 2;
    goto done;
  }

  // Set up input
  if(argc > optind){
    // Read from specified file
    if((input = fopen(argv[optind],"r")) == NULL){
      fprintf(stderr,"%s: Can't read %s; %s\n",argv[0],argv[optind],strerror(errno));
      exitcode = 1;
      goto done;
    }
  } else {
    input = stdin; // Read from standard input
  }
  // Find out what we're reading from
  if(isatty(fileno(input))){
    fprintf(stderr,"%s: Can't read from terminal\n",argv[0]);
    exitcode = 1;
    goto done;
  }
  if(fstat(fileno(input),&statbuf) == -1){
    fprintf(stderr,"%s: fstat of input failed: %s\n",argv[0],strerror(errno));
    exitcode = 1;
    goto done;
  }
  if(S_ISDIR(statbuf.st_mode)){
    fprintf(stderr,"%s: Can't read a directory\n",argv[0]);
    exitcode = 1;
    goto done;
  }
  if(S_ISREG(statbuf.st_mode)){
    if(!Quiet){
      long long nsamples;

      nsamples = statbuf.st_size / sizeof(struct sample);
      fprintf(stderr,"%s: demodulating %'lld bytes; %'lld samples; %'.2lf sec @ %'.1lf Hz\n",argv[0],
	      (long long)statbuf.st_size,
	      nsamples,nsamples/Samprate,Samprate);
    }
  }
  while(1){
    // Load up FFT input buffer. Ignore remainder and quit if we don't have enough
    if(!Flip_samples){
      int i,n;
      struct sample sample;

      for(i=0; i<Fftsize; i++){
	n = fread(&sample,sizeof(sample),1,input);
	if(n < 1){
	  exitcode = 0;
	  goto done;
	}
	buffer[i] = sample.i + Q*sample.q;
      }
    } else { // I & Q samples swapped (invert spectrum)
      int i,n;
      struct sample sample;

      for(i=0; i<Fftsize; i++){
	n = fread(&sample,sizeof(sample),1,input);
	if(n < 1){
	  exitcode = 0;
	  goto done;
	}
	buffer[i] = sample.q + Q*sample.i;
      }
    }
    // Apply Doppler chirp if specified
    if(Doppler_rate != 0){
      complex double lophase;         // Doppler adjustment oscillator phase
      complex double lofreq;          // Doppler adjustment, rad/sample (frequency)
      int i;

      lofreq = 1;   // Start LO at 0 Hz
      lophase = 1;  //         and 0 phase
      for(i=0; i<Fftsize; i++){
	buffer[i] *= conj(lophase);   // Shift *down* by LO
	lofreq *= loaccel;            // Double integration: increase LO frequency
	lophase *= lofreq;            //                     and LO phase
      }
    }

    // Search for carrier frequency with FFT
    // This admittedly doesn't account for scalloping losses, so we might not actually get the peak bin
    {
      int firstbin,lastbin,i;
      int peak,next,prev;
      double energy,maxenergy,ap,dp,am,dm,d;

      fftw_execute(ff);

      // Limit search range only when locked and tracking
      // This minimizes the risk of a noise spike pulling us off frequency from a weak carrier
      if(Search_width != 0 && cn0 > cn0_threshold){
	// Clip search range limits to edges of Nyquist band
	if(Carrier_search_freq - Search_width <= -Samprate/2)
	  firstbin = 0;
	else {
	  firstbin = (Carrier_search_freq - Search_width)/Binsize;
	  if(firstbin < 0)
	    firstbin += Fftsize;
	}	
	if(Carrier_search_freq + Search_width >= Samprate/2)
	  lastbin = Fftsize/2 - 1;
	else {
	  lastbin = (Carrier_search_freq + Search_width)/Binsize;
	  if(lastbin < 0)
	    lastbin += Fftsize;
	}
      } else {
	firstbin = 0;
	lastbin = Fftsize;
      }
      if(firstbin > lastbin){
	// Swap, will happen when range straddles 0 Hz
	int tmp;

	tmp = firstbin;
	firstbin = lastbin;
	lastbin = tmp;
      }
      assert(0 <= firstbin && firstbin <= Fftsize && 0 <= lastbin && lastbin <= Fftsize);

      // Look for bin with most energy
      peak = -1;
      maxenergy = 0;
      for(i=firstbin; i<lastbin; i++){
	energy = cenergy(spectrum[i]);
	if(energy >= maxenergy){ // so it always succeeds at least the first time and peak can't be untouched
	  // The most energetic frequency bin contains the carrier
	  maxenergy = energy;
	  peak = i;
	}
      }
      assert(peak != -1);
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
    }
      
    {
      // Spin down to baseband with arbitrary phase
      double cstep,carrier_sq,carrier_sum;
      double complex cpstep,carrier,dc;
      int i;

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
      cn0 = 10*log10(0.5*Samprate*carrier_sum/(carrier_sq-carrier_sum)); // Update C/N0 estimate in dB

      if(cn0 > cn0_threshold)
	Carrier_search_freq = carrier_freq; // Center this frequency in search window

      if(!Quiet)
	fprintf(stderr,"%s: sample %'lld (%'.3lf sec, %s); carrier %'.1lf Hz; C/No = %'.2lf dB%s\n",argv[0],
		total_samples, total_samples/Samprate, format_hms(total_samples/Samprate),carrier_freq,
		cn0,cn0 >= cn0_threshold ? " locked" : "");
      for(i=0; i<Fftsize; i++){
	short s;
	// Drop 3 db to ensure clipping can't occur
	// Worst case is an input sample of abs(I) = abs(Q) = 32767
	// that gets rotated to the Q axis
	s = (short) (cimag(buffer[i]) * M_SQRT1_2);
	fwrite(&s,sizeof(s),1,stdout);
      }
    }
    fflush(stdout);
    total_samples += Fftsize;
  }
  done:; // Clean up and exit
  if(input != NULL)
    fclose(input);

  if(ff != NULL)
    fftw_destroy_plan(ff);

  if(buffer != NULL)
    fftw_free(buffer);

  if(spectrum != NULL)
    fftw_free(spectrum);

  exit(exitcode);
}
