// Process baseband PM receiver output from ISEE-3/ICE, find frame sync, find symbol sync, and decode
// Phil Karn, KA9Q, May 2014


#include <stdio.h>
#include <linux/limits.h>
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

#define FRAMEBITS 1024   // 1024 bits per minor frame
#define SYNCBITS 34      // Last 34 bits of encoded tail + sync are usable
#define SYMRATE 1024.475    // Nominal symrate

#define SYNC_FAIL (-1234567890) // An unlikely value

int Verbose;
int Synclen;
double *Syncvector;
double Framesamples;
double Symbolsamples;
double Samprate = 250000;    // Downconverted after PM demodulation by sox
double Symrate = SYMRATE;      // Symbol rate

// Inverse error functions, for estimating Eb/No
double erf1(double z);
double erfc1(double z);

// FFT correlator stuff
fftw_plan Corr_ff1,Corr_ff2,Corr_ffr;
double *Corr_input,*Corr_vector;
double complex *Corr_vector_transform;
double complex *Corr_data_transform;
double *Corr_result;
int Corr_size;

// Set up sync vector for search correlator
// Use current estimates of sample rate
void generate_sync(double symrate){
  int ind,k;
  unsigned char data[10];
  unsigned char symbols[2*8*10];

  // The last 5 bytes (40 bits) of each 1024-bit minor frame are constants given below
  // Run them through the encoder and use the last 34 symbols as the sync vector
  // Those are the only invariant symbols from the encoder, after the user data has been flushed
  memset(data,0,sizeof(data));
  data[0] = 0x12;
  data[1] = 0xfc;
  data[2] = 0x81;
  data[3] = 0x9f;
  data[4] = 0xbe;
  encode(symbols,data,10,0);
#if 1
  for(k=0;k<80;k++)
    printf(" %d",symbols[k]);
  printf("\n");
#endif
  // Update numbers based on current sample rate
  Symrate = symrate;
  Synclen = SYNCBITS * Symbolsamples + 1; // Fudge to prevent off-by-one overrun

  printf("Symbol rate: %'.3lf Hz; samples/sym: %'.3lf; samples/frame: %'.1lf; samples in sync: %'d\n",
	  symrate,Samprate/symrate,2*FRAMEBITS*Symbolsamples,Synclen);

  if(Syncvector != NULL)
    free(Syncvector);
  if((Syncvector = malloc(sizeof(*Syncvector) * (int)Synclen)) == NULL){
    fprintf(stderr,"Can't malloc sync vector\n");
    exit(3);
  }
  // Manchester encode last 34 symbols in sequence
  ind = 0;
  for(k=0;k<SYNCBITS;k++){
    // First half of manchester symbol
    for(;ind < (k+0.5) * Symbolsamples;ind++)
      Syncvector[ind] = symbols[k+80-SYNCBITS] ? -1 : 1;

    // Second half
    for(;ind < (k+1) * Symbolsamples;ind++)
      Syncvector[ind] = symbols[k+80-SYNCBITS] ? 1 : -1;
  }
  assert(ind <= Synclen);
#if 0
  for(k=0;k<Synclen;k++){
    putchar(Syncvector[k] == 1 ? '+':'-');
  }
  putchar('\n');
  printf("sync vector done\n");
#endif

  // Set up FFT correlator
  Corr_size = 2*Framesamples;
  Corr_size = 1 << 20; // hack!! change to round Corr_size up to next power of 2
  Corr_input = fftw_alloc_real(Corr_size);
  assert(Corr_input != NULL);
  Corr_vector = fftw_alloc_real(Corr_size);
  assert(Corr_vector != NULL);
  Corr_result = fftw_alloc_real(Corr_size);
  assert(Corr_result != NULL);
  Corr_vector_transform = fftw_alloc_complex(Corr_size);
  assert(Corr_vector_transform != NULL);
  Corr_data_transform = fftw_alloc_complex(Corr_size);
  assert(Corr_data_transform != NULL);
  
  fftw_import_system_wisdom();
  Corr_ffr = fftw_plan_dft_c2r_1d(Corr_size,Corr_data_transform,Corr_result,FFTW_ESTIMATE);
  assert(Corr_ffr != NULL);
  Corr_ff2 = fftw_plan_dft_r2c_1d(Corr_size,Corr_input,Corr_data_transform,FFTW_ESTIMATE);
  assert(Corr_ff2 != NULL);
  Corr_ff1 = fftw_plan_dft_r2c_1d(Corr_size,Corr_vector,Corr_vector_transform,FFTW_ESTIMATE);
  
  // Load up the sync vector
  for(k=0;k<Synclen;k++)
    Corr_vector[k] = Syncvector[k];
  
  // Zero pad
  for(;k<Corr_size;k++)
    Corr_vector[k] = 0;
  
  fftw_execute(Corr_ff1); // Compute transform of sync vector
  // Take complex conjugate so we don't have to do it every time
  for(k=0;k<Corr_size;k++)
    Corr_vector_transform[k] = conj(Corr_vector_transform[k]);
}


// Sync vector correlation with FFT - Thanks Mario, DL5MLO
int fft_sync_search(short *samples,int low,int high,int samp){
  int i;
  int peakindex;
  double maxpeak;
  int nonzero;

  // Load up the data, pad with zeroes
  nonzero = 1;
  for(i=0;i<Framesamples;i++){
    Corr_input[i] = samples[i];
    if(samples[i] != 0)
      nonzero = 0;
  }
  if(nonzero)
    return SYNC_FAIL; // Buffer is all zeros

  for(;i < Corr_size;i++)
    Corr_input[i] = 0;

  fftw_execute(Corr_ff2);

  // Multiply transforms (complex conj has already been done in sync_setup)
  for(i=0;i<Corr_size;i++)
    Corr_data_transform[i] *= Corr_vector_transform[i];
    
  // Inverse transform
  fftw_execute(Corr_ffr);

  if(samp != -1){
    // Dump to plot file
    FILE *plot;
    char plotname[PATH_MAX];

    // Don't overwrite  if it already exists
    snprintf(plotname,PATH_MAX,"sync.%d.plot",samp);
    plot = fopen(plotname,"w");
    fprintf(plot,"signed double\n");
    for(i=0;i<Corr_size;i++){
      fprintf(plot,"dot %d %lf\n",i,Corr_result[i]);
    }
    fclose(plot);
  }
  // Now look at result for biggest peak
  peakindex = -1;
  maxpeak = 0;

  assert(low >= 0 && high >= 0);
  if(high > Corr_size)
    high = Corr_size;
  for(i=low;i<high;i++){
    if(Corr_result[i] > maxpeak){
      maxpeak = Corr_result[i];
      peakindex = i;
    }
  }
  if(maxpeak == 0){
    // All 0's, squelch closed, fail
    return SYNC_FAIL;
  }
  if(peakindex > Corr_size/2)
    peakindex = Corr_size - peakindex;

  return peakindex;
}


int main(int argc,char *argv[]){
  unsigned char data[FRAMEBITS/8]; // One minor frame
  unsigned char symbols[2*FRAMEBITS];
  unsigned char nsymbols[2*FRAMEBITS];
  int i,fd;
  char *filename;
  struct stat statbuf;
  off_t length;
  short *samples;
  int  nsamples;
  int start;
  int startsync = SYNC_FAIL,endsync = SYNC_FAIL;
  void *vd; // Viterbi decoder handle
  int symerrors;
  int firstsample;
  int frame;
  int begin;
  int maxmetric,minmetric;
  int ind;
  char *locale;
  double clock_tolerance = 5; // Maximum allowable clock offset, samples/frame
  
  if((locale = getenv("LANG")) != NULL)
    setlocale(LC_ALL,locale);
  else
    setlocale(LC_ALL,"en_US.utf8");

  begin = 0;
  while((i = getopt(argc,argv,"c:o:r:t:")) != EOF){
    switch(i){
    case 't':
      clock_tolerance = atof(optarg);
      break;
    case 'c':
      Symrate = atof(optarg);
      break;
    case 'r':
      Samprate = atof(optarg);
      break;
    case 'o':
      begin = atoi(optarg); // Starting sample, default 0
      break;
    }
  }
  Symbolsamples = Samprate/Symrate;
  Framesamples = Symbolsamples * 2*FRAMEBITS;

  // Read baseband samples, downsampled to Samprate (important!) and look for sync vector 
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
  if((samples = mmap(NULL,length,PROT_READ,MAP_SHARED,fd,0)) == MAP_FAILED){
    fprintf(stderr,"mmap(%s,%lld) failed: %s\n",filename,
	    (long long)length,strerror(errno));
    close(fd);
    exit(1);
  }
  nsamples = length / sizeof(*samples);
  printf("%s: %'d samples, %'.3lf seconds @ %.1lf Hz\n",filename,
	  nsamples,nsamples/Samprate,Samprate);
  generate_sync(Symrate);
  vd = create_viterbi224(FRAMEBITS);
  if(vd == NULL){
    fprintf(stderr,"Can't create viterbi decoder\n");
    exit(2);
  }
  
  for(frame=1;begin + Framesamples < nsamples;frame++){
    int low,high;

    // Look for starting frame sync
    if(startsync == SYNC_FAIL){
      // No endsync from previous frame, find initial sync
      while(1){
	startsync = fft_sync_search(&samples[begin],0,(int)Framesamples,begin);
	if(startsync != SYNC_FAIL){
	  startsync += begin;
	  break;
	}
	begin += Framesamples; // Skip a frame time and try again
	printf("Start sync search failure, skip to %'d\n",begin);
      }
    }
    assert(startsync != SYNC_FAIL);
    // Look for ending frame sync
    // Start halfway through the frame, but actually search the middle
    // of the result, which straddles the next sync. Only search in the
    // allowable range to keep the clock from being driven into the weeds by noise
    while(1){
      start = startsync + Framesamples/2;
      low = 0.5 * Framesamples - clock_tolerance; // very tight tolerance, +/- 5 samples
      high = 0.5 * Framesamples + clock_tolerance;
      endsync = fft_sync_search(&samples[start],low,high,-1);
      if(endsync != SYNC_FAIL){
	endsync += start;
	break;
      }
      // Go back to look for another start sync
      begin = startsync + Framesamples;
      printf("End sync search failure, skip to %'d\n",begin);
      startsync = SYNC_FAIL;
      goto again;
    }
    // Got both start and end sync, proceed
    i = startsync/Samprate;
    printf("Frame %'d @ sample %'d (%'d:%02d)\n",frame,startsync,i/60,i % 60);

    // Update estimate of sample clock rate et al
#if 0 // Hack - not let it change
    Framesamples = endsync - startsync; // Actual samples in current frame (int)
    Symbolsamples = Framesamples/(2*FRAMEBITS); // Actual samples in symbol (float)
    Symrate = Samprate / Symbolsamples;
#endif
    printf("Symbol rate: %'.3lf Hz; samples/sym: %'.3lf; samples/frame: %'.1lf\n",
	  Symrate,Symbolsamples,Framesamples);
#if 0    
    // If clock estimate has changed a lot, regenerate the sync vector with the new clock
    if(abs(SYNCBITS*Samprate/Symrate - Synclen) > Symbolsamples/10.){
      printf("Clock changed; regenerate sync\n");
      generate_sync(Symrate);
    }
#endif
    // startsync points to the first symbol in the 34-bit sync sequence,
    // so we want to skip past it to the first symbol in the new frame
    firstsample = SYNCBITS*Symbolsamples + startsync;
    for(i=0; i < 2*FRAMEBITS; i++){
      double sum;
      int midpoint,last;
      
      ind = firstsample + i * Symbolsamples;
      midpoint = firstsample + (i+0.5)*Symbolsamples;
      last = firstsample + (i+1.0)*Symbolsamples;
      // Integrate bit, remove manchester
      sum = 0;
      for(;ind < midpoint;ind++)   // first half of symbol
	sum -= samples[ind];
      for(; ind < last;ind++)
	sum += samples[ind];
      
      sum += 128; // Offset-128 for Viterbi decoder
      symbols[i] = (sum > 255) ? 255 : ((sum < 0) ? 0 : sum); // Clip to range 0-255
    }
    // We start with the encoder having just transmitted these five fixed bytes,
    // 3 bytes of coder dump and 2 bytes of sync:
    // 12 fc 81 9f be
    init_viterbi224(vd,0x819fbe);
    update_viterbi224_blk(vd,symbols,FRAMEBITS);
    chainback_viterbi224(vd,data,FRAMEBITS,0x819fbe);
    maxmetric = max_metric_viterbi224(vd);
    minmetric = min_metric_viterbi224(vd);

    for(i=0; i<128; i++){
      printf("%02x",data[i]);
      if((i % 16) == 15)
	putchar('\n');
      else
	putchar(' ');
    }
    // Re-encode and compare to received symbols to count channel errors
    encode(nsymbols,data,FRAMEBITS/8,0x819fbe);
    symerrors = 0;
    for(i=0;i<2*FRAMEBITS;i++){
      if(nsymbols[i] != (symbols[i] > 128)){
	symerrors++;
#if 0
	printf("sym error %d: %d != %d\n",i,nsymbols[i],symbols[i]);
#endif
      }
    }
    printf("Viterbi path metric range %'d - %'d, diff %'d\n",minmetric,maxmetric,maxmetric-minmetric);
    if(symerrors){
      double esn0;

      esn0 = erfc1(2.*symerrors/(2*FRAMEBITS)); // Amplitude ratio
      esn0 *= esn0; // square to get power ratio
      // ebn0 = 2 * esn0 for rate 1/2 code
      printf("re-encode symbol errors: %'d/%'d; estimated Eb/No = %.2lf dB\n",
	       symerrors,2*FRAMEBITS,10*log10(2*esn0));
    } else {
      printf("No re-encode symbol errors; estimated Eb/No > %.2lf dB\n",10.5); // hack; 7.5 dB has a BER < 1/2048
    }
    putchar('\n');
    fflush(stdout);
    startsync = endsync;
  again:;
  }
  delete_viterbi224(vd);

  exit(0);
}

// Inverse erfc function, used to estimate SNR from BER
double erf1(double z){
  int k,m;
  long double sum;
  static long double c[100];

  if(c[0] == 0){
    // Initialize C's
    c[0] = 1;
    for(k=1;k<100;k++){
      c[k] = 0;
      for(m=0;m<k;m++){
	assert(c[m] != 0 && c[k-1-m] != 0);
	c[k] += c[m] * c[k-1-m] / ((m+1)*(2*m+1));
      }
      //      printf("c[%d] = %Lg\n",k,c[k]);
    }
  }

  sum = 0;
  for(k=0;k<100;k++)
    sum += (c[k] / (2*k+1)) * powl(z/M_2_SQRTPIl,2*k+1);

  return sum;
}

// Inverse complementary error function
double erfc1(double z){
  return erf1(1-z);
}
