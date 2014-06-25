// Demodulate and decode ISEE-3 frames
// Phil Karn, KA9Q, June 2014

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
#include <assert.h>
#include <locale.h>
#include "code.h"
#include "viterbi224.h"
#include "timeformat.h"

#define FRAMEBITS 1024   // 1024 bits per minor frame
#define FRAMESYMBOLS (2*FRAMEBITS)
#define SYNCBITS 34      // Last 34 bits of encoded tail + sync are usable
#define SYMRATE 1024.467    // Nominal symrate for 512 bps mode

#define SYNCWORD 0x12fc819fbeLL

int Verbose;
double Symbolsamples;
double Samprate = 250000;   // Sample rate of PM demodulated symbols
double Symrate = SYMRATE;   // Symbol rate

double trial_demod(short *samples,int firstsample,double symbolsamples,int symbols);

// 34 bit encoded sync symbols

int sync_vector[] = {
  0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1,
  1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};


int main(int argc,char *argv[]){
  int i,fd;
  char *filename;
  struct stat statbuf;
  off_t length;
  short *samples;
  long long nsamples;
  int firstsample;
  char *locale;
  double maxenergy;
  int symphase;
  int nochange;
  double clock_incr;
  int phase_incr;
  double window = 1; // Seconds of symbols to examine
  int ind;
  double scount;
  void *vd;
  unsigned char vdsyms[2];
  int renorms = 0;
  int symbols = 0;
  double symbol_register[FRAMESYMBOLS];
  double sync_peak_even,sync_peak_odd,sync_sum;
  int k;
  int vd_phase;
  double integrator;
  int frame = 1;
  unsigned long long decoded_data[16];
  double gain;
  size_t offset = 0;

  memset(symbol_register,0,sizeof(symbol_register));
  memset(decoded_data,0,sizeof(decoded_data));
  if((locale = getenv("LANG")) != NULL)
    setlocale(LC_ALL,locale);
  else
    setlocale(LC_ALL,"en_US.utf8");

  offset = 0;
  Samprate = 250000;
  Symrate = SYMRATE;
  while((i = getopt(argc,argv,"c:r:o:s:")) != EOF){
    switch(i){
    case 'o':
      offset = atoi(optarg);    // Skip to sample
      break;
    case 'c':
      Symrate = atof(optarg);   // Initial estimate, will be updated
      break;
    case 'r':
      Samprate = atof(optarg);  // Fixed, assumed exact
      break;
    case 's':
      Symrate = atof(optarg);   // Symbol rate
      break;
    }
  }
Symbolsamples = Samprate/Symrate;

  // Read baseband samples, downsampled to Samprate (important!)                                                                                      
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
  printf("%s: %'lld samples; %'.3lf sec (%s) @ %'.1lf Hz\n",
	 filename,nsamples,nsamples/Samprate,format_hms(nsamples/Samprate),Samprate);
  vd = create_viterbi224(250); // Must be bigger than the traceback depth
  assert(vd != NULL);
  init_viterbi224(vd,0); // See what happens

  // Look at blocks of data 'window' seconds wide
  // A 1 second window gives a SNR of 10*log10(Symrate) = ~30 dB over the symbol SNR
  for(firstsample = offset + Symbolsamples/2;
      firstsample + FRAMESYMBOLS*Symbolsamples < nsamples;
      firstsample += FRAMESYMBOLS*Symbolsamples){

    // Full search of symbol phase at current clock estimate
    maxenergy = 0;
    symphase = 0; // Avoid uninitialized variable warning
    for(i = -Symbolsamples/2; i < Symbolsamples/2; i++){
      double energy;
      
      energy = trial_demod(samples,firstsample+i,Symbolsamples,FRAMESYMBOLS); // demod over one window
      if(energy > maxenergy){
	maxenergy = energy;
	symphase = i;
      }
    }
    firstsample += symphase;
    symphase = 0;
    // Fine-tune the symbol clock and phase estimates
    clock_incr = 0.5 * Symbolsamples/(window*Samprate); // +/-0.5 sample at end of window
    phase_incr = 1;                       // +/-1 sample throughout window
    for(nochange = 0; nochange < 2;){
      double energy;
      
      // Adjust clock, see if it improves
      if((energy = trial_demod(samples,firstsample,Symbolsamples + clock_incr,FRAMESYMBOLS)) > maxenergy){
	maxenergy = energy;
	Symbolsamples += clock_incr;
	Symrate = Samprate / Symbolsamples;
	nochange = 0;        // we made progress, don't stop
      } else if((energy = trial_demod(samples,firstsample,Symbolsamples - clock_incr,FRAMESYMBOLS)) > maxenergy){
	maxenergy = energy;
	Symbolsamples -= clock_incr;
	Symrate = Samprate/Symbolsamples;
	clock_incr = -clock_incr; // Try this same direction first next time
	nochange = 0;
      } else
	nochange++;          // no improvement
      
      // See if changing phase helps
      if((energy = trial_demod(samples,firstsample + phase_incr,Symbolsamples,FRAMESYMBOLS)) > maxenergy){
	maxenergy = energy;
	firstsample += phase_incr;
	nochange = 0;
      } else if((energy = trial_demod(samples,firstsample - phase_incr,Symbolsamples,FRAMESYMBOLS)) > maxenergy){
	maxenergy = energy;
	firstsample += phase_incr;
	phase_incr = -phase_incr; // Try this direction first next time
	nochange = 0;
      } else
	nochange++;           // If neither a change in clock or phase helped, we'll fall through the loop
    }
    printf("Frame %'d starting at sample %'d (%'.3lf sec, %s): clock %'.4lf Hz; %'.4lf samp/sym; energy %.3lf dB\n",
	    frame,firstsample,firstsample/Samprate,format_hms(firstsample/Samprate),Symrate,Symbolsamples,10*log10(maxenergy));

    fflush(stdout);
    // Demodulate using parameters
    ind = firstsample;       // Index of first sample to be integrated
    scount = ind + 0.5 * Symbolsamples; // Sample at middle of first symbol; note integer truncation
    sync_peak_even = sync_peak_odd = 0;
    for(i=0;i<FRAMESYMBOLS;i++){
      // Integrate first half of symbol
      integrator = 0;        // reset integrator
      for(;ind < scount;ind++)
	integrator -= samples[ind];

      // Integrate second half of symbol
      scount += 0.5 * Symbolsamples;   
      for(; ind < scount;ind++)
	integrator += samples[ind];
      scount += 0.5 * Symbolsamples; // update for next symbol

      symbols++;
      // Accumulate in shift register
      symbol_register[i] = integrator;
      // Run sync correlator
      sync_sum = 0;
      for(k=0;k<34;k++){
	if(sync_vector[k])
	  sync_sum += symbol_register[(FRAMESYMBOLS + i - 34 + k)% FRAMESYMBOLS];
	else
	  sync_sum -= symbol_register[(FRAMESYMBOLS + i - 34 + k)% FRAMESYMBOLS];
      }	  
      if((i % 2) == 0){
	if(sync_sum > sync_peak_even)
	  sync_peak_even = sync_sum;
      } else {
	if(sync_sum > sync_peak_odd)
	  sync_peak_odd = sync_sum;
      }
    }
    // Use the best sync peak just to determine Viterbi decoder phase
    vd_phase = sync_peak_even < sync_peak_odd;
    //    printf("vd phase = %d\n",vd_phase);
    
    gain = 75./sqrt(maxenergy); // Hack
    for(i=0;i<FRAMESYMBOLS;i++){
      // Integrator now has soft decision sym, apply gain and bias for viterbi decoder and clip
      integrator = gain * symbol_register[i] + 128;
      if(integrator > 255)
	integrator = 255;
      else if(integrator < 0)
	integrator = 0;

      vdsyms[vd_phase] = integrator;
      if(vd_phase == 1){
	int bit;
	renorms += update_viterbi224_blk(vd,vdsyms,1);

	// Chain back and decode a bit
	// This is far enough back for the paths to merge, we don't
	// need to waste time finding the best path
	bit = decodebit_viterbi224(vd,200,0); // Cannot be larger than value given to create_viterbi
	// 1024-bit shift register with decoded bits
	for(k=0;;k++){
	  decoded_data[k] <<= 1;
	  if(k == 15)
	    break;
	  decoded_data[k] |= (decoded_data[k+1] >> 63);
	}
	decoded_data[15] |= bit;
	
	//	printf("symbols %6u renorms %02d data %016llx\n",i,renorms,decoded_data[15]);
	if((decoded_data[15] & 0xffffffffffL) == SYNCWORD){
	  int k;
	  
	  for(k=0;k<16;k++){
	    int n;
	    
	    for(n=56;n >= 0;n -= 8){
	      printf("%02llx",(decoded_data[k] >> n) & 0xff);
	      if(n == 0 && (k % 2) == 1)
		putchar('\n');
	      else
		putchar(' ');
	    }
	  }
	}
      }
      vd_phase = !vd_phase;
    }
    putchar('\n');
    //    printf("sync peak even %lf odd %lf\n",sync_peak_even,sync_peak_odd);
    frame++;
    fflush(stdout);
  }
  exit(0);
}

double trial_demod(short *samples,int firstsample,double symbolsamples,int symbols){
  double energy,integrator;
  int ind,i;
  double scount;

  energy = 0;
  ind = firstsample;       // Index of first sample to be integrated
  scount = ind + 0.5 * symbolsamples; // Sample at middle of first symbol; note integer truncation

  for(i=0;i<symbols;i++){
    integrator = 0;        // reset integrator
    // Integrate first half of symbol
    for(;ind < scount;ind++)
      integrator -= samples[ind];

    // Integrate second half of symbol
    scount += 0.5 * symbolsamples;   
    for(; ind < scount;ind++)
      integrator += samples[ind];

    scount += 0.5 * symbolsamples; // update for next symbol

    integrator *= integrator; // symbol amplitude -> energy
    energy += integrator;     // Accumulate energy
  }
  //  return energy / (ind - firstsample); // Normalize for number of samples used
  return energy / symbols;
}


