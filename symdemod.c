// Track symbol timing, emit demodulated symbols as 8-bit excess-128 soft decision values
// Phil Karn, KA9Q, June 2014

#include <stdio.h>
#include <linux/limits.h>
#include <unistd.h>
#include <stdlib.h>
#define __USE_GNU   1
#include <math.h>
#include <sys/types.h>
#include <string.h>
#include <fcntl.h>
#include <errno.h>
#include <stdint.h>
#include <getopt.h>
#include <assert.h>
#include <locale.h>
#include "timeformat.h"

#define FRAMEBITS 1024   // 1024 bits per minor frame
#define FRAMESYMBOLS (2*FRAMEBITS)
#define SYMRATE 1024.467    // Nominal symrate for 512 bps mode


int Verbose;
double Symbolsamples;
int Samprate;              // Sample rate of incoming demodulated PM
double Symrate;            // Symbol rate
int Symbolclocks;          // Clocks per symbol: 1 for 512 bps/1024 sps and up; 8 for  64 bps/128 sps
int Quiet;

double trial_demod(double *samples,int firstsample,double symbolsamples,int symbols,double gain);

int main(int argc,char *argv[]){
  int i;
  double *samples;
  int firstsample;
  char *locale;
  int symbols = 0;
  int nsamples = 0;
  int fullwater;
  int total_samples = 0;
  double window;

  if((locale = getenv("LANG")) != NULL)
    setlocale(LC_ALL,locale);
  else
    setlocale(LC_ALL,"en_US.utf8");

  Samprate = 250000;
  Symrate = SYMRATE;
  Symbolclocks = 1;
  window = 1; // Seconds of symbols to examine
  while((i = getopt(argc,argv,"w:c:r:qC:")) != EOF){
    switch(i){
    case 'w':                  // Period (sec) over which to estimate clock params
      window = atof(optarg);
      break;
    case 'q':
      Quiet = 1;
      break;
    case 'c':
      Symrate = atof(optarg);   // Initial estimate, will be updated
      if(Symrate < 1000){       // Really 1024, but just in case the clock is really slow
	Symbolclocks = rint(1024/Symrate); // Automatically change to 1024 Hz subcarrier
	fprintf(stderr,"Symbolclocks set to %d\n",Symbolclocks);
      }
      break;
    case 'r':
      Samprate = atoi(optarg);  // Fixed, assumed exact
      break;
    case 'C':
      Symbolclocks = atoi(optarg); // Override automatic setting
      break;
    }
  }

  Symbolsamples = Samprate/Symrate;
  
  // Allocate room for 1.5 frames at nominal symbol rate
  //  fullwater = FRAMESYMBOLS * Symbolsamples + Samprate/2;   // 1/2 sec of excess - hack
  fullwater = FRAMESYMBOLS * 1.5 * Symbolsamples; // 1.5 frames
  samples = malloc(fullwater * sizeof(*samples));
  assert(samples != NULL);

  firstsample = Symbolsamples/2;
  while(1){
    double maxenergy;
    int nochange;
    double clock_incr;
    int symphase;
    double gain;
    int phase_incr;
    int nsymbols;

    if(firstsample >= Samprate){ // Don't purge less than a second
      // purge old samples & slide down 
      int slide;
      
      slide = firstsample - 2*Symbolsamples; // leave a little slop in case symbol timing jumps down
      memmove(samples,&samples[slide],sizeof(*samples)*(nsamples - slide));
      nsamples -= slide;
      firstsample -= slide;
      total_samples += slide;
    }
    // Replenish input buffer
    while(nsamples < fullwater){
      int cnt,r;
      
      cnt = fullwater - nsamples;
      r = read(0,&samples[nsamples],sizeof(*samples) * cnt);
      if(r <= 0){
	// EOF or error
	break;
      }
      nsamples += r/sizeof(*samples);
    }
    if(nsamples < Symbolsamples * FRAMESYMBOLS)
      break;      // Insufficient data; all done
    
    // Full search of symbol phase at current clock estimate
    // Look at blocks of data 'window' seconds wide
    maxenergy = 0;
    symphase = 0; // Avoid uninitialized variable warning
    nsymbols = window*Symrate; // Integer number of symbols to demod
    for(i = -Symbolsamples/2; i < Symbolsamples/2; i++){
      double energy;
      
      energy = trial_demod(samples,firstsample+i,Symbolsamples,nsymbols,0.); // demod over one window
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
      if((energy = trial_demod(samples,firstsample,Symbolsamples + clock_incr,nsymbols,0.)) > maxenergy){
	maxenergy = energy;
	Symbolsamples += clock_incr;
	Symrate = Samprate / Symbolsamples;
	nochange = 0;        // we made progress, don't stop
      } else if((energy = trial_demod(samples,firstsample,Symbolsamples - clock_incr,nsymbols,0.)) > maxenergy){
	maxenergy = energy;
	Symbolsamples -= clock_incr;
	Symrate = Samprate/Symbolsamples;
	clock_incr = -clock_incr; // Try this same direction first next time
	nochange = 0;
      } else
	nochange++;          // no improvement
      
      // See if changing phase helps
      if((energy = trial_demod(samples,firstsample + phase_incr,Symbolsamples,nsymbols,0.)) > maxenergy){
	maxenergy = energy;
	firstsample += phase_incr;
	nochange = 0;
      } else if((energy = trial_demod(samples,firstsample - phase_incr,Symbolsamples,nsymbols,0.)) > maxenergy){
	maxenergy = energy;
	firstsample += phase_incr;
	phase_incr = -phase_incr; // Try this direction first next time
	nochange = 0;
      } else
	nochange++;           // If neither a change in clock or phase helped, we'll fall through the loop
    }
    assert(firstsample >= 0);
    if(!Quiet)
      fprintf(stderr,"%s: sample %'d (%'.3lf sec, %s) symbol %'d: clock %'.4lf Hz; %'.4lf samp/sym; energy %.3lf dB\n",
	      argv[0],
	      firstsample+total_samples,
	      (double)(firstsample+total_samples)/Samprate,
	      format_hms((double)(firstsample+total_samples)/Samprate),
	      symbols,
	      Symrate,Symbolsamples,10*log10(maxenergy));
    

    // Demodulate using parameters
    gain = 100./sqrt(maxenergy); // Hack
    trial_demod(samples,firstsample,Symbolsamples,nsymbols,gain);
    firstsample += nsymbols * Symbolsamples;
    symbols += window * Symrate;
    fflush(stdout); // Keep the shell pipeline going
  }
  exit(0);
}

// Demodulate block using specified phase and clock rate, return total demodulated energy
// If gain == 0, it's a trial demod to measure energy at a clock/phase hypothesis
// If gain != 0, demodulate the symbols for real and output
double trial_demod(double *samples,int firstsample,double symbolsamples,int symbols,double gain){
  double energy,integrator,scount;
  int ind,i;
  double halfclock;
  
  energy = 0;
  ind = firstsample;       // Index of first sample to be integrated
  halfclock = (0.5 / Symbolclocks) * symbolsamples; // Width of half clock cycle
  scount = ind + halfclock; // Sample at middle of first symbol; note integer truncation

  for(i=0;i<symbols;i++){
    int j;

    integrator = 0;        // reset integrator
    // For 1024 sps and up, there's one clock per symbol. For 128 sps/64 bps, there are 8 clocks/symbol
    // For 32 sps/16 bps, there are 32 clocks/symbol
    for(j=1; j<=Symbolclocks; j++){
      // Integrate first half of clock
      for(;ind < scount; ind++)
	integrator -= samples[ind];
      
      // Integrate second half of clock
      scount += halfclock;
      for(; ind < scount; ind++)
	integrator += samples[ind];
      scount += halfclock; // update for next clock
    }

    // Integrator now has soft decision sym, 
    if(gain != 0){
      double scaled;

      // apply gain and offset-128 bias for viterbi decoder and clip
      scaled = gain * integrator + 128;
      if(scaled > 255)
	scaled = 255;
      else if(scaled < 0)
	scaled = 0;
      
      putchar((unsigned char)scaled);
    }
    integrator *= integrator; // symbol amplitude -> energy
    energy += integrator;     // Accumulate energy
  }
  //  return energy / (ind - firstsample); // Normalize for number of samples used
  return energy / symbols;
}
