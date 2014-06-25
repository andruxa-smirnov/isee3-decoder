// Track symbol timing, emit demodulated symbols as 8-bit excess-128 soft decision values
// Phil Karn, KA9Q, June 2014

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <fenv.h>
#define __USE_GNU   1
#include <math.h>
#include <string.h>
#include <assert.h>
#include <locale.h>
#include "timeformat.h"

#define FRAMEBITS 1024   // 1024 bits per minor frame
#define FRAMESYMBOLS (2*FRAMEBITS)
#define NOMINALCLOCK 1024.0
#define ACTUALCLOCK 1024.545058    // Measured clock rate @ 128 sps; Doppler is only about 0.01 Hz of this

int Verbose;
double Symbolsamples;      // Samples per symbol (updated with Symrate frequency)
int Samprate;              // Sample rate of incoming demodulated PM (assumed constant and exact)
double Symrate;            // Symbol rate, Hz
int Symbolclocks;          // Clocks per symbol: 1 for 512 bps/1024 sps and up; 8 for  64 bps/128 sps
int Clocktrack;            // Track clock frequency if 1, use fixed value if 0
int Quiet;                 // Suppress debug messages if 1

double trial_demod(short *samples,int firstsample,double symbolsamples,int symbols,double gain);
double timesearch(int *offset,short *samples,int firstsample,double symbolsamples,int nsymbols);

int main(int argc,char *argv[]){
  char *locale;
  short *samples;    // Input sample buffer
  int firstsample;   // Index of first unprocessed sample in buffer for this window
  int nsamples = 0;  // Actual number of samples in input buffer (irrespective of firstsample)
  int fullwater;     // Number of samples to keep in input buffer
  double window;     // Amount to process per iteration, seconds
  int nsymbols;      // number of symbols per iteration
  long long total_samples = 0; // Total samples processed
  long long total_symbols = 0;   // Count of total symbols processed
  int i;

  if((locale = getenv("LANG")) != NULL)
    setlocale(LC_ALL,locale);
  else
    setlocale(LC_ALL,"en_US.utf8");

  fesetround(FE_TONEAREST); // Set rounding mode for nearbyint()

  // Set defaults, read options
  Samprate = 250000;
  Symrate = ACTUALCLOCK;
  Symbolclocks = 1;
  window = 1.0; // Seconds of symbols to examine
  Clocktrack = 0;  // Don't track clock by default
  while((i = getopt(argc,argv,"w:c:r:qtC:")) != EOF){
    switch(i){
    case 't':
      Clocktrack = 1;
      break;
    case 'w':                  // Period (sec) over which to estimate clock params
      window = atof(optarg);
      break;
    case 'q':
      Quiet = 1;
      break;
    case 'c':
      if(!strchr(optarg,'.')){
	// If no decimal given, scale to measured clock rate
	Symrate = atof(optarg) * ACTUALCLOCK / NOMINALCLOCK;
      } else {
	Symrate = atof(optarg);   // Use exact user-specified value
      }
      if(Symrate < 1000){       // Really 1024, but just in case the clock is really slow
	Symbolclocks = rint(NOMINALCLOCK/Symrate); // Automatically change to 1024 Hz subcarrier
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
  if(!Quiet)
    fprintf(stderr,"%s: sample rate %'d Hz; estimation window %.3lf sec; clocks/symbol %d; symbol rate %.3lf Hz; tracking %s\n",
	    argv[0],Samprate,window,Symbolclocks,Symrate,Clocktrack ? "on" : "off");
  Symbolsamples = Samprate/Symrate;
  fullwater = window * 2.0 * Samprate; // 2 full windows' worth of samples
  samples = malloc(fullwater * sizeof(*samples));
  assert(samples != NULL);
  nsymbols = window*Symrate; // Integer number of symbols to demod
  firstsample = Symbolsamples/2;

  while(1){
    double maxenergy;    // Maximum energy found so far in a trial demod
    double gain;         // Gain passed to actual demodulation step to scale for Viterbi decoder
    int symphase;        // Symbol timing adjustment, samples

    if(firstsample >= window * Samprate){ // Don't purge less than a window
      // purge old samples & slide down remaining samples
      int slide;         // How many samples to remove
      
      slide = firstsample - 2*Symbolsamples; // leave a little slop in case symbol timing jumps down
      if(slide > nsamples)
	slide = nsamples;
      memmove(samples,&samples[slide],sizeof(*samples)*(nsamples - slide));
      nsamples -= slide;
      firstsample -= slide;
      total_samples += slide;
    }
    // Replenish input buffer to full water mark
    while(nsamples < fullwater){
      int cnt;

      cnt = sizeof(*samples) * (fullwater - nsamples); // bytes to read
      cnt = read(0,&samples[nsamples],cnt); // Actual bytes  read
      if(cnt <= 0)
	break;	// EOF or error

      nsamples += cnt/sizeof(*samples);
    }
    if(nsamples < window * Samprate)
      break;      // Insufficient data; all done
    
    assert(firstsample < nsamples);
    // Full search of symbol phase at current clock estimate
    // Look at blocks of data 'window' seconds wide
    maxenergy = timesearch(&symphase,samples,firstsample,Symbolsamples,nsymbols);
    firstsample += symphase;

    if(Clocktrack){
      // Fine-tune the symbol clock and phase estimates
      double clock_incr;   // How much to change the clock in one step
      int nochange;        // Count of trials that have resulted in no improvement
      int phase_incr;      // How many samples to change the symbol timing in one step

      clock_incr = 0.5 * Symbolsamples/(window*Samprate); // +/-0.5 sample at end of window
      phase_incr = 1;                       // +/-1 sample throughout window
      for(nochange = 0; nochange < 2;){
	double energy;
	
	// Adjust clock, see if it improves total demodulated energy
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
	
	// See if changing phase helps improve total demodulated energy
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
      // Update in case Symrate has changed a lot, but defer until now to avoid upsetting energy calculation
      nsymbols = window*Symrate;
    }
    assert(firstsample >= 0);
    if(!Quiet)
      fprintf(stderr,"%s: sample %'lld (%'.3lf sec, %s) symbol %'lld: clock %'.4lf Hz; %'.4lf samp/sym; timing adj %+d samples; energy %.3lf dB\n",
	      argv[0],
	      firstsample+total_samples,
	      (double)(firstsample+total_samples)/Samprate,
	      format_hms((double)(firstsample+total_samples)/Samprate),
	      total_symbols,
	      Symrate,
	      Symbolsamples,
	      symphase,
	      10*log10(maxenergy));
    

    // Demodulate using parameters
    gain = 100./sqrt(maxenergy); // Hack
    trial_demod(samples,firstsample,Symbolsamples,nsymbols,gain);
    firstsample += nsymbols * Symbolsamples;
    total_symbols += nsymbols;
    fflush(stdout); // Keep the shell pipeline going
  }
  exit(0);
}

// Demodulate block using specified phase and clock rate, return total demodulated energy
// If gain == 0, it's a trial demod to measure energy at a clock/phase hypothesis
// If gain != 0, demodulate the symbols for real and output
double trial_demod(short *samples,int firstsample,double symbolsamples,int nsymbols,double gain){
  double energy,scount;
  int ind,i,scount_int;
  double halfclock;
  
  energy = 0;
  ind = firstsample;       // Index of first sample to be integrated
  // Because there are a fractional number of samples per symbol, we keep track
  // of things with floating point numbers. But the actual integration is done
  // over an integral number of samples, rounded to the nearest one, that can fluctuate
  // between symbols as fractional samples accumulate to whole samples.
  // This can cause frequency aliasing problems if the signal really has energy
  // all the way up to the Nyquist bandwidth, but this seems unlikely at high sample rates.
  halfclock = (0.5 / Symbolclocks) * symbolsamples; // Width of half clock cycle
  scount = ind + halfclock; // First sample in next part of symbol
  scount_int = nearbyint(scount);

  for(i=0;i<nsymbols;i++){
    int j;
    long integrator;

    integrator = 0;        // reset integrator
    // For 1024 sps and up, there's one clock per symbol. For 128 sps/64 bps, there are 8 clocks/symbol
    // For 32 sps/16 bps, there are 32 clocks/symbol
    for(j=0; j<Symbolclocks; j++){
      // Integrate first half of clock cycle
      for(;ind < scount_int; ind++)
	integrator -= samples[ind];
      
      // Integrate second half of clock cycle
      scount += halfclock;
      scount_int = nearbyint(scount);
      for(; ind < scount_int; ind++)
	integrator += samples[ind];
      scount += halfclock; // update for next clock
      scount_int = nearbyint(scount);
    }
    // Integrator now has soft decision symbol
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
    // symbol amplitude -> energy
    energy += (long long)integrator * integrator;
  }
  return energy / nsymbols; //  average energy per symbol
}


// With a given clock, search all timing offsets for maximum energy in a somewhat efficient manner
double timesearch(int *symphase,short *samples,int firstsample,double symbolsamples,int nsymbols){
  int i,j,ind,offset,sp;
  long symbols[nsymbols];
  double maxenergy,halfclock,scount,energy;
  int switchpoints[nsymbols*2*Symbolclocks];

  assert(symphase != NULL && samples != NULL);

  // Determine relative waveform transition points
  halfclock = (0.5 / Symbolclocks) * symbolsamples;
  scount = halfclock;

  // Do first offset at -1/2 symbol
  offset = -Symbolsamples/2;
  memset(symbols,0,sizeof(symbols));
  sp = 0;
  energy = 0;
  ind = firstsample + offset;
  assert(ind >= 0);
  for(i=0;i<nsymbols;i++){
    symbols[i] = 0;
    for(j=0; j<Symbolclocks; j++){
      // Integrate first half of clock cycle
      switchpoints[sp] = nearbyint(scount); // Remember for later adjustments
      scount += halfclock;
      for(;ind < switchpoints[sp] + firstsample + offset; ind++)
	symbols[i] -= samples[ind];

      sp++;
      // Integrate second half of clock cycle
      switchpoints[sp] = nearbyint(scount);
      scount += halfclock;
      for(; ind < switchpoints[sp] + firstsample + offset; ind++)
	symbols[i] += samples[ind];

      sp++;
    }
    energy += (long long)symbols[i] * symbols[i];
  }
  assert(sp <= sizeof(switchpoints)/sizeof(*switchpoints));
  // So far the only one is the best one
  maxenergy = energy;
  *symphase = offset;
	 
  // Now do the other offsets by updating the symbols already computed and recomputing energy
  for(offset++ ; offset < Symbolsamples/2; offset++){
    sp = 0;
    energy = 0;

    for(i=0;i<nsymbols;i++){
      for(j=0; j<Symbolclocks; j++){
	// Remove previous first sample
	if(sp == 0)
	  symbols[i] += samples[firstsample + offset -1];
	else
	  symbols[i] += samples[firstsample + offset + switchpoints[sp-1]-1];

	// Flip sign of sample in the middle from + to -
	symbols[i] -= 2 * samples[firstsample + offset + switchpoints[sp]-1];
	sp++;

	// Include one more sample to the right
	symbols[i] += samples[firstsample + offset + switchpoints[sp]-1];
	sp++;
      }
      energy += (long long)symbols[i] * symbols[i];
    }
    assert(sp <= sizeof(switchpoints)/sizeof(*switchpoints));
    if(energy > maxenergy){
      // Broke the old record
      maxenergy = energy;
      *symphase = offset;
    }
  }
  return maxenergy / nsymbols;
}
