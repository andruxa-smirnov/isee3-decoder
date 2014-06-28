// Experimental Fano decoder for ISEE-3/ICE
// Phil Karn, KA9Q, June 2014

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#define __USE_GNU   1
#include <math.h>
#include <sys/types.h>
#include <string.h>
#include <errno.h>
#include <stdint.h>
#include <getopt.h>
#include <assert.h>
#include <locale.h>
#include "code.h"
#include "viterbi224.h"
#include "fano.h"
#include "timeformat.h"

#define FRAMEBITS 1024   // 1024 bits per minor frame
#define FRAMESYMBOLS (2*FRAMEBITS) // rate 1/2 code
#define SYNCBITS 34      // Last 34 bits of encoded tail + sync are usable
#define SYNCWORD 0x12fc819fbeLL // Last 5 bytes of every data frame

int Verbose;
double Symrate;
int Fano;
double Fano_scale;
int Fano_delta;
unsigned long Fano_maxcycles;

// 34 bit encoded sync symbols
int sync_vector[] = {
  0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1,
  1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

static inline int parity(unsigned long long x){
  return __builtin_parityll(x);
}


int main(int argc,char *argv[]){
  int i;
  char *locale;
  unsigned char symbols[2*FRAMESYMBOLS];
  void *vd = NULL;
  unsigned char data[FRAMEBITS/8];
  int sync_start;
  int symcount;
  int lock;
  long long total_symbols = 0;
  int mettab[2][256];
  long long frames = 1;

  if((locale = getenv("LANG")) != NULL)
    setlocale(LC_ALL,locale);
  else
    setlocale(LC_ALL,"en_US.utf8");

  memset(symbols,0,sizeof(symbols));

  Symrate = 1024;
  Fano_scale = 8;
  Fano_delta = 4 * Fano_scale;
  Fano_maxcycles = 100;
  while((i = getopt(argc,argv,"fr:s:m:d:")) != EOF){
    switch(i){
    case 'd':
      Fano_delta = atoi(optarg);
      break;
    case 'm':
      Fano_maxcycles = atol(optarg);
      break;
    case 's':
      Fano_scale = atof(optarg);
      break;
    case 'f':
      Fano = 1;
      break;
    case 'r':
      Symrate = atof(optarg);
      break;
    }
  }
  if(!Fano){
    printf("%s: using Viterbi decoding\n",argv[0]);
    vd = create_viterbi224(FRAMEBITS);
    assert(vd != NULL);
  } else {
    double noise_amp,sig_amp,total_amp;
    double est_esn0;

    // Ideally we'd have the independent signal and noise amplitudes, but we only
    // have the total amplitude as scaled by symdemod
    // Compute signal and noise amplitudes assuming operation at decoder threshold, Eb/No = 3dB
    total_amp = 100.; // Total amplitude set by symdemod
    est_esn0 = 1.0; // Es/No = 0 dB; Eb/No = 3 dB
    noise_amp = total_amp/sqrt(1 + 2*est_esn0);
    sig_amp = noise_amp * sqrt(2*est_esn0);

    printf("%s: using Fano decoding: delta %'d; scale %'.1lf; maxcycles %'lu; signal %.1lf; noise %.1lf\n",
	   argv[0],Fano_delta,Fano_scale,Fano_maxcycles,sig_amp,noise_amp);
    gen_met(mettab,sig_amp,noise_amp,0.5,Fano_scale);
  }
  sync_start = -1;
  symcount = 0;
  lock = 0;
  while(1){
    int c,i,decode_result;

    // Refill buffer with 1 frame + 1 sync length
    while(symcount < FRAMESYMBOLS+SYNCBITS){
      if((c = getchar()) == EOF)
	goto done;
      symbols[symcount++] = c;
    }
    if(!lock){ // If the last frame decoded we don't need to do this
      int record_sum;

      record_sum = -1000000;
      for(i=0;i<FRAMESYMBOLS;i++){
	// Run sync correlator over 1 frame
	int sync_sum,k;
	
	sync_sum = 0;
	for(k=0;k<SYNCBITS;k++){
	  int sym;
	  
	  sym = (int)symbols[i + k] - 128;
	  sync_sum += sync_vector[k] ? sym : -sym;
	}	  
	if(sync_sum > record_sum){
	  record_sum = sync_sum;
	  sync_start = i;
	}
      }
    }      
    // sync_start now points to the start of the sync of the frame before the current one
    // Make sure we have one frame of symbols past the sync
    while(symcount < sync_start + FRAMESYMBOLS+SYNCBITS){
      if((c = getchar()) == EOF)
	goto done;
      symbols[symcount++] = c;
    }
    if(!Fano){
      // Attempt a decode starting after the sync
      init_viterbi224(vd,SYNCWORD & 0xffffff); // Known starting state after sync
      update_viterbi224_blk(vd,&symbols[sync_start+SYNCBITS],FRAMEBITS);
      chainback_viterbi224(vd,data,FRAMEBITS,SYNCWORD & 0xffffff);
      decode_result = FRAMEBITS; // Viterbi always "succeeds"
    } else {
      // Try a Fano decode
      unsigned long metric;
      unsigned long cycles;

      memset(data,0,sizeof(data)); // Wipe out previous data
      decode_result = fano(&metric,&cycles,data,&symbols[sync_start+SYNCBITS],FRAMEBITS,mettab,Fano_delta,100,
	       SYNCWORD & 0xffffff,SYNCWORD & 0xffffff);
#if 0
      printf("Fano returns %d metric %'ld cycles %'ld\n",decode_result,metric,cycles);
#endif
    }
    // If the decoder succeeded, see if the decoded frame ends with the 5-byte sync sequence
    // The last 23 bits from the Viterbi decoder will always be there since we forced them,
    // but the 17 bits before them are actual data
    // A decoder failure always indicates non-lock
    lock = 0;
    if(decode_result == FRAMEBITS){
      long long lastword;
      int i;

      lastword = 0;
      for(i=123;i<128;i++)
	lastword = (lastword << 8) | data[i];
      
      if(lastword == SYNCWORD)
	lock = 1;      // Good frame
    }
    // Dump data
    // This will dump 0's on a failed Fano frame. We could get clever and show X's or something.
    {
      long long start_time;

      start_time = total_symbols + sync_start + SYNCBITS;
      printf("Frame %'llu at symbol %'lld (%s) - %s:\n",
	     frames,start_time,format_hms(start_time/Symrate),
	     lock ? "good" : "bad");
      for(i=0;i<FRAMEBITS/8;i++){
	printf("%02x",data[i]);
	if((i % 16) == 15)
	  putchar('\n');
	else
	  putchar(' ');
      }
      frames++;
      putchar('\n'); // Blank line between frames
      fflush(stdout);
    }
    // Purge frame we just decoded
    {
      int adjust;

      adjust = sync_start + FRAMESYMBOLS;
      assert(symcount >= adjust);
      memmove(symbols,&symbols[adjust],
	      sizeof(*symbols)*(symcount - adjust));
      symcount -= adjust;
      sync_start = 0;
      total_symbols += adjust;
      // The symbol buffer now starts with the 34-bit sync sequence of the frame we just decoded
    }
  }
done:;
  if(vd != NULL)
    delete_viterbi224(vd);

  exit(0);
}

