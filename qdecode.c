// Experimental - QLI decoder, great for clean symbol stream
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

#define FRAMEBITS 1024   // 1024 bits per minor frame
#define FRAMESYMBOLS (2*FRAMEBITS)
#define SYNCBITS 34      // Last 34 bits of encoded tail + sync are usable

#define SYNCWORD 0x12fc819fbeLL

#define SYMBOLBUFSIZE 4096 // bigger than it has to be. Must be greater than decoder delay

int Verbose;
int Quiet;
int Dontflip;

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
  unsigned char vdsyms[2];
  unsigned char oldsymbols[SYMBOLBUFSIZE];
  int symbols = 0;
  int counts[1000];
  int k;
  int sync_sum;
  int sync_count = 0;
  int peak_inphase_sync;
  int peak_outphase_sync;


  memset(counts,0,sizeof(counts));
  for(i=0;i<SYMBOLBUFSIZE;i += 2){
    oldsymbols[i] = G1FLIP ? 255 : 0;
    oldsymbols[i+1] = G2FLIP ? 255 : 0;

  }

  if((locale = getenv("LANG")) != NULL)
    setlocale(LC_ALL,locale);
  else
    setlocale(LC_ALL,"en_US.utf8");

  symbols = 0; // low bit controls decoder symbol phasing
  Dontflip = 0;
  while((i = getopt(argc,argv,"Fpq")) != EOF){
    switch(i){
    case 'F':
      Dontflip = 1; // Keep starting phase
      break;
    case 'q':
      Quiet = 1;
      break;
    case 'p':
      symbols = 1; // start with opposite decoder phase
      break;
    }
  }
  peak_inphase_sync = peak_outphase_sync = -1000000;

  while(1){
    int c;

    if((c = getchar()) == EOF)
      break;
    // Save for later re-encode comparison
    oldsymbols[symbols] = c;
    vdsyms[symbols % 2] = c;

    if(!Dontflip){
      sync_sum = 0;
      
      for(k=0;k<34;k++){
	if(sync_vector[k])
	  sync_sum += oldsymbols[(SYMBOLBUFSIZE + symbols + k - 33) % SYMBOLBUFSIZE] - 128;
	else
	  sync_sum -= oldsymbols[(SYMBOLBUFSIZE + symbols + k - 33) % SYMBOLBUFSIZE] - 128;
      }	  
      
      if((symbols % 2) == 0){
	if(sync_sum > peak_outphase_sync){
	  peak_outphase_sync = sync_sum;
	}
      } else {
	if(sync_sum > peak_inphase_sync){
	  peak_inphase_sync = sync_sum;
	}	
	if(++sync_count >= FRAMESYMBOLS){
	  // See who had the strongest sync last frame
	  sync_count = 0;
	  if(peak_outphase_sync > peak_inphase_sync){
	    // Flip phase
	    if(!Quiet)
	      fprintf(stderr,"%s: flipping phase\n",argv[0]);
	    if((symbols % 2) == 0)
	      symbols++;
	    else
	      symbols--;
	  }
	  peak_inphase_sync = peak_outphase_sync = -1000000;
	}
      }	
    }
    if((symbols % 2) == 1){
      int bit;
      bit = (vdsyms[0] > 128) ^ (vdsyms[1] > 128) ^ 1;
      putchar(bit  ? '1' : '0');
      fflush(stdout);
    }
    symbols = (symbols + 1) % SYMBOLBUFSIZE;
  }
  exit(0);
}
