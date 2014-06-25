// ISEE-3/ICE Viterbi decoder
// Phil Karn, KA9Q, June 2014

#define __USE_GNU   1
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
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
  void *vd;
  unsigned char vdsyms[2];
  int startup_delay;
  int decode_delay = 200;
  unsigned long long bits = 0;
  unsigned char oldsymbols[SYMBOLBUFSIZE];
  int symbols = 0;
  unsigned long long symerrs = 0;
  unsigned long long re_encoder = 0;
  int sync_sum,k;
  int status_interval = 1024;
  int sync_count = 0;
  int peak_inphase_sync,peak_outphase_sync;

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
  while((i = getopt(argc,argv,"d:pi:qF")) != EOF){
    switch(i){
    case 'F':
      Dontflip = 1; // Force phase
      break;
    case 'q':
      Quiet = 1;
      break;
    case 'p':
      symbols = 1; // start with opposite decoder phase
      break;
    case 'i':
      status_interval = atoi(optarg);
      break;
    case 'd':
      decode_delay = atoi(optarg);
      break;
    }
  }
  if(decode_delay < 24){
    fprintf(stderr,"%s: decoder delay too small, using 200\n",argv[0]);
    decode_delay = 200;
  } else if(decode_delay > 1024){
    fprintf(stderr,"%s: Warning; excessive decode delay; 1MB/bit needed\n",argv[0]);
  }
  startup_delay = decode_delay;

  vd = create_viterbi224(decode_delay+1); // Must be bigger than the traceback depth
  assert(vd != NULL);
  init_viterbi224(vd,0); // See what happens


  peak_inphase_sync = peak_outphase_sync = -1000000;

  while(1){
    int c;

    if((c = getchar()) == EOF)
      break;
    // Save for later re-encode comparison
    oldsymbols[symbols] = c;
    vdsyms[symbols % 2] = c;

    if(!Dontflip){ // Turn all this off if the user wants to keep the starting phase
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
      int s1,s2;
      update_viterbi224_blk(vd,vdsyms,1);

      // Chain back and decode a bit
      // This is far enough back for the paths to merge, we don't
      // need to waste time finding the best path
      // suppress the first decode_delay bits of output
      if(startup_delay == 0){
	bit = decodebit_viterbi224(vd,decode_delay,0); // Cannot be larger than value given to create_viterbi
	putchar(bit  ? '1' : '0');
	fflush(stdout);
	re_encoder = (re_encoder << 1) | bit;
      } else {
	startup_delay--;
      } 
      assert(POLY1 != POLY2);
      s1 = G1FLIP ^ parity(re_encoder & POLY1);
      s2 = G2FLIP ^ parity(re_encoder & POLY2);

      // Tally corrected channel symbols
#if 0
      {
	int c1,c2;
	c1 = oldsymbols[(SYMBOLBUFSIZE + symbols - 2*(decode_delay+K-2) - 1) % SYMBOLBUFSIZE];
	c2 = oldsymbols[(SYMBOLBUFSIZE + symbols - 2*(decode_delay+K-2)) % SYMBOLBUFSIZE];
	if((c1 > 128) ^ s1 || (c2 > 128) ^ s2){
	  fprintf(stderr,"%s: sym error bits %llu %d/%d   %d/%d\n",argv[0],bits,c1,s1,c2,s2);
	}
      }
#endif
      if(startup_delay == 0)
	// Compare re-encoded symbols with hard-sliced versions of those actually received
	symerrs += (s1 ^ (oldsymbols[(SYMBOLBUFSIZE + symbols - 2*(decode_delay+K-2) - 1) % SYMBOLBUFSIZE] > 128))
	  + (s2 ^ (oldsymbols[(SYMBOLBUFSIZE + symbols - 2*(decode_delay+K-2)) % SYMBOLBUFSIZE] > 128));
      

      if(!Quiet && status_interval != 0 && (++bits % status_interval) == 0){
	fprintf(stderr,"%s: bits %'llu; symerrs %'llu/%d %'.3lg%%\n",argv[0],
		bits,symerrs,2*status_interval,100.*symerrs/(2.*status_interval));
	symerrs = 0;
      }
    }
    symbols = (symbols + 1) % SYMBOLBUFSIZE;
  }
  exit(0);
}
