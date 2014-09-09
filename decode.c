// Fano/Viterbi decoder for ISEE-3/ICE
// Phil Karn, KA9Q, July 2014

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
#define SYNCWORD 0x12fc819fbeULL // Last 5 bytes of every data frame

int Verbose;
int Viterbi_enabled;
double Symrate;
int Fano_enabled;
double Fano_scale;
int Fano_delta;
unsigned long Fano_maxcycles;
int No_bad_frames;
int Persistent;

// 34 bit encoded sync symbols
int sync_vector[] = {
  0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1,
  1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

int main(int argc,char *argv[]){
  int i;
  char *locale;
  unsigned char symbols[2*FRAMESYMBOLS];
  void *vd = NULL;
  unsigned char data[FRAMEBITS/8];
  int sync_start;
  int symcount;
  int lock;
  unsigned long long total_symbols = 0;
  int mettab[2][256];
  long long frames = 1;
  enum { NONE=0, VITERBI, FANO } decoder;

  // I've always hated large numbers displayed without commas. You go blind counting digits.
  // The format strings include the (') flag so setting the locale will cause all large
  // numbers to be displayed with the local thousands separator.
  // In the US (default), this is a comma (,)
  if((locale = getenv("LANG")) != NULL)
    setlocale(LC_ALL,locale);
  else
    setlocale(LC_ALL,"en_US.utf8");

  // Set defaults
  No_bad_frames = 0;
  Symrate = 1024;
  Fano_scale = 8;
  Fano_maxcycles = 100;
  Fano_delta = 4 * Fano_scale;
  Viterbi_enabled = 1;
  Fano_enabled = 1;
  Persistent = 0;

  while((i = getopt(argc,argv,"nFVvr:s:m:d:p")) != EOF){
    switch(i){
    case 'p':
      Persistent = 1;
      break;
    case 'n':
      No_bad_frames = 1;
      break;
    case 'F':
      Viterbi_enabled = 0;
      break;
    case 'V':
      Fano_enabled = 0;
      break;
    case 'v':
      Verbose++;
      break;
    case 'r':
      Symrate = atof(optarg);
      break;
    case 's':
      Fano_scale = atof(optarg);
      break;
    case 'm':
      Fano_maxcycles = atol(optarg);
      break;
    case 'd':
      Fano_delta = atoi(optarg);
      break;
    default:
      printf("usage: %s [-F] [-V] [-v] [-r symrate] [-s fano_scale] [-m fano_maxcycles] [-d fano_delta]\n",
	     argv[0]);
    }
  }
  printf("%s: Fano %s; Viterbi %s\n",argv[0],
	 Fano_enabled ? "enabled" : "disabled",
	 Viterbi_enabled ? "enabled" : "disabled");

  if(No_bad_frames)
    printf("%s: Not displaying bad frames\n",argv[0]);

  if(!Fano_enabled && !Viterbi_enabled){
    printf("%s: Specify only one of -F or -V\n",argv[0]);
    exit(1);
  }
  if(Fano_enabled){
    // Set up Fano decoder
    double noise_amp,sig_amp,total_amp;
    double est_esn0;
    
    // Ideally we'd have the independent signal and noise amplitudes, but we only
    // have the total amplitude as scaled by symdemod
    // Compute signal and noise amplitudes assuming operation at decoder threshold, Eb/No = 3dB
    total_amp = 100.; // Total amplitude set by symdemod
    est_esn0 = 1.0; // Es/No = 0 dB; Eb/No = 3 dB
    noise_amp = total_amp/sqrt(1 + 2*est_esn0);
    sig_amp = noise_amp * sqrt(2*est_esn0);
    
    printf("%s: Fano decoder params: delta %'d; scale %'.1lf; maxcycles %'lu; signal %.1lf; noise %.1lf\n",
	   argv[0],Fano_delta,Fano_scale,Fano_maxcycles,sig_amp,noise_amp);
    gen_met(mettab,sig_amp,noise_amp,0.5,Fano_scale);
  } else if(Viterbi_enabled){
    // Only Viterbi is enabled, set up decoder just once since we'll use its memory constantly
    // When we do Fano with Viterbi fallback, the memory is allocated only when actually needed
    if((vd = create_viterbi224(FRAMEBITS)) == NULL){
      printf("%s: Fano decoding disabled but cannot alloc %d * 1 MB + 2 * 16 MB = %.1lf GB RAM for Viterbi decoder!\nTry -F for Fano only or the default of Fano with Viterbi fallback to allocate memory only when the Viterbi decoder is actually used.\n",
	     argv[0],FRAMEBITS,(FRAMEBITS+32)/1024.);
      exit(2);
    }
  }
  sync_start = -1;
  symcount = 0;
  lock = 0;
  memset(symbols,0,sizeof(symbols));
  while(1){
    int decode_result;

    // Refill buffer with 1 frame + 1 sync length
    if(symcount < FRAMESYMBOLS+SYNCBITS){
      int n,cnt;

      cnt = FRAMESYMBOLS+SYNCBITS-symcount;
      n = fread(&symbols[symcount],sizeof(*symbols),cnt,stdin);
      if(n < cnt)
	goto done;
      symcount += n;
    }
    if(!lock){ // If the last frame decoded we don't need to do this; its sync should be right at 0.
      int record_sum,i;

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
      // sync_start now points to the start of the sync of the frame before the current one
      // Make sure we have one frame of symbols past the sync
      if(symcount < sync_start + FRAMESYMBOLS+SYNCBITS){
	int n,cnt;

	cnt = sync_start + FRAMESYMBOLS+SYNCBITS-symcount;
	n = fread(&symbols[symcount],sizeof(*symbols),cnt,stdin);
	if(n < cnt)
	  goto done;
	symcount += n;
      }
    }
    decoder = NONE;
    decode_result = 0;
    if(Fano_enabled){
      // Try a Fano decode
      unsigned long metric,cycles;

      decoder = FANO;      
      memset(data,0,sizeof(data)); // Wipe out previous data
      decode_result = fano(&metric,&cycles,data,&symbols[sync_start+SYNCBITS],FRAMEBITS,mettab,Fano_delta,100,
			   SYNCWORD & 0xffffff,SYNCWORD & 0xffffff);

#if 0
      printf("%s: Fano returns %d; metric %'ld; cycles %'ld\n",argv[0],decode_result,metric,cycles);
#endif
    }
    if(Viterbi_enabled && (!Fano_enabled || ((Persistent || lock) && decode_result != FRAMEBITS))){
      // Decode with Viterbi if:
      // 1. Viterbi is enabled AND either
      // 2. Fano is disabled OR
      // 3. Fano tried and failed, AND the previous frame decoded with either Fano or Viterbi, OR
      // 4  the -p (persistent) option is selected
      // Alloc memory for decoder if not already done
      if(vd == NULL && (vd = create_viterbi224(FRAMEBITS)) == NULL){
	printf("%s: Cannot malloc %d * 1 MB + 2 * 16 MB = %.1lf GB RAM for Viterbi decoder! If this repeats, try Fano only (-F)\n",
	       argv[0],FRAMEBITS,(FRAMEBITS+32)/1024.);
      } else {
	init_viterbi224(vd,SYNCWORD & 0xffffff); // Known starting state after sync
	update_viterbi224_blk(vd,&symbols[sync_start+SYNCBITS],FRAMEBITS);
	chainback_viterbi224(vd,data,FRAMEBITS,SYNCWORD & 0xffffff);
	decoder = VITERBI;
	decode_result = FRAMEBITS;
	if(Fano_enabled){
	  // Don't hold onto all that RAM if we only use Viterbi as a fallback to Fano
	  delete_viterbi224(vd);
	  vd = NULL;
	}
      }
    }
    // If a decoder was tried and thinks it succeeded (Viterbi always does),
    // see if the decoded frame ends with the 5-byte sync sequence.
    // The last 23 bits will always be there since we force them in the decoders
    // but the 17 bits before them are actual data.
    lock = 0;
    if(decode_result == FRAMEBITS){
      unsigned long long lastword;
      int i;

      lastword = 0;
      for(i=123;i<128;i++)
	lastword = (lastword << 8) | data[i];
      
      if(lastword == SYNCWORD)
	lock = 1;      // Good frame
    }
    // Dump data
    if(lock || !No_bad_frames){
      unsigned long long start_symbol;
      int i;

      start_symbol = total_symbols + sync_start + SYNCBITS;
      printf("Frame %'llu at symbol %'llu (%s) with %s %s\n",
	     frames,start_symbol,format_hms(start_symbol/Symrate),
	     decoder == VITERBI ? "Viterbi" : decoder == FANO ? "Fano" : "None",
	     !lock ? "(bad)" : "");
      for(i=0;i<FRAMEBITS/8;i++){
	printf("%02x",data[i]);
	if((i % 16) == 15)
	  putchar('\n');
	else
	  putchar(' ');
      }
      putchar('\n'); // Blank line between frames
      fflush(stdout);
    }
    frames++; // Count good and bad frames
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
  if(vd != NULL){
    delete_viterbi224(vd);
    vd = NULL;
  }
  exit(0);
}
