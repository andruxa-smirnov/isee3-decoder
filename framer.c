// Frame  ISSE-3 frames
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
#include <stdint.h>
#include <getopt.h>
#include <assert.h>
#include <locale.h>
#include "code.h"
#include "viterbi224.h"

#define FRAMEBITS 1024   // 1024 bits per minor frame
#define FRAMESYMBOLS (2*FRAMEBITS)
#define SYNCBITS 34      // Last 34 bits of encoded tail + sync are usable
#define SYMRATE 1024.467    // Nominal symrate for 512 bps mode

#define SYNCWORD 0x12fc819fbeLL

int Verbose;
double Symbolsamples;
double Samprate = 250000;   // Sample rate of PM demodulated symbols
double Symrate = SYMRATE;   // Symbol rate

double trial_demod(double *samples,int firstsample,double symbolsamples,int symbols);
void hms(int *hours, int *minutes, double *seconds,double t);
char *format_hms(double t);

// 34 bit encoded sync symbols

int sync_vector[] = {
  0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1,
  1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};


int main(int argc,char *argv[]){
  int i;
  char *locale;
  int k;
  unsigned long long frames = 1;
  unsigned long long decoded_data[16];
  char c;
  unsigned long long bits = 0;
  int bitrate = 512;

  if((locale = getenv("LANG")) != NULL)
    setlocale(LC_ALL,locale);
  else
    setlocale(LC_ALL,"en_US.utf8");

  memset(decoded_data,0,sizeof(decoded_data));
  while((i = getopt(argc,argv,"r:")) != EOF){
    switch(i){
    case 'r':
      bitrate = atoi(optarg);
      break;
    }
  }
  // Read viterbi decoder  output, accumulate into frames
  bits = 0;
  while(1){
    if((c = getchar()) == EOF)
      break;
    
    c = (c == '1') ? 1 : 0;
    // 1024-bit shift register with decoded bits
    for(k=0;;k++){
      decoded_data[k] <<= 1;
      if(k == 15)
	break;
      decoded_data[k] |= (decoded_data[k+1] >> 63);
    }
    decoded_data[15] |= c;
    
    if((decoded_data[15] & 0xffffffffffL) == SYNCWORD){
      int k;
      
      printf("Frame %'llu at bit %'llu (%s)\n",frames,bits,format_hms((double)bits/bitrate));
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
      frames++;
      putchar('\n');
      fflush(stdout);
    }
    bits++;
  }
  exit(0);
}

