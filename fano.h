int fano(unsigned long *metric, unsigned long *cycles,
	unsigned char *data,const unsigned char *symbols,
	unsigned int nbits,int mettab[2][256],int delta,
	 unsigned long maxcycles,unsigned long long encstate);
int encode(unsigned char *symbols,const unsigned char *data,unsigned int nbytes,unsigned long long encstate);
void gen_met(int mettab[2][256],double signal,double noise,double bias,double scale);


