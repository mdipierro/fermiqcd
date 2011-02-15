#include "assert.h"

#define MSK1 0xdfffffefU
#define MSK2 0xddfecb7fU
#define MSK3 0xbffaffffU
#define MSK4 0xbffffff6U
#define PARITY1 0x00000001U
#define PARITY2 0x00000000U
#define PARITY3 0x00000000U
#define PARITY4 0x13c9e684U

class mdp_prng_sfmt {
private:
  bool initialized;
  static const unsigned int MEXP=19937;
  static const unsigned int N=156;
  static const unsigned int N32=624;
  static const unsigned int POS1=122;
  static const unsigned int SL1=18;
  static const unsigned int SL2=1;
  static const unsigned int SR1=11;
  static const unsigned int SR2=1;
  unsigned int idx;  
  struct W128_T { unsigned int u[4]; };
  typedef struct W128_T w128_t;
  w128_t sfmt[N];
  unsigned int *psfmt32;
  void period_certification(void) {
    static unsigned int parity[4] = {PARITY1, PARITY2, PARITY3, PARITY4};
    unsigned int inner = 0;
    unsigned int i, j;
    unsigned int work;
    
    for (i = 0; i < 4; i++)
      inner ^= psfmt32[i] & parity[i];
    for (i = 16; i > 0; i >>= 1)
      inner ^= inner >> i;
    inner &= 1;
    /* check OK */
    if (inner == 1) { return; }
    /* check NG, and modification */
    for (i = 0; i < 4; i++) {
      work = 1;
      for (j = 0; j < 32; j++) {
	if ((work & parity[i]) != 0) {
	  psfmt32[i] ^= work;
	  return;
	}
	work = work << 1;
      }
    }
  }
  void rshift128(w128_t *out, w128_t const *in, int shift) {
    unsigned int mdp_int th, tl, oh, ol;
    
    th = ((unsigned int mdp_int)in->u[3] << 32) | ((unsigned int mdp_int)in->u[2]);
    tl = ((unsigned int mdp_int)in->u[1] << 32) | ((unsigned int mdp_int)in->u[0]);
    
    oh = th >> (shift * 8);
    ol = tl >> (shift * 8);
    ol |= th << (64 - shift * 8);
    out->u[1] = (unsigned int)(ol >> 32);
    out->u[0] = (unsigned int)ol;
    out->u[3] = (unsigned int)(oh >> 32);
    out->u[2] = (unsigned int)oh;
  }
  void lshift128(w128_t *out, w128_t const *in, int shift) {
    unsigned int mdp_int th, tl, oh, ol;

    th = ((unsigned int mdp_int)in->u[3] << 32) | ((unsigned int mdp_int)in->u[2]);
    tl = ((unsigned int mdp_int)in->u[1] << 32) | ((unsigned int mdp_int)in->u[0]);

    oh = th << (shift * 8);
    ol = tl << (shift * 8);
    oh |= tl >> (64 - shift * 8);
    out->u[1] = (unsigned int)(ol >> 32);
    out->u[0] = (unsigned int)ol;
    out->u[3] = (unsigned int)(oh >> 32);
    out->u[2] = (unsigned int)oh;
  }

  void gen_rand_all(void) {
    int i;
    w128_t *r1, *r2;
    
    r1 = &sfmt[N - 2];
    r2 = &sfmt[N - 1];
    for (i = 0; i < N - POS1; i++) {
      do_recursion(&sfmt[i], &sfmt[i], &sfmt[i + POS1], r1, r2);
      r1 = r2;
      r2 = &sfmt[i];
    }
    for (; i < N; i++) {
      do_recursion(&sfmt[i], &sfmt[i], &sfmt[i + POS1 - N], r1, r2);
      r1 = r2;
      r2 = &sfmt[i];
    }
  }
  
  void do_recursion(w128_t *r, w128_t *a, w128_t *b, w128_t *c,
		    w128_t *d) {
    w128_t x;
    w128_t y;
    
    lshift128(&x, a, SL2);
    rshift128(&y, c, SR2);
    r->u[0] = a->u[0] ^ x.u[0] ^ ((b->u[0] >> SR1) & MSK1) ^ y.u[0] 
      ^ (d->u[0] << SL1);
    r->u[1] = a->u[1] ^ x.u[1] ^ ((b->u[1] >> SR1) & MSK2) ^ y.u[1] 
      ^ (d->u[1] << SL1);
    r->u[2] = a->u[2] ^ x.u[2] ^ ((b->u[2] >> SR1) & MSK3) ^ y.u[2] 
      ^ (d->u[2] << SL1);
    r->u[3] = a->u[3] ^ x.u[3] ^ ((b->u[3] >> SR1) & MSK4) ^ y.u[3] 
      ^ (d->u[3] << SL1);
  }
  unsigned int gen_rand32() {
    unsigned int r;    
    assert(initialized);
    if (idx >= N32) {
      gen_rand_all();
      idx = 0;
    }
    r = psfmt32[idx++];
    return r;
  }

public:
  void initialize(unsigned int seed) {
    unsigned int i;
    psfmt32=(unsigned int *)&(sfmt[0].u[0]);
    psfmt32[0] = seed;
    for (i = 1; i < N32; i++) {
      psfmt32[i] = 1812433253UL * (psfmt32[i - 1] 
				   ^ (psfmt32[i - 1] >> 30))
	+ i;
    }
    idx = N32;
    period_certification();
    initialized = 1;
  }
  
  
  float plain() {
    unsigned int v=gen_rand32();
    return v * (1.0/4294967295.0); 
  }
};

