#include "stdio.h"

class mdp_complex {
public:
  float re, im;
  float &real() {
    return re;
  }
  float &imag() {
    return im;
  } 
};

#define ALIGN16 __attribute__ ((aligned (16)))
#define ALIGN64 __attribute__ ((aligned (64)))
#define _ASM __asm__ __volatile__

typedef struct { float c1,c2,c3,c4;  } _sse_float  ALIGN16;  
typedef struct { _sse_float c1,c2,c3;} _sse_vector ALIGN16;
typedef struct { int c1,c2,c3,c4;}     _sse_int    ALIGN16;
typedef struct { double c1,c2; }       _sse_double ALIGN16;

typedef struct {mdp_complex c11,c12,c13,c21,c22,c23,c31,c32,c33; } _sse_su3;
typedef struct {mdp_complex c1,c2,c3; } _sse_su3_vector;
typedef struct {_sse_su3_vector c1,c2,c3,c4; } _sse_spinor;

static _sse_float _sse_float_sgn12 __attribute__ ((unused)) = {-1.0f,-1.0f,1.0f,1.0f};
static _sse_float _sse_float_sgn13 __attribute__ ((unused)) = {-1.0f,1.0f,-1.0f,1.0f};
static _sse_float _sse_float_sgn14 __attribute__ ((unused)) = {-1.0f,1.0f,1.0f,-1.0f};
static _sse_float _sse_float_sgn23 __attribute__ ((unused)) = {1.0f,-1.0f,-1.0f,1.0f};
static _sse_float _sse_float_sgn24 __attribute__ ((unused)) = {1.0f,-1.0f,1.0f,-1.0f};
static _sse_float _sse_float_sgn34 __attribute__ ((unused)) = {1.0f,1.0f,-1.0f,-1.0f};
static _sse_int   _sse_double_sgn  __attribute__ ((unused)) = {0x0,0x80000000,0x0,0x0};
static _sse_int   _sse_double_sgn2 __attribute__ ((unused)) = {0x0,0x0,0x0,0x80000000};

#define _sse_float_su3_multiply(u) { \
 _ASM ("movss %0, %%xmm3 \n\t" \
      "movss %1, %%xmm6 \n\t" \
      "movss %2, %%xmm4 \n\t" \
      "movss %3, %%xmm7 \n\t" \
      "movss %4, %%xmm5 \n\t" \
      "shufps $0x0, %%xmm3, %%xmm3 \n\t" \
      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
      "shufps $0x0, %%xmm4, %%xmm4 \n\t" \
      "mulps %%xmm0, %%xmm3 \n\t" \
      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
      "mulps %%xmm1, %%xmm6 \n\t" \
      "shufps $0x0, %%xmm5, %%xmm5 \n\t" \
      "mulps %%xmm0, %%xmm4 \n\t" \
      "addps %%xmm6, %%xmm3 \n\t" \
      "mulps %%xmm2, %%xmm7 \n\t" \
      "mulps %%xmm0, %%xmm5 \n\t" \
      "addps %%xmm7, %%xmm4 \n\t" \
      "movss %5, %%xmm6 \n\t" \
      "movss %6, %%xmm7 \n\t" \
      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
      "mulps %%xmm1, %%xmm6 \n\t" \
      "mulps %%xmm2, %%xmm7 \n\t" \
      "addps %%xmm6, %%xmm5 \n\t" \
      "addps %%xmm7, %%xmm3 \n\t" \
      "movss %7, %%xmm6 \n\t" \
      "movss %8, %%xmm7 \n\t" \
      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
      "mulps %%xmm1, %%xmm6 \n\t" \
      "mulps %%xmm2, %%xmm7 \n\t" \
      "addps %%xmm6, %%xmm4 \n\t" \
      "addps %%xmm7, %%xmm5" \
      : \
      : \
      "m" ((u).c11.real()), \
      "m" ((u).c12.real()), \
      "m" ((u).c21.real()), \
      "m" ((u).c23.real()), \
      "m" ((u).c31.real()), \
      "m" ((u).c32.real()), \
      "m" ((u).c13.real()), \
      "m" ((u).c22.real()), \
      "m" ((u).c33.real())); \
_ASM ("movss %0, %%xmm6 \n\t" \
      "movss %1, %%xmm7 \n\t" \
      "shufps $0xb1, %%xmm0, %%xmm0 \n\t" \
      "shufps $0xb1, %%xmm1, %%xmm1 \n\t" \
      "shufps $0xb1, %%xmm2, %%xmm2 \n\t" \
      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
      "mulps %9, %%xmm0 \n\t" \
      "mulps %9, %%xmm1 \n\t" \
      "mulps %9, %%xmm2 \n\t" \
      "mulps %%xmm0, %%xmm6 \n\t" \
      "mulps %%xmm1, %%xmm7 \n\t" \
      "addps %%xmm6, %%xmm3 \n\t" \
      "addps %%xmm7, %%xmm4 \n\t" \
      "movss %2, %%xmm6 \n\t" \
      "movss %3, %%xmm7 \n\t" \
      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
      "mulps %%xmm2, %%xmm6 \n\t" \
      "mulps %%xmm0, %%xmm7 \n\t" \
      "addps %%xmm6, %%xmm5 \n\t" \
      "addps %%xmm7, %%xmm4 \n\t" \
      "movss %4, %%xmm6 \n\t" \
      "movss %5, %%xmm7 \n\t" \
      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
      "mulps %%xmm1, %%xmm6 \n\t" \
      "mulps %%xmm0, %%xmm7 \n\t" \
      "addps %%xmm6, %%xmm3 \n\t" \
      "addps %%xmm7, %%xmm5 \n\t" \
      "movss %6, %%xmm0 \n\t" \
      "movss %7, %%xmm6 \n\t" \
      "movss %8, %%xmm7 \n\t" \
      "shufps $0x0, %%xmm0, %%xmm0 \n\t" \
      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
      "mulps %%xmm2, %%xmm0 \n\t" \
      "mulps %%xmm1, %%xmm6 \n\t" \
      "mulps %%xmm2, %%xmm7 \n\t" \
      "addps %%xmm0, %%xmm3 \n\t" \
      "addps %%xmm6, %%xmm5 \n\t" \
      "addps %%xmm7, %%xmm4" \
      : \
      : \
      "m" ((u).c11.imag()), \
      "m" ((u).c22.imag()), \
      "m" ((u).c33.imag()), \
      "m" ((u).c21.imag()), \
      "m" ((u).c12.imag()), \
      "m" ((u).c31.imag()), \
      "m" ((u).c13.imag()), \
      "m" ((u).c32.imag()), \
      "m" ((u).c23.imag()), \
      "m" (_sse_float_sgn13)); }

// //////////////////////////////////////////////////////////////////////////  
// Multiplies a pair sl,sh of su3 vectors with an su3 matrix u^dagger, 
// assuming sl and sh are in the low and high words of xmm0,xmm1,xmm2
//
// On output the result is in xmm3,xmm4,xmm5 and the registers 
// xmm0,xmm1,xmm2 are changed
// //////////////////////////////////////////////////////////////////////////  

#define _sse_float_su3_inverse_multiply(u) { \
_ASM ("movss %0, %%xmm3 \n\t" \
      "movss %1, %%xmm6 \n\t" \
      "movss %2, %%xmm4 \n\t" \
      "movss %3, %%xmm7 \n\t" \
      "movss %4, %%xmm5 \n\t" \
      "shufps $0x0, %%xmm3, %%xmm3 \n\t" \
      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
      "shufps $0x0, %%xmm4, %%xmm4 \n\t" \
      "mulps %%xmm0, %%xmm3 \n\t" \
      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
      "mulps %%xmm1, %%xmm6 \n\t" \
      "shufps $0x0, %%xmm5, %%xmm5 \n\t" \
      "mulps %%xmm0, %%xmm4 \n\t" \
      "addps %%xmm6, %%xmm3 \n\t" \
      "mulps %%xmm2, %%xmm7 \n\t" \
      "mulps %%xmm0, %%xmm5 \n\t" \
      "addps %%xmm7, %%xmm4 \n\t" \
      "movss %5, %%xmm6 \n\t" \
      "movss %6, %%xmm7 \n\t" \
      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
      "mulps %%xmm1, %%xmm6 \n\t" \
      "mulps %%xmm2, %%xmm7 \n\t" \
      "addps %%xmm6, %%xmm5 \n\t" \
      "addps %%xmm7, %%xmm3 \n\t" \
      "movss %7, %%xmm6 \n\t" \
      "movss %8, %%xmm7 \n\t" \
      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
      "mulps %%xmm1, %%xmm6 \n\t" \
      "mulps %%xmm2, %%xmm7 \n\t" \
      "addps %%xmm6, %%xmm4 \n\t" \
      "addps %%xmm7, %%xmm5" \
      : \
      : \
      "m" ((u).c11.real()), \
      "m" ((u).c21.real()), \
      "m" ((u).c12.real()), \
      "m" ((u).c32.real()), \
      "m" ((u).c13.real()), \
      "m" ((u).c23.real()), \
      "m" ((u).c31.real()), \
      "m" ((u).c22.real()), \
      "m" ((u).c33.real())); \
_ASM ("movss %0, %%xmm6 \n\t" \
      "movss %1, %%xmm7 \n\t" \
      "shufps $0xb1, %%xmm0, %%xmm0 \n\t" \
      "shufps $0xb1, %%xmm1, %%xmm1 \n\t" \
      "shufps $0xb1, %%xmm2, %%xmm2 \n\t" \
      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
      "mulps %9, %%xmm0 \n\t" \
      "mulps %9, %%xmm1 \n\t" \
      "mulps %9, %%xmm2 \n\t" \
      "mulps %%xmm0, %%xmm6 \n\t" \
      "mulps %%xmm1, %%xmm7 \n\t" \
      "addps %%xmm6, %%xmm3 \n\t" \
      "addps %%xmm7, %%xmm4 \n\t" \
      "movss %2, %%xmm6 \n\t" \
      "movss %3, %%xmm7 \n\t" \
      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
      "mulps %%xmm2, %%xmm6 \n\t" \
      "mulps %%xmm0, %%xmm7 \n\t" \
      "addps %%xmm6, %%xmm5 \n\t" \
      "addps %%xmm7, %%xmm4 \n\t" \
      "movss %4, %%xmm6 \n\t" \
      "movss %5, %%xmm7 \n\t" \
      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
      "mulps %%xmm1, %%xmm6 \n\t" \
      "mulps %%xmm0, %%xmm7 \n\t" \
      "addps %%xmm6, %%xmm3 \n\t" \
      "addps %%xmm7, %%xmm5 \n\t" \
      "movss %6, %%xmm0 \n\t" \
      "movss %7, %%xmm6 \n\t" \
      "movss %8, %%xmm7 \n\t" \
      "shufps $0x0, %%xmm0, %%xmm0 \n\t" \
      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
      "mulps %%xmm2, %%xmm0 \n\t" \
      "mulps %%xmm1, %%xmm6 \n\t" \
      "mulps %%xmm2, %%xmm7 \n\t" \
      "addps %%xmm0, %%xmm3 \n\t" \
      "addps %%xmm6, %%xmm5 \n\t" \
      "addps %%xmm7, %%xmm4" \
      : \
      : \
      "m" ((u).c11.imag()), \
      "m" ((u).c22.imag()), \
      "m" ((u).c33.imag()), \
      "m" ((u).c12.imag()), \
      "m" ((u).c21.imag()), \
      "m" ((u).c13.imag()), \
      "m" ((u).c31.imag()), \
      "m" ((u).c23.imag()), \
      "m" ((u).c32.imag()), \
      "m" (_sse_float_sgn24)); }


int main() {
  _sse_su3 u;
  _sse_float_su3_multiply(u);
  _sse_float_su3_multiply(u);
  return 0;
}
