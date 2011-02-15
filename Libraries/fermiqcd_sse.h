/////////////////////////////////////////////////////////////////
/// @file fermiqcd_sse.h
/// @version 2009-12-21
/// @author Martin Luesher and Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Basic actions for Wilson Fermions optimized in assembler
///
//////////////////////////////////////////////////////////////////

#pragma fPIC

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

// //////////////////////////////////////////////////////////////////////////
// Cache manipulation macros (float)
// //////////////////////////////////////////////////////////////////////////


#define _sse_float_prefetch_spinor(addr) \
_ASM ("prefetcht0 %0 \n\t" \
      "prefetcht0 %1" \
      : \
      : \
      "m" (*(((char*)(((unsigned int)(addr))&~0x7f)))), \
      "m" (*(((char*)(((unsigned int)(addr))&~0x7f))+128)))
     
#define _sse_float_prefetch_su3(addr) \
_ASM ("prefetcht0 %0 \n\t" \
      "prefetcht0 %1" \
      : \
      : \
      "m" (*(((char*)(((unsigned int)(addr))&~0x7f)))), \
      "m" (*(((char*)(((unsigned int)(addr))&~0x7f))+128)))

// //////////////////////////////////////////////////////////////////////////
//
// Macros for su3 vectors 
//
// Most of these macros operate on pairs of su3 vectors that are stored
// in the low and high words of xmm0,xmm1,xmm2 or xmm3,xmm4,xmm5. For example,
//
// xmm0 -> sl.c1.real(),sl.c1.imag(),sh.c1.real(),sh.c1.imag()
// xmm1 -> sl.c2.real(),sl.c2.imag(),sh.c2.real(),sh.c2.imag()
// xmm2 -> sl.c3.real(),sl.c3.imag(),sh.c3.real(),sh.c3.imag()
//
// (where sl and sh are of type su3_vector). This can also be interpreted as
// an _sse_vector s that is stored in these registers according to
//
// xmm0 -> s.c1.c1,s.c1.c2,s.c1.c3,s.c1.c4
// xmm1 -> s.c2.c1,s.c2.c2,s.c2.c3,s.c2.c4
// xmm2 -> s.c3.c1,s.c3.c2,s.c3.c3,s.c3.c4
//
// The load and store macros can be used to move data in either format
// from and to the xmm registers
//
// //////////////////////////////////////////////////////////////////////////

// //////////////////////////////////////////////////////////////////////////
//
// Operations for SU3 color linear algebra used in mdp_matrix (float)
//
// //////////////////////////////////////////////////////////////////////////

// //////////////////////////////////////////////////////////////////////////
// Loads two su3 vectors sl and sh to the low and high words of xmm0,xmm1,xmm2
// //////////////////////////////////////////////////////////////////////////

#define _sse_float_pair_load(sl,sh) \
_ASM ("movlps %0, %%xmm0 \n\t" \
      "movlps %1, %%xmm1 \n\t" \
      "movlps %2, %%xmm2 \n\t" \
      "movhps %3, %%xmm0 \n\t" \
      "movhps %4, %%xmm1 \n\t" \
      "movhps %5, %%xmm2 " \
       : \
       : \
       "m" ((sl).c1), \
       "m" ((sl).c2), \
       "m" ((sl).c3), \
       "m" ((sh).c1), \
       "m" ((sh).c2), \
       "m" ((sh).c3))

// //////////////////////////////////////////////////////////////////////////
// Loads two su3 vectors sl and sh to the low and high words of xmm3,xmm4,xmm5
// //////////////////////////////////////////////////////////////////////////  

#define _sse_float_pair_load_up(sl,sh) \
_ASM ("movlps %0, %%xmm3 \n\t" \
      "movlps %1, %%xmm4 \n\t" \
      "movlps %2, %%xmm5 \n\t" \
      "movhps %3, %%xmm3 \n\t" \
      "movhps %4, %%xmm4 \n\t" \
      "movhps %5, %%xmm5" \
      : \
      : \
      "m" ((sl).c1), \
      "m" ((sl).c2), \
      "m" ((sl).c3), \
      "m" ((sh).c1), \
      "m" ((sh).c2), \
      "m" ((sh).c3))


// //////////////////////////////////////////////////////////////////////////  
// Stores the low and high words of xmm0,xmm1,xmm2 to the su3 vectors rl and rh
// //////////////////////////////////////////////////////////////////////////  

#define _sse_float_pair_store(rl,rh) \
_ASM ("movlps %%xmm0, %0 \n\t" \
      "movlps %%xmm1, %1 \n\t" \
      "movlps %%xmm2, %2 \n\t" \
      "movhps %%xmm0, %3 \n\t" \
      "movhps %%xmm1, %4 \n\t" \
      "movhps %%xmm2, %5" \
      : \
      "=m" ((rl).c1), \
      "=m" ((rl).c2), \
      "=m" ((rl).c3), \
      "=m" ((rh).c1), \
      "=m" ((rh).c2), \
      "=m" ((rh).c3))

// //////////////////////////////////////////////////////////////////////////  
// Stores the low and high words of xmm3,xmm4,xmm5 to the su3 vectors rl and rh
// //////////////////////////////////////////////////////////////////////////  

#define _sse_float_pair_store_up(rl,rh) \
_ASM ("movlps %%xmm3, %0 \n\t" \
      "movlps %%xmm4, %1 \n\t" \
      "movlps %%xmm5, %2 \n\t" \
      "movhps %%xmm3, %3 \n\t" \
      "movhps %%xmm4, %4 \n\t" \
      "movhps %%xmm5, %5" \
      : \
      "=m" ((rl).c1), \
      "=m" ((rl).c2), \
      "=m" ((rl).c3), \
      "=m" ((rh).c1), \
      "=m" ((rh).c2), \
      "=m" ((rh).c3))

// //////////////////////////////////////////////////////////////////////////  
// Loads the components s.c1,s.c2,s.c3 of an _sse_float_vector s to xmm0,xmm1,xmm2
// //////////////////////////////////////////////////////////////////////////  

#define _sse_float_vector_load(s) \
_ASM ("movaps %0, %%xmm0 \n\t" \
      "movaps %1, %%xmm1 \n\t" \
      "movaps %2, %%xmm2" \
      : \
      : \
      "m" ((s).c1), \
      "m" ((s).c2), \
      "m" ((s).c3))

#define _sse_float_vector_load_up(s) \
_ASM ("movaps %0, %%xmm3 \n\t" \
      "movaps %1, %%xmm4 \n\t" \
      "movaps %2, %%xmm5" \
      : \
      : \
      "m" ((s).c1), \
      "m" ((s).c2), \
      "m" ((s).c3))
// //////////////////////////////////////////////////////////////////////////  
// Stores xmm0,xmm1,xmm2 to the components r.c1,r.c2,r.c3 of an _sse_float_vector r 
// //////////////////////////////////////////////////////////////////////////  

#define _sse_float_vector_store(r) \
_ASM ("movaps %%xmm0, %0 \n\t" \
      "movaps %%xmm1, %1 \n\t" \
      "movaps %%xmm2, %2" \
      : \
      "=m" ((r).c1), \
      "=m" ((r).c2), \
      "=m" ((r).c3))

// //////////////////////////////////////////////////////////////////////////  
// Multiplies xmm0,xmm1,xmm2 with a constant _sse_float c
// //////////////////////////////////////////////////////////////////////////  

#define _sse_float_vector_mul(c) \
_ASM ("mulps %0, %%xmm0 \n\t" \
      "mulps %0, %%xmm1 \n\t" \
      "mulps %0, %%xmm2" \
      : \
      : \
      "m" (c))

// //////////////////////////////////////////////////////////////////////////  
// multiplies xmm0, xmm1, xmm2 for complex=a*(b+I)
// (written by Massimo Di Pierro)
// //////////////////////////////////////////////////////////////////////////  
/* deprecated

#define _sse_float_vector_mulc(a,b) \
_ASM ("mulps %0, %%xmm0 \n\t" \
      "mulps %0, %%xmm1 \n\t" \
      "mulps %0, %%xmm2 \n\t" \
      "movaps %%xmm0, %%xmm3 \n\t" \
      "movaps %%xmm1, %%xmm4 \n\t" \
      "movaps %%xmm2, %%xmm5 \n\t" \
      "mulps %1, %%xmm0 \n\t" \
      "mulps %1, %%xmm1 \n\t" \
      "mulps %1, %%xmm2 \n\t" \
      "shufps $0xb1, %%xmm3, %%xmm3 \n\t" \
      "shufps $0xb1, %%xmm4, %%xmm4 \n\t" \
      "shufps $0xb1, %%xmm5, %%xmm5 \n\t" \
      "mulps %2, %%xmm3 \n\t" \
      "mulps %2, %%xmm4 \n\t" \
      "mulps %2, %%xmm5 \n\t" \
      "addps %%xmm3, %%xmm0 \n\t" \
      "addps %%xmm4, %%xmm1 \n\t" \
      "addps %%xmm5, %%xmm2" \
      : \
      : \
      "m" (y), \
      "m" (x/y), \
      "m" (_sse_float_sgn13))
*/
// //////////////////////////////////////////////////////////////////////////  
// Adds xmm3,xmm4,xmm5 to xmm1,xmm2,xmm3
// (modified by Massimo Di Pierro to work with g++ instead of gcc)
// //////////////////////////////////////////////////////////////////////////  

#ifdef SSE2FIX
#define _sse_float_vector_add() \
_ASM ("addps %xmm3, %xmm0 \n\t" \
      "addps %xmm4, %xmm1 \n\t" \
      "addps %xmm5, %xmm2 \n\t" \
      : \
      : )
#else
#define _sse_float_vector_add() \
_ASM ("addps %%xmm3, %%xmm0 \n\t" \
      "addps %%xmm4, %%xmm1 \n\t" \
      "addps %%xmm5, %%xmm2 \n\t" \
      : \
      : )
#endif

// //////////////////////////////////////////////////////////////////////////  
// Subtracts xmm3,xmm4,xmm5 from xmm1,xmm2,xmm3
// (modified by Massimo Di Pierro to work with g++ instead of gcc)
// //////////////////////////////////////////////////////////////////////////  

#ifdef SSE2FIX
#define _sse_float_vector_sub() \
_ASM ("subps %xmm3, %xmm0 \n\t" \
      "subps %xmm4, %xmm1 \n\t" \
      "subps %xmm5, %xmm2" \
      : \
      :)
#else
#define _sse_float_vector_sub() \
_ASM ("subps %%xmm3, %%xmm0 \n\t" \
      "subps %%xmm4, %%xmm1 \n\t" \
      "subps %%xmm5, %%xmm2" \
      : \
      :)
#endif

// //////////////////////////////////////////////////////////////////////////  
// Multiplies the high words xmm3,xmm4,xmm5 with -1 and adds these registers// to xmm0,xmm1,xmm2
// //////////////////////////////////////////////////////////////////////////  

#define _sse_float_vector_addsub() \
_ASM ("mulps %0, %%xmm3 \n\t" \
      "mulps %0, %%xmm4 \n\t" \
      "mulps %0, %%xmm5 \n\t" \
      "addps %%xmm3, %%xmm0 \n\t" \
      "addps %%xmm4, %%xmm1 \n\t" \
      "addps %%xmm5, %%xmm2" \
      : \
      : \
      "m" (_sse_float_sgn34))

// //////////////////////////////////////////////////////////////////////////  
// Multiplies a pair sl,sh of su3 vectors with an su3 matrix u,
// assuming sl and sh are in the low and high words of xmm0,xmm1,xmm2
//
// On output the result is in xmm3,xmm4,xmm5 and the registers 
// xmm0,xmm1,xmm2 are changed
// //////////////////////////////////////////////////////////////////////////  

#define _sse_float_su3_multiply(u) { \
_ASM ("movss %0, %%xmm3 \n\t" \
      "movss %1, %%xmm6 \n\t" \
      "movss %2, %%xmm4 \n\t" \
      "movss %3, %%xmm7 \n\t" \
      "movss %4, %%xmm5 " \
      : \
      : \
      "m" ((u).c11.real()), \
      "m" ((u).c12.real()), \
      "m" ((u).c21.real()), \
      "m" ((u).c23.real()), \
      "m" ((u).c31.real())); \
_ASM ("shufps $0x0, %%xmm3, %%xmm3 \n\t" \
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
      "movss %0, %%xmm6 \n\t" \
      "movss %1, %%xmm7 \n\t" \
      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
      "mulps %%xmm1, %%xmm6 \n\t" \
      "mulps %%xmm2, %%xmm7 \n\t" \
      "addps %%xmm6, %%xmm5 \n\t" \
      "addps %%xmm7, %%xmm3 \n\t" \
      "movss %2, %%xmm6 \n\t" \
      "movss %3, %%xmm7 \n\t" \
      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
      "mulps %%xmm1, %%xmm6 \n\t" \
      "mulps %%xmm2, %%xmm7 \n\t" \
      "addps %%xmm6, %%xmm4 \n\t" \
      "addps %%xmm7, %%xmm5" \
      : \
      : \
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
      "mulps %4, %%xmm0 \n\t" \
      "mulps %4, %%xmm1 \n\t" \
      "mulps %4, %%xmm2 \n\t" \
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
      "addps %%xmm7, %%xmm4 " \
      : \
      : \
      "m" ((u).c11.imag()), \
      "m" ((u).c22.imag()), \
      "m" ((u).c33.imag()), \
      "m" ((u).c21.imag()), \
      "m" (_sse_float_sgn13)); \
_ASM ("movss %0, %%xmm6 \n\t" \
      "movss %1, %%xmm7 \n\t" \
      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
      "mulps %%xmm1, %%xmm6 \n\t" \
      "mulps %%xmm0, %%xmm7 \n\t" \
      "addps %%xmm6, %%xmm3 \n\t" \
      "addps %%xmm7, %%xmm5 \n\t" \
      "movss %2, %%xmm0 \n\t" \
      "movss %3, %%xmm6 \n\t" \
      "movss %4, %%xmm7 \n\t" \
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
      "m" ((u).c12.imag()), \
      "m" ((u).c31.imag()), \
      "m" ((u).c13.imag()), \
      "m" ((u).c32.imag()), \
      "m" ((u).c23.imag())); }

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
      "movss %4, %%xmm5 " \
      : \
      : \
      "m" ((u).c11.real()), \
      "m" ((u).c21.real()), \
      "m" ((u).c12.real()), \
      "m" ((u).c32.real()), \
      "m" ((u).c13.real())); \
_ASM ("shufps $0x0, %%xmm3, %%xmm3 \n\t" \
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
      "movss %0, %%xmm6 \n\t" \
      "movss %1, %%xmm7 \n\t" \
      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
      "mulps %%xmm1, %%xmm6 \n\t" \
      "mulps %%xmm2, %%xmm7 \n\t" \
      "addps %%xmm6, %%xmm5 \n\t" \
      "addps %%xmm7, %%xmm3 \n\t" \
      "movss %2, %%xmm6 \n\t" \
      "movss %3, %%xmm7 \n\t" \
      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
      "mulps %%xmm1, %%xmm6 \n\t" \
      "mulps %%xmm2, %%xmm7 \n\t" \
      "addps %%xmm6, %%xmm4 \n\t" \
      "addps %%xmm7, %%xmm5 " \
      : \
      : \
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
      "mulps %4, %%xmm0 \n\t" \
      "mulps %4, %%xmm1 \n\t" \
      "mulps %4, %%xmm2 \n\t" \
      "mulps %%xmm0, %%xmm6 \n\t" \
      "mulps %%xmm1, %%xmm7 \n\t" \
      "addps %%xmm6, %%xmm3 \n\t" \
      "addps %%xmm7, %%xmm4 \n\t" \
      "movss %2, %%xmm6 \n\t" \
      "movss %3, %%xmm7 " \
      : \
      : \
      "m" ((u).c11.imag()), \
      "m" ((u).c22.imag()), \
      "m" ((u).c33.imag()), \
      "m" ((u).c12.imag()), \
      "m" (_sse_float_sgn24)); \
_ASM ("shufps $0x0, %%xmm6, %%xmm6 \n\t" \
      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
      "mulps %%xmm2, %%xmm6 \n\t" \
      "mulps %%xmm0, %%xmm7 \n\t" \
      "addps %%xmm6, %%xmm5 \n\t" \
      "addps %%xmm7, %%xmm4 \n\t" \
      "movss %0, %%xmm6 \n\t" \
      "movss %1, %%xmm7 \n\t" \
      "shufps $0x0, %%xmm6, %%xmm6 \n\t" \
      "shufps $0x0, %%xmm7, %%xmm7 \n\t" \
      "mulps %%xmm1, %%xmm6 \n\t" \
      "mulps %%xmm0, %%xmm7 \n\t" \
      "addps %%xmm6, %%xmm3 \n\t" \
      "addps %%xmm7, %%xmm5 \n\t" \
      "movss %2, %%xmm0 \n\t" \
      "movss %3, %%xmm6 \n\t" \
      "movss %4, %%xmm7 \n\t" \
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
      "m" ((u).c21.imag()), \
      "m" ((u).c13.imag()), \
      "m" ((u).c31.imag()), \
      "m" ((u).c23.imag()), \
      "m" ((u).c32.imag())); }

// //////////////////////////////////////////////////////////////////////////  
//
// stuff used for optimized gamma matrix algebra (float)
//
// //////////////////////////////////////////////////////////////////////////  


// //////////////////////////////////////////////////////////////////////////  
// Multiplies the low words xmm3,xmm4,xmm5 with -1 and adds these registers
// to xmm0,xmm1,xmm2
// //////////////////////////////////////////////////////////////////////////  

#define _sse_float_vector_subadd() \
_ASM ("mulps %0, %%xmm3 \n\t" \
      "mulps %0, %%xmm4 \n\t" \
      "mulps %0, %%xmm5 \n\t" \
      "addps %%xmm3, %%xmm0 \n\t" \
      "addps %%xmm4, %%xmm1 \n\t" \
      "addps %%xmm5, %%xmm2" \
      : \
      : \
      "m" (_sse_float_sgn12))

// //////////////////////////////////////////////////////////////////////////  
// Multiplies xmm3,xmm4,xmm5 with i and adds them to xmm1,xmm2,xmm3
// //////////////////////////////////////////////////////////////////////////  

#define _sse_float_vector_i_add() \
_ASM ("shufps $0xb1, %%xmm3, %%xmm3 \n\t" \
      "shufps $0xb1, %%xmm4, %%xmm4 \n\t" \
      "shufps $0xb1, %%xmm5, %%xmm5 \n\t" \
      "mulps %0, %%xmm3 \n\t" \
      "mulps %0, %%xmm4 \n\t" \
      "mulps %0, %%xmm5 \n\t" \
      "addps %%xmm3, %%xmm0 \n\t" \
      "addps %%xmm4, %%xmm1 \n\t" \
      "addps %%xmm5, %%xmm2" \
      : \
      : \
      "m" (_sse_float_sgn13))

// //////////////////////////////////////////////////////////////////////////  
// Multiplies xmm3,xmm4,xmm5 with i and subtracts them from xmm1,xmm2,xmm3
// //////////////////////////////////////////////////////////////////////////  

#define _sse_float_vector_i_sub() \
_ASM ("shufps $0xb1, %%xmm3, %%xmm3 \n\t" \
      "shufps $0xb1, %%xmm4, %%xmm4 \n\t" \
      "shufps $0xb1, %%xmm5, %%xmm5 \n\t" \
      "mulps %0, %%xmm3 \n\t" \
      "mulps %0, %%xmm4 \n\t" \
      "mulps %0, %%xmm5 \n\t" \
      "addps %%xmm3, %%xmm0 \n\t" \
      "addps %%xmm4, %%xmm1 \n\t" \
      "addps %%xmm5, %%xmm2" \
      : \
      : \
      "m" (_sse_float_sgn24))

// //////////////////////////////////////////////////////////////////////////  
// Exchanges the high and low words of xmm3,xmm4,xmm5, multiplies them with i
// and adds the result to xmm1,xmm2,xmm3
// //////////////////////////////////////////////////////////////////////////  

#define _sse_float_vector_xch_i_add() \
_ASM ("shufps $0x1b, %%xmm3, %%xmm3 \n\t" \
      "shufps $0x1b, %%xmm4, %%xmm4 \n\t" \
      "shufps $0x1b, %%xmm5, %%xmm5 \n\t" \
      "mulps %0, %%xmm3 \n\t" \
      "mulps %0, %%xmm4 \n\t" \
      "mulps %0, %%xmm5 \n\t" \
      "addps %%xmm3, %%xmm0 \n\t" \
      "addps %%xmm4, %%xmm1 \n\t" \
      "addps %%xmm5, %%xmm2" \
      : \
      : \
      "m" (_sse_float_sgn13))

// //////////////////////////////////////////////////////////////////////////  
// Exchanges the high and low words of xmm3,xmm4,xmm5, multiplies them with i
// and subtracts the result from xmm1,xmm2,xmm3
// //////////////////////////////////////////////////////////////////////////  

#define _sse_float_vector_xch_i_sub() \
_ASM ("shufps $0x1b, %%xmm3, %%xmm3 \n\t" \
      "shufps $0x1b, %%xmm4, %%xmm4 \n\t" \
      "shufps $0x1b, %%xmm5, %%xmm5 \n\t" \
      "mulps %0, %%xmm3 \n\t" \
      "mulps %0, %%xmm4 \n\t" \
      "mulps %0, %%xmm5 \n\t" \
      "addps %%xmm3, %%xmm0 \n\t" \
      "addps %%xmm4, %%xmm1 \n\t" \
      "addps %%xmm5, %%xmm2" \
      : \
      : \
      "m" (_sse_float_sgn24))

// //////////////////////////////////////////////////////////////////////////  
// Multiplies the low and high words of xmm3,xmm4,xmm5 with i and -i
// respectively and adds these registers to xmm1,xmm2,xmm3
// //////////////////////////////////////////////////////////////////////////  

#define _sse_float_vector_i_addsub() \
_ASM ("shufps $0xb1, %%xmm3, %%xmm3 \n\t" \
      "shufps $0xb1, %%xmm4, %%xmm4 \n\t" \
      "shufps $0xb1, %%xmm5, %%xmm5 \n\t" \
      "mulps %0, %%xmm3 \n\t" \
      "mulps %0, %%xmm4 \n\t" \
      "mulps %0, %%xmm5 \n\t" \
      "addps %%xmm3, %%xmm0 \n\t" \
      "addps %%xmm4, %%xmm1 \n\t" \
      "addps %%xmm5, %%xmm2" \
      : \
      : \
      "m" (_sse_float_sgn14))

// //////////////////////////////////////////////////////////////////////////  
// Multiplies the low and high words of xmm3,xmm4,xmm5 with -i and i
// respectively and adds these registers to xmm1,xmm2,xmm3
// //////////////////////////////////////////////////////////////////////////  

#define _sse_float_vector_i_subadd() \
_ASM ("shufps $0xb1, %%xmm3, %%xmm3 \n\t" \
      "shufps $0xb1, %%xmm4, %%xmm4 \n\t" \
      "shufps $0xb1, %%xmm5, %%xmm5 \n\t" \
      "mulps %0, %%xmm3 \n\t" \
      "mulps %0, %%xmm4 \n\t" \
      "mulps %0, %%xmm5 \n\t" \
      "addps %%xmm3, %%xmm0 \n\t" \
      "addps %%xmm4, %%xmm1 \n\t" \
      "addps %%xmm5, %%xmm2" \
      : \
      : \
      "m" (_sse_float_sgn23))

// //////////////////////////////////////////////////////////////////////////  
// Exchanges the high and low words in xmm3,xmm4,xmm5 
// //////////////////////////////////////////////////////////////////////////  

#ifdef SSE2FIX 
#define _sse_float_vector_xch() \
_ASM ("shufps $0x4e, %xmm3, %xmm3 \n\t" \
      "shufps $0x4e, %xmm4, %xmm4 \n\t" \
      "shufps $0x4e, %xmm5, %xmm5" \
      : \
      :)
#else
#define _sse_float_vector_xch() \
_ASM ("shufps $0x4e, %%xmm3, %%xmm3 \n\t" \
      "shufps $0x4e, %%xmm4, %%xmm4 \n\t" \
      "shufps $0x4e, %%xmm5, %%xmm5" \
      : \
      :)
#endif



// //////////////////////////////////////////////////////////////////////////
// Cache manipulation macros (double)
// //////////////////////////////////////////////////////////////////////////

#define _sse_double_prefetch_16(addr) \
_ASM ("prefetcht0 %0" \
      : \
      : "m" (*(addr)))

#define _sse_double_prefetch_spinor(addr) \
_ASM ("prefetcht0 %0 \n\t" \
      "prefetcht0 %1" \
      : \
      : \
      "m" (*(((char*)(((unsigned int)(addr))&~0x7f)))), \
      "m" (*(((char*)(((unsigned int)(addr))&~0x7f))+128)))

#define _sse_double_prefetch_nta_spinor(addr) \
_ASM ("prefetchnta %0 \n\t" \
      "prefetchnta %1" \
      : \
      : \
      "m" (*(((char*)(((unsigned int)(addr))&~0x7f)))), \
      "m" (*(((char*)(((unsigned int)(addr))&~0x7f))+128)))

#define _sse_double_prefetch_su3(addr) \
_ASM ("prefetcht0 %0 \n\t" \
      "prefetcht0 %1" \
      : \
      : \
      "m" (*(((char*)(((unsigned int)(addr))&~0x7f)))), \
      "m" (*(((char*)(((unsigned int)(addr))&~0x7f))+128)))


// //////////////////////////////////////////////////////////////////////////
//
// Operations for SU3 color linear algebra used in mdp_matrix (double)
//
// //////////////////////////////////////////////////////////////////////////

// //////////////////////////////////////////////////////////////////////////
// Loads an su3 vector s to xmm0,xmm1,xmm2
// //////////////////////////////////////////////////////////////////////////

#define _sse_double_load(s) \
_ASM ("movapd %0, %%xmm0 \n\t" \
      "movapd %1, %%xmm1 \n\t" \
      "movapd %2, %%xmm2" \
      : \
      : \
      "m" ((s).c1), \
      "m" ((s).c2), \
      "m" ((s).c3))

#define _sse_double_load_123(c1, c2, c3) \
_ASM ("movapd %0, %%xmm0 \n\t" \
      "movapd %1, %%xmm1 \n\t" \
      "movapd %2, %%xmm2" \
      : \
      : \
      "m" (c1), \
      "m" (c2), \
      "m" (c3))


// //////////////////////////////////////////////////////////////////////////
// Loads an su3 vector s to xmm3,xmm4,xmm5
// //////////////////////////////////////////////////////////////////////////  

#define _sse_double_load_up(s) \
_ASM ("movapd %0, %%xmm3 \n\t" \
      "movapd %1, %%xmm4 \n\t" \
      "movapd %2, %%xmm5" \
      : \
      : \
      "m" ((s).c1), \
      "m" ((s).c2), \
      "m" ((s).c3))

#define _sse_double_load_up_123(c1, c2, c3) \
_ASM ("movapd %0, %%xmm3 \n\t" \
      "movapd %1, %%xmm4 \n\t" \
      "movapd %2, %%xmm5" \
      : \
      : \
      "m" (c1), \
      "m" (c2), \
      "m" (c3))

// //////////////////////////////////////////////////////////////////////////
// Stores xmm0,xmm1,xmm2 to the components r.c1,r.c2,r.c3 of an su3 vector
// //////////////////////////////////////////////////////////////////////////

#define _sse_double_store(r) \
_ASM ("movapd %%xmm0, %0 \n\t" \
      "movapd %%xmm1, %1 \n\t" \
      "movapd %%xmm2, %2" \
      : \
      "=m" ((r).c1), \
      "=m" ((r).c2), \
      "=m" ((r).c3))

#define _sse_double_store_123(c1, c2, c3) \
_ASM ("movapd %%xmm0, %0 \n\t" \
      "movapd %%xmm1, %1 \n\t" \
      "movapd %%xmm2, %2" \
      : \
      "=m" (c1), \
      "=m" (c2), \
      "=m" (c3))

// //////////////////////////////////////////////////////////////////////////
// Stores xmm3,xmm4,xmm5 to the components r.c1,r.c2,r.c3 of an su3 vector
// //////////////////////////////////////////////////////////////////////////

#define _sse_double_store_up(r) \
_ASM ("movapd %%xmm3, %0 \n\t" \
      "movapd %%xmm4, %1 \n\t" \
      "movapd %%xmm5, %2" \
      : \
      "=m" ((r).c1), \
      "=m" ((r).c2), \
      "=m" ((r).c3))

#define _sse_double_store_up_123(c1, c2, c3) \
_ASM ("movapd %%xmm3, %0 \n\t" \
      "movapd %%xmm4, %1 \n\t" \
      "movapd %%xmm5, %2" \
      : \
      "=m" (c1), \
      "=m" (c2), \
      "=m" (c3))

// //////////////////////////////////////////////////////////////////////////
// Multiplies xmm0,xmm1,xmm2 with a constant _sse_double c
// //////////////////////////////////////////////////////////////////////////

#define _sse_double_vector_mul(c) \
_ASM ("mulpd %0, %%xmm0 \n\t" \
      "mulpd %0, %%xmm1 \n\t" \
      "mulpd %0, %%xmm2" \
      : \
      : \
      "m" (c))

// //////////////////////////////////////////////////////////////////////////
// multiplies xmm0, xmm1, xmm2 for complex=a*(b+I)
// //////////////////////////////////////////////////////////////////////////
/* deprecated
#define _sse_double_vector_mulc(a,b) \
_ASM ("mulpd %0, %%xmm0 \n\t" \
      "mulpd %0, %%xmm1 \n\t" \
      "mulpd %0, %%xmm2 \n\t" \
      "movapd %%xmm0, %%xmm3 \n\t" \
      "movapd %%xmm1, %%xmm4 \n\t" \
      "movapd %%xmm2, %%xmm5 \n\t" \
      "mulpd %1, %%xmm0 \n\t" \
      "mulpd %1, %%xmm1 \n\t" \
      "mulpd %1, %%xmm2 \n\t" \
      "shufpd $0x1, %%xmm3, %%xmm3 \n\t" \
      "shufpd $0x1, %%xmm4, %%xmm4 \n\t" \
      "shufpd $0x1, %%xmm5, %%xmm5 \n\t" \
      "xorpd %2, %%xmm3 \n\t" \
      "xorpd %2, %%xmm4 \n\t" \
      "xorpd %2, %%xmm5 \n\t" \
      "addpd %%xmm3, %%xmm0 \n\t" \
      "addpd %%xmm4, %%xmm1 \n\t" \
      "addpd %%xmm5, %%xmm2" \
      : \
      : \
      "m" (a), \
      "m" (b), \
      "m" (_sse_double_sgn))      
*/
// //////////////////////////////////////////////////////////////////////////
// multiplies xmm0, xmm1, xmm2 for complex=x+I*y
// //////////////////////////////////////////////////////////////////////////

#define _sse_double_vector_mul_complex(x,y) \
_ASM ("movapd %%xmm0, %%xmm3 \n\t" \
      "movapd %%xmm1, %%xmm4 \n\t" \
      "movapd %%xmm2, %%xmm5 \n\t" \
      "mulpd %1, %%xmm3 \n\t" \
      "mulpd %1, %%xmm4 \n\t" \
      "mulpd %1, %%xmm5 \n\t" \
      "shufpd $0x1, %%xmm3, %%xmm3 \n\t" \
      "shufpd $0x1, %%xmm4, %%xmm4 \n\t" \
      "shufpd $0x1, %%xmm5, %%xmm5 \n\t" \
      "xorpd %2, %%xmm3 \n\t" \
      "xorpd %2, %%xmm4 \n\t" \
      "xorpd %2, %%xmm5 \n\t" \
      "mulpd %0, %%xmm0 \n\t" \
      "mulpd %0, %%xmm1 \n\t" \
      "mulpd %0, %%xmm2 \n\t" \
      "addpd %%xmm0, %%xmm3 \n\t" \
      "addpd %%xmm1, %%xmm4 \n\t" \
      "addpd %%xmm2, %%xmm5" \
      : \
      : \
      "m" (x), \
      "m" (y), \
      "m" (_sse_double_sgn))      

// //////////////////////////////////////////////////////////////////////////
// Adds xmm3,xmm4,xmm5 to xmm1,xmm2,xmm3
// //////////////////////////////////////////////////////////////////////////

#ifdef SSE2FIX
#define _sse_double_vector_add() \
_ASM ("addpd %xmm3, %xmm0 \n\t" \
      "addpd %xmm4, %xmm1 \n\t" \
      "addpd %xmm5, %xmm2" \
      : \
      :)
#else
#define _sse_double_vector_add() \
_ASM ("addpd %%xmm3, %%xmm0 \n\t" \
      "addpd %%xmm4, %%xmm1 \n\t" \
      "addpd %%xmm5, %%xmm2" \
      : \
      :)
#endif

// //////////////////////////////////////////////////////////////////////////  
// Subtracts xmm3,xmm4,xmm5 from xmm1,xmm2,xmm3
// //////////////////////////////////////////////////////////////////////////

#ifdef SSE2FIX
#define _sse_double_vector_sub() \
_ASM ("subpd %xmm3, %xmm0 \n\t" \
      "subpd %xmm4, %xmm1 \n\t" \
      "subpd %xmm5, %xmm2" \
      : \
      :)
#else
#define _sse_double_vector_sub() \
_ASM ("subpd %%xmm3, %%xmm0 \n\t" \
      "subpd %%xmm4, %%xmm1 \n\t" \
      "subpd %%xmm5, %%xmm2" \
      : \
      :)
#endif

// //////////////////////////////////////////////////////////////////////////
// Multiplies an su3 vector s with an su3 matrix u, assuming s is
// stored in  xmm0,xmm1,xmm2
// On output the result is in xmm3,xmm4,xmm5 and the registers 
// xmm0,xmm1,xmm2 are changed
// //////////////////////////////////////////////////////////////////////////

#define _sse_double_su3_multiply(u) { \
_ASM ("movsd %0, %%xmm3 \n\t" \
      "movsd %1, %%xmm6 \n\t" \
      "movsd %2, %%xmm4 \n\t" \
      "movsd %3, %%xmm7 \n\t" \
      "movsd %4, %%xmm5 " \
      : \
      : \
      "m" ((u).c11.real()), \
      "m" ((u).c12.real()), \
      "m" ((u).c21.real()), \
      "m" ((u).c23.real()), \
      "m" ((u).c31.real())); \
_ASM ("unpcklpd %%xmm3, %%xmm3 \n\t" \
      "unpcklpd %%xmm6, %%xmm6 \n\t" \
      "unpcklpd %%xmm4, %%xmm4 \n\t" \
      "mulpd %%xmm0, %%xmm3 \n\t" \
      "unpcklpd %%xmm7, %%xmm7 \n\t" \
      "mulpd %%xmm1, %%xmm6 \n\t" \
      "unpcklpd %%xmm5, %%xmm5 \n\t" \
      "mulpd %%xmm0, %%xmm4 \n\t" \
      "addpd %%xmm6, %%xmm3 \n\t" \
      "mulpd %%xmm2, %%xmm7 \n\t" \
      "mulpd %%xmm0, %%xmm5 \n\t" \
      "addpd %%xmm7, %%xmm4 \n\t" \
      "movsd %0, %%xmm6 \n\t" \
      "movsd %1, %%xmm7 \n\t" \
      "unpcklpd %%xmm6, %%xmm6 \n\t" \
      "unpcklpd %%xmm7, %%xmm7 \n\t" \
      "mulpd %%xmm1, %%xmm6 \n\t" \
      "mulpd %%xmm2, %%xmm7 \n\t" \
      "addpd %%xmm6, %%xmm5 \n\t" \
      "addpd %%xmm7, %%xmm3 \n\t" \
      "movsd %2, %%xmm6 \n\t" \
      "movsd %3, %%xmm7 \n\t" \
      "unpcklpd %%xmm6, %%xmm6 \n\t" \
      "unpcklpd %%xmm7, %%xmm7 \n\t" \
      "mulpd %%xmm1, %%xmm6 \n\t" \
      "mulpd %%xmm2, %%xmm7 \n\t" \
      "addpd %%xmm6, %%xmm4 \n\t" \
      "addpd %%xmm7, %%xmm5 " \
      : \
      : \
      "m" ((u).c32.real()), \
      "m" ((u).c13.real()), \
      "m" ((u).c22.real()), \
      "m" ((u).c33.real())); \
_ASM ("movsd %0, %%xmm6 \n\t" \
      "movsd %1, %%xmm7 \n\t" \
      "shufpd $0x1, %%xmm0, %%xmm0 \n\t" \
      "shufpd $0x1, %%xmm1, %%xmm1 \n\t" \
      "shufpd $0x1, %%xmm2, %%xmm2 \n\t" \
      "unpcklpd %%xmm6, %%xmm6 \n\t" \
      "unpcklpd %%xmm7, %%xmm7 \n\t" \
      "xorpd %4, %%xmm0 \n\t" \
      "xorpd %4, %%xmm1 \n\t" \
      "xorpd %4, %%xmm2 \n\t" \
      "mulpd %%xmm0, %%xmm6 \n\t" \
      "mulpd %%xmm1, %%xmm7 \n\t" \
      "addpd %%xmm6, %%xmm3 \n\t" \
      "addpd %%xmm7, %%xmm4 \n\t" \
      "movsd %2, %%xmm6 \n\t" \
      "movsd %3, %%xmm7 " \
      : \
      : \
      "m" ((u).c11.imag()), \
      "m" ((u).c22.imag()), \
      "m" ((u).c33.imag()), \
      "m" ((u).c21.imag()), \
      "m" (_sse_double_sgn)); \
_ASM ("unpcklpd %%xmm6, %%xmm6 \n\t" \
      "unpcklpd %%xmm7, %%xmm7 \n\t" \
      "mulpd %%xmm2, %%xmm6 \n\t" \
      "mulpd %%xmm0, %%xmm7 \n\t" \
      "addpd %%xmm6, %%xmm5 \n\t" \
      "addpd %%xmm7, %%xmm4 \n\t" \
      "movsd %0, %%xmm6 \n\t" \
      "movsd %1, %%xmm7 \n\t" \
      "unpcklpd %%xmm6, %%xmm6 \n\t" \
      "unpcklpd %%xmm7, %%xmm7 \n\t" \
      "mulpd %%xmm1, %%xmm6 \n\t" \
      "mulpd %%xmm0, %%xmm7 \n\t" \
      "addpd %%xmm6, %%xmm3 \n\t" \
      "addpd %%xmm7, %%xmm5 \n\t" \
      "movsd %2, %%xmm0 \n\t" \
      "movsd %3, %%xmm6 \n\t" \
      "movsd %4, %%xmm7 \n\t" \
      "unpcklpd %%xmm0, %%xmm0 \n\t" \
      "unpcklpd %%xmm6, %%xmm6 \n\t" \
      "unpcklpd %%xmm7, %%xmm7 \n\t" \
      "mulpd %%xmm2, %%xmm0 \n\t" \
      "mulpd %%xmm1, %%xmm6 \n\t" \
      "mulpd %%xmm2, %%xmm7 \n\t" \
      "addpd %%xmm0, %%xmm3 \n\t" \
      "addpd %%xmm6, %%xmm5 \n\t" \
      "addpd %%xmm7, %%xmm4 " \
      : \
      : \
      "m" ((u).c12.imag()), \
      "m" ((u).c31.imag()), \
      "m" ((u).c13.imag()), \
      "m" ((u).c32.imag()), \
      "m" ((u).c23.imag())); }

// //////////////////////////////////////////////////////////////////////////
// Multiplies an su3 vector s with an su3 matrix u^dagger, assuming s is
// stored in  xmm0,xmm1,xmm2
//
// On output the result is in xmm3,xmm4,xmm5 and the registers 
// xmm0,xmm1,xmm2 are changed
// //////////////////////////////////////////////////////////////////////////

#define _sse_double_su3_inverse_multiply(u) { \
_ASM ("movsd %0, %%xmm3 \n\t" \
      "movsd %1, %%xmm6 \n\t" \
      "movsd %2, %%xmm4 \n\t" \
      "movsd %3, %%xmm7 \n\t" \
      "movsd %4, %%xmm5 " \
      : \
      : \
      "m" ((u).c11.real()), \
      "m" ((u).c21.real()), \
      "m" ((u).c12.real()), \
      "m" ((u).c32.real()), \
      "m" ((u).c13.real())); \
_ASM ("unpcklpd %%xmm3, %%xmm3 \n\t" \
      "unpcklpd %%xmm6, %%xmm6 \n\t" \
      "unpcklpd %%xmm4, %%xmm4 \n\t" \
      "mulpd %%xmm0, %%xmm3 \n\t" \
      "unpcklpd %%xmm7, %%xmm7 \n\t" \
      "mulpd %%xmm1, %%xmm6 \n\t" \
      "unpcklpd %%xmm5, %%xmm5 \n\t" \
      "mulpd %%xmm0, %%xmm4 \n\t" \
      "addpd %%xmm6, %%xmm3 \n\t" \
      "mulpd %%xmm2, %%xmm7 \n\t" \
      "mulpd %%xmm0, %%xmm5 \n\t" \
      "addpd %%xmm7, %%xmm4 \n\t" \
      "movsd %0, %%xmm6 \n\t" \
      "movsd %1, %%xmm7 \n\t" \
      "unpcklpd %%xmm6, %%xmm6 \n\t" \
      "unpcklpd %%xmm7, %%xmm7 \n\t" \
      "mulpd %%xmm1, %%xmm6 \n\t" \
      "mulpd %%xmm2, %%xmm7 \n\t" \
      "addpd %%xmm6, %%xmm5 \n\t" \
      "addpd %%xmm7, %%xmm3 \n\t" \
      "movsd %2, %%xmm6 \n\t" \
      "movsd %3, %%xmm7 \n\t" \
      "unpcklpd %%xmm6, %%xmm6 \n\t" \
      "unpcklpd %%xmm7, %%xmm7 \n\t" \
      "mulpd %%xmm1, %%xmm6 \n\t" \
      "mulpd %%xmm2, %%xmm7 \n\t" \
      "addpd %%xmm6, %%xmm4 \n\t" \
      "addpd %%xmm7, %%xmm5" \
      : \
      : \
      "m" ((u).c23.real()), \
      "m" ((u).c31.real()), \
      "m" ((u).c22.real()), \
      "m" ((u).c33.real())); \
_ASM ("movsd %0, %%xmm6 \n\t" \
      "movsd %1, %%xmm7 \n\t" \
      "xorpd %4, %%xmm0 \n\t" \
      "xorpd %4, %%xmm1 \n\t" \
      "xorpd %4, %%xmm2 \n\t" \
      "unpcklpd %%xmm6, %%xmm6 \n\t" \
      "unpcklpd %%xmm7, %%xmm7 \n\t" \
      "shufpd $0x1, %%xmm0, %%xmm0 \n\t" \
      "shufpd $0x1, %%xmm1, %%xmm1 \n\t" \
      "shufpd $0x1, %%xmm2, %%xmm2 \n\t" \
      "mulpd %%xmm0, %%xmm6 \n\t" \
      "mulpd %%xmm1, %%xmm7 \n\t" \
      "addpd %%xmm6, %%xmm3 \n\t" \
      "addpd %%xmm7, %%xmm4 \n\t" \
      "movsd %2, %%xmm6 \n\t" \
      "movsd %3, %%xmm7 " \
      : \
      : \
      "m" ((u).c11.imag()), \
      "m" ((u).c22.imag()), \
      "m" ((u).c33.imag()), \
      "m" ((u).c12.imag()), \
      "m" (_sse_double_sgn)); \
_ASM ("unpcklpd %%xmm6, %%xmm6 \n\t" \
      "unpcklpd %%xmm7, %%xmm7 \n\t" \
      "mulpd %%xmm2, %%xmm6 \n\t" \
      "mulpd %%xmm0, %%xmm7 \n\t" \
      "addpd %%xmm6, %%xmm5 \n\t" \
      "addpd %%xmm7, %%xmm4 \n\t" \
      "movsd %0, %%xmm6 \n\t" \
      "movsd %1, %%xmm7 \n\t" \
      "unpcklpd %%xmm6, %%xmm6 \n\t" \
      "unpcklpd %%xmm7, %%xmm7 \n\t" \
      "mulpd %%xmm1, %%xmm6 \n\t" \
      "mulpd %%xmm0, %%xmm7 \n\t" \
      "addpd %%xmm6, %%xmm3 \n\t" \
      "addpd %%xmm7, %%xmm5 \n\t" \
      "movsd %2, %%xmm0 \n\t" \
      "movsd %3, %%xmm6 \n\t" \
      "movsd %4, %%xmm7 \n\t" \
      "unpcklpd %%xmm0, %%xmm0 \n\t" \
      "unpcklpd %%xmm6, %%xmm6 \n\t" \
      "unpcklpd %%xmm7, %%xmm7 \n\t" \
      "mulpd %%xmm2, %%xmm0 \n\t" \
      "mulpd %%xmm1, %%xmm6 \n\t" \
      "mulpd %%xmm2, %%xmm7 \n\t" \
      "addpd %%xmm0, %%xmm3 \n\t" \
      "addpd %%xmm6, %%xmm5 \n\t" \
      "addpd %%xmm7, %%xmm4 " \
      : \
      : \
      "m" ((u).c21.imag()), \
      "m" ((u).c13.imag()), \
      "m" ((u).c31.imag()), \
      "m" ((u).c23.imag()), \
      "m" ((u).c32.imag())); }

// //////////////////////////////////////////////////////////////////////////  
//
// stuff used for optimized gamma matrix algebra (float)
//
// //////////////////////////////////////////////////////////////////////////  

// //////////////////////////////////////////////////////////////////////////
// Multiplies xmm3,xmm4,xmm5 with i
// //////////////////////////////////////////////////////////////////////////

#define _sse_double_vector_i_mul() \
_ASM ("shufpd $0x1, %%xmm3, %%xmm3 \n\t" \
      "shufpd $0x1, %%xmm4, %%xmm4 \n\t" \
      "shufpd $0x1, %%xmm5, %%xmm5 \n\t" \
      "xorpd %0, %%xmm3 \n\t" \
      "xorpd %0, %%xmm4 \n\t" \
      "xorpd %0, %%xmm5" \
      : \
      : \
      "m" (_sse_double_sgn))

// //////////////////////////////////////////////////////////////////////////
// Multiplies xmm3,xmm4,xmm5 with -i
// //////////////////////////////////////////////////////////////////////////

#define _sse_double_vector_minus_i_mul() \
_ASM ("xorpd %0, %%xmm3 \n\t" \
      "xorpd %0, %%xmm4 \n\t" \
      "xorpd %0, %%xmm5 \n\t" \
      "shufpd $0x1, %%xmm3, %%xmm3 \n\t" \
      "shufpd $0x1, %%xmm4, %%xmm4 \n\t" \
      "shufpd $0x1, %%xmm5, %%xmm5" \
      : \
      : \
      "m" (_sse_double_sgn))


// //////////////////////////////////////////////////////////////////////////  
//
// stuff used for optimized vector operations
//
// //////////////////////////////////////////////////////////////////////////  

// //////////////////////////////////////////////////////////////////////////
// r+=s[0].c1+s[0].c2+...+s[7].c1+s[7].c2;
// //////////////////////////////////////////////////////////////////////////

#define _sse_double_add_norm_square_16(r,c) { \
_ASM ("movapd %0, %%xmm0 \n\t" \
      "movapd %1, %%xmm1 \n\t" \
      "movapd %2, %%xmm2 \n\t" \
      "movapd %3, %%xmm3" \
      : \
      : \
      "m" (*((r))), \
      "m" (*((r)+1)), \
      "m" (*((r)+2)), \
      "m" (*((r)+3))); \
_ASM ("movapd %0, %%xmm4 \n\t" \
      "movapd %1, %%xmm5 \n\t" \
      "movapd %2, %%xmm6 \n\t" \
      "movapd %3, %%xmm7 \n\t" \
      "mulpd %%xmm0, %%xmm0 \n\t" \
      "mulpd %%xmm1, %%xmm1 \n\t" \
      "mulpd %%xmm2, %%xmm2 \n\t" \
      "mulpd %%xmm3, %%xmm3 \n\t" \
      "mulpd %%xmm4, %%xmm4 \n\t" \
      "mulpd %%xmm5, %%xmm5 \n\t" \
      "mulpd %%xmm6, %%xmm6 \n\t" \
      "mulpd %%xmm7, %%xmm7 \n\t" \
      "addpd %%xmm0, %%xmm1 \n\t" \
      "addpd %%xmm2, %%xmm3 \n\t" \
      "addpd %%xmm4, %%xmm5 \n\t" \
      "addpd %%xmm6, %%xmm7 \n\t" \
      "addpd %%xmm1, %%xmm3 \n\t" \
      "addpd %%xmm5, %%xmm7 \n\t" \
      "addpd %%xmm3, %%xmm7" \
      : \
      : \
      "m" (*((r)+4)), \
      "m" (*((r)+5)), \
      "m" (*((r)+6)), \
      "m" (*((r)+7))); \
_ASM ("movapd %0,     %%xmm1 \n\t" \
      "addpd  %%xmm1, %%xmm7 \n\t" \
      "movapd %%xmm7, %0" \
      : \
      "=m" (c)); }

// //////////////////////////////////////////////////////////////////////////
// r[0].c1=s[0].c1, r[0].2=s[0].c2+...+r[7].c2=s[7].c2;
// //////////////////////////////////////////////////////////////////////////

#define _sse_double_add_real_scalar_product_16(r,s,c) { \
_ASM ("movapd %0, %%xmm0 \n\t" \
      "movapd %1, %%xmm1 \n\t" \
      "movapd %2, %%xmm2 \n\t" \
      "movapd %3, %%xmm3 \n\t" \
      : \
      : \
      "m" (*((r))), \
      "m" (*((r)+1)), \
      "m" (*((r)+2)), \
      "m" (*((r)+3))); \
_ASM ("mulpd %0, %%xmm0 \n\t"		\
      "mulpd %1, %%xmm1 \n\t" \
      "mulpd %2, %%xmm2 \n\t" \
      "mulpd %3, %%xmm3 \n\t" \
      : \
      : \
      "m" (*((s))), \
      "m" (*((s)+1)), \
      "m" (*((s)+2)), \
      "m" (*((s)+3))); \
_ASM ("movapd %0, %%xmm4 \n\t" \
      "movapd %1, %%xmm5 \n\t" \
      "movapd %2, %%xmm6 \n\t" \
      "movapd %3, %%xmm7 \n\t" \
      : \
      : \
      "m" (*((r)+4)), \
      "m" (*((r)+5)), \
      "m" (*((r)+6)), \
      "m" (*((r)+7)));	      \
_ASM ("mulpd %0, %%xmm4 \n\t"		\
      "mulpd %1, %%xmm5 \n\t" \
      "mulpd %2, %%xmm6 \n\t" \
      "mulpd %3, %%xmm7 \n\t" \
      : \
      : \
      "m" (*((s)+4)), \
      "m" (*((s)+5)), \
      "m" (*((s)+6)), \
      "m" (*((s)+7))); \
_ASM ("addpd %%xmm0, %%xmm1 \n\t" \
      "addpd %%xmm2, %%xmm3 \n\t" \
      "addpd %%xmm4, %%xmm5 \n\t" \
      "addpd %%xmm6, %%xmm7 \n\t" \
      "addpd %%xmm1, %%xmm3 \n\t" \
      "addpd %%xmm5, %%xmm7 \n\t" \
      "addpd %%xmm3, %%xmm7 \n\t" \
      "movapd %0, %%xmm1    \n\t" \
      "addpd %%xmm1, %%xmm7 \n\t" \
      "movapd %%xmm7, %0    \n\t" \
      : \
      "=m" (c)); }

#define _sse_double_add_imag_scalar_product_16(r,s,c) { \
_ASM ("movapd %0, %%xmm0 \n\t" \
      "movapd %1, %%xmm1 \n\t" \
      "movapd %2, %%xmm2 \n\t" \
      "movapd %3, %%xmm3 \n\t" \
      : \
      : \
      "m" (*((r))), \
      "m" (*((r)+1)), \
      "m" (*((r)+2)), \
      "m" (*((r)+3))); \
_ASM ("shufpd $0x1, %%xmm0, %%xmm0 \n\t"	\
      "shufpd $0x1, %%xmm1, %%xmm1 \n\t" \
      "shufpd $0x1, %%xmm2, %%xmm2 \n\t" \
      "shufpd $0x1, %%xmm3, %%xmm3 \n\t" \
      "mulpd %0, %%xmm0 \n\t" \
      "mulpd %1, %%xmm1 \n\t" \
      "mulpd %2, %%xmm2 \n\t" \
      "mulpd %3, %%xmm3 \n\t" \
      : \
      : \
      "m" (*((s))), \
      "m" (*((s)+1)), \
      "m" (*((s)+2)), \
      "m" (*((s)+3))); \
_ASM ("movapd %0, %%xmm4 \n\t" \
      "movapd %1, %%xmm5 \n\t" \
      "movapd %2, %%xmm6 \n\t" \
      "movapd %3, %%xmm7 \n\t" \
      : \
      : \
      "m" (*((r)+4)), \
      "m" (*((r)+5)), \
      "m" (*((r)+6)), \
      "m" (*((r)+7)));			 \
_ASM ("shufpd $0x1, %%xmm4, %%xmm4 \n\t"	\
      "shufpd $0x1, %%xmm5, %%xmm5 \n\t" \
      "shufpd $0x1, %%xmm6, %%xmm6 \n\t" \
      "shufpd $0x1, %%xmm7, %%xmm7 \n\t" \
      "mulpd %0, %%xmm4 \n\t" \
      "mulpd %1, %%xmm5 \n\t" \
      "mulpd %2, %%xmm6 \n\t" \
      "mulpd %3, %%xmm7 \n\t" \
      : \
      : \
      "m" (*((s)+4)), \
      "m" (*((s)+5)), \
      "m" (*((s)+6)), \
      "m" (*((s)+7))); \
_ASM ("addpd %%xmm0, %%xmm1 \n\t" \
      "addpd %%xmm2, %%xmm3 \n\t" \
      "addpd %%xmm4, %%xmm5 \n\t" \
      "addpd %%xmm6, %%xmm7 \n\t" \
      "addpd %%xmm1, %%xmm3 \n\t" \
      "addpd %%xmm5, %%xmm7 \n\t" \
      "addpd %%xmm3, %%xmm7 \n\t" \
      "movapd %0, %%xmm1    \n\t" \
      "addpd %%xmm1, %%xmm7 \n\t" \
      "movapd %%xmm7, %0    \n\t" \
      : \
      "=m" (c)); }

#define _sse_double_hermitian_su3(r,s) { \
_ASM ("movapd %0, %%xmm0 \n\t"\
      "xorpd  %3, %%xmm0 \n\t" \
      "movapd %1, %%xmm1 \n\t"\
      "xorpd  %3, %%xmm1 \n\t" \
      "movapd %2, %%xmm2 \n\t"\
      "xorpd  %3, %%xmm2 \n\t" \
      : \
      : \
      "m" (*((s))), \
      "m" (*((s)+4)), \
      "m" (*((s)+8)), \
      "m" (_sse_double_sgn2)); \
_ASM ("movapd %%xmm0, %0 \n\t"\
      "movapd %%xmm1, %1 \n\t"\
      "movapd %%xmm2, %2 \n\t"\
      : \
      "=m" (*((r))), \
      "=m" (*((r)+4)), \
      "=m" (*((r)+8))); \
_ASM ("movapd %0, %%xmm0 \n\t"\
      "xorpd  %3, %%xmm0 \n\t" \
      "movapd %1, %%xmm1 \n\t"\
      "xorpd  %3, %%xmm1 \n\t" \
      "movapd %2, %%xmm2 \n\t"\
      "xorpd  %3, %%xmm2 \n\t" \
      : \
      : \
      "m" (*((s)+1)), \
      "m" (*((s)+2)), \
      "m" (*((s)+5)), \
      "m" (_sse_double_sgn2)); \
_ASM ("movapd %%xmm0, %0 \n\t"\
      "movapd %%xmm1, %1 \n\t"\
      "movapd %%xmm2, %2 \n\t"\
      : \
      "=m" (*((r)+3)), \
      "=m" (*((r)+6)), \
      "=m" (*((r)+7))); \
_ASM ("movapd %0, %%xmm0 \n\t"\
      "xorpd  %3, %%xmm0 \n\t" \
      "movapd %1, %%xmm1 \n\t"\
      "xorpd  %3, %%xmm1 \n\t" \
      "movapd %2, %%xmm2 \n\t"\
      "xorpd  %3, %%xmm2 \n\t" \
      : \
      : \
      "m" (*((s)+3)), \
      "m" (*((s)+6)), \
      "m" (*((s)+7)), \
      "m" (_sse_double_sgn2)); \
_ASM ("movapd %%xmm0, %0 \n\t"\
      "movapd %%xmm1, %1 \n\t"\
      "movapd %%xmm2, %2 \n\t"\
      : \
      "=m" (*((r)+1)), \
      "=m" (*((r)+2)), \
      "=m" (*((r)+5))); } \

// //////////////////////////////////////////////////////////////////////////
// r[0].c1=s[0].c1, r[0].2=s[0].c2+...+r[7].c2=s[7].c2;
// //////////////////////////////////////////////////////////////////////////

#define _sse_double_copy_16(r,s) { \
_ASM ("movapd %0, %%xmm0 \n\t" \
      "movapd %1, %%xmm1 \n\t" \
      "movapd %2, %%xmm2 \n\t" \
      "movapd %3, %%xmm3 \n\t" \
      : \
      : \
      "m" (*((s))), \
      "m" (*((s)+1)), \
      "m" (*((s)+2)), \
      "m" (*((s)+3)));	       \
_ASM ("movapd %0, %%xmm4 \n\t"	\
      "movapd %1, %%xmm5 \n\t" \
      "movapd %2, %%xmm6 \n\t" \
      "movapd %3, %%xmm7 \n\t" \
      : \
      : \
      "m" (*((s)+4)), \
      "m" (*((s)+5)), \
      "m" (*((s)+6)), \
      "m" (*((s)+7))); \
_ASM ("movapd %%xmm0, %0 \n\t" \
      "movapd %%xmm1, %1 \n\t" \
      "movapd %%xmm2, %2 \n\t" \
      "movapd %%xmm3, %3 \n\t" \
      : \
      "=m" (*((r))), \
      "=m" (*((r)+1)), \
      "=m" (*((r)+2)), \
      "=m" (*((r)+3)));	       \
_ASM ("movapd %%xmm4, %0 \n\t"	\
      "movapd %%xmm5, %1 \n\t" \
      "movapd %%xmm6, %2 \n\t" \
      "movapd %%xmm7, %3 \n\t" \
      : \
      "=m" (*((r)+4)), \
      "=m" (*((r)+5)), \
      "=m" (*((r)+6)), \
      "=m" (*((r)+7))); }

// //////////////////////////////////////////////////////////////////////////
// r[0].c1+=s[0].c1, r[0].2+=s[0].c2 ... r[7].c2+=s[7].c2;
// //////////////////////////////////////////////////////////////////////////

#define _sse_double_add_16(r,s) { \
_ASM ("movapd %0, %%xmm0 \n\t" \
     "movapd %1, %%xmm1 \n\t" \
     "movapd %2, %%xmm2 \n\t" \
     "movapd %3, %%xmm3 \n\t" \
     : \
     : \
     "m" (*((s))), \
     "m" (*((s)+1)), \
     "m" (*((s)+2)), \
     "m" (*((s)+3)));	      \
_ASM ("movapd %0, %%xmm4 \n\t"		\
     "movapd %1, %%xmm5 \n\t" \
     "movapd %2, %%xmm6 \n\t" \
     "movapd %3, %%xmm7 \n\t" \
     : \
     : \
     "m" (*((s)+4)), \
     "m" (*((s)+5)), \
     "m" (*((s)+6)), \
     "m" (*((s)+7))); \
_ASM ("addpd %0, %%xmm0 \n\t" \
      "addpd %1, %%xmm1 \n\t" \
      "addpd %2, %%xmm2 \n\t" \
      "addpd %3, %%xmm3 \n\t" \
      : \
      : \
      "m" (*((r))), \
      "m" (*((r)+1)), \
      "m" (*((r)+2)), \
      "m" (*((r)+3)));	      \
_ASM ("addpd %0, %%xmm4 \n\t"		\
      "addpd %1, %%xmm5 \n\t" \
      "addpd %2, %%xmm6 \n\t" \
      "addpd %3, %%xmm7 \n\t" \
      : \
      : \
      "m" (*((r)+4)), \
      "m" (*((r)+5)), \
      "m" (*((r)+6)), \
      "m" (*((r)+7))); \
_ASM ("movapd %%xmm0, %0 \n\t" \
      "movapd %%xmm1, %1 \n\t" \
      "movapd %%xmm2, %2 \n\t" \
      "movapd %%xmm3, %3 \n\t" \
      :			       \
      "=m" (*((r))),	       \
      "=m" (*((r)+1)),	       \
      "=m" (*((r)+2)),	       \
      "=m" (*((r)+3)));	       \
_ASM ("movapd %%xmm4, %0 \n\t"		\
      "movapd %%xmm5, %1 \n\t"			\
      "movapd %%xmm6, %2 \n\t"			\
      "movapd %%xmm7, %3 \n\t"			\
      :						\
      "=m" (*((r)+4)),				\
      "=m" (*((r)+5)),				\
      "=m" (*((r)+6)),				\
      "=m" (*((r)+7))); }

// //////////////////////////////////////////////////////////////////////////
// r[0].c1-=s[0].c1, r[0].2-=s[0].c2 ... r[7].c2-=s[7].c2; CHECK MAY BE WRONG
// //////////////////////////////////////////////////////////////////////////

#define _sse_double_sub_16(r,s) { \
_ASM ("movapd %0, %%xmm0 \n\t" \
     "movapd %1, %%xmm1 \n\t" \
     "movapd %2, %%xmm2 \n\t" \
     "movapd %3, %%xmm3 \n\t" \
     : \
     : \
     "m" (*((s))), \
     "m" (*((s)+1)), \
     "m" (*((s)+2)), \
      "m" (*((s)+3)));	      \
_ASM ("movapd %0, %%xmm4 \n\t"		\
     "movapd %1, %%xmm5 \n\t" \
     "movapd %2, %%xmm6 \n\t" \
     "movapd %3, %%xmm7 \n\t" \
     : \
     : \
     "m" (*((s)+4)), \
     "m" (*((s)+5)), \
     "m" (*((s)+6)), \
     "m" (*((s)+7))); \
_ASM ("subpd %0, %%xmm0 \n\t" \
      "subpd %1, %%xmm1 \n\t" \
      "subpd %2, %%xmm2 \n\t" \
      "subpd %3, %%xmm3 \n\t" \
      : \
      : \
      "m" (*((r))), \
      "m" (*((r)+1)), \
      "m" (*((r)+2)), \
      "m" (*((r)+3)));	      \
_ASM ("subpd %0, %%xmm4 \n\t"		\
      "subpd %1, %%xmm5 \n\t" \
      "subpd %2, %%xmm6 \n\t" \
      "subpd %3, %%xmm7 \n\t" \
      : \
      : \
      "m" (*((r)+4)), \
      "m" (*((r)+5)), \
      "m" (*((r)+6)), \
      "m" (*((r)+7))); \
_ASM ("movapd %%xmm0, %0 \n\t" \
     "movapd %%xmm1, %1 \n\t" \
     "movapd %%xmm2, %2 \n\t" \
     "movapd %%xmm3, %3 \n\t" \
     : \
     "=m" (*((r))), \
     "=m" (*((r)+1)), \
     "=m" (*((r)+2)), \
      "=m" (*((r)+3)));	      \
_ASM ("movapd %%xmm4, %0 \n\t"		\
     "movapd %%xmm5, %1 \n\t" \
     "movapd %%xmm6, %2 \n\t" \
     "movapd %%xmm7, %3 \n\t" \
     : \
     "=m" (*((r)+4)), \
     "=m" (*((r)+5)), \
     "=m" (*((r)+6)), \
     "=m" (*((r)+7))); }

// //////////////////////////////////////////////////////////////////////////
// r[0].c1=c*s[0].c1, r[0].2=c*s[0].c2 ... r[7].c2=c*s[7].c2;
// //////////////////////////////////////////////////////////////////////////

#define  _sse_double_add_multiply_16(r,c,s) { \
_ASM ("movapd %0, %%xmm0 \n\t" \
     "movapd %1, %%xmm1 \n\t" \
     "movapd %2, %%xmm2 \n\t" \
     "movapd %3, %%xmm3 \n\t" \
     : \
     : \
     "m" (*((s))), \
     "m" (*((s)+1)), \
     "m" (*((s)+2)), \
      "m" (*((s)+3)));	      \
_ASM ("movapd %0, %%xmm4 \n\t"		\
     "movapd %1, %%xmm5 \n\t" \
     "movapd %2, %%xmm6 \n\t" \
     "movapd %3, %%xmm7 \n\t" \
     "mulpd %4, %%xmm0 \n\t" \
     "mulpd %4, %%xmm1 \n\t" \
     "mulpd %4, %%xmm2 \n\t" \
     "mulpd %4, %%xmm3 \n\t" \
     "mulpd %4, %%xmm4 \n\t" \
     "mulpd %4, %%xmm5 \n\t" \
     "mulpd %4, %%xmm6 \n\t" \
     "mulpd %4, %%xmm7 \n\t" \
     : \
     : \
     "m" (*((s)+4)), \
     "m" (*((s)+5)), \
     "m" (*((s)+6)), \
     "m" (*((s)+7)), \
     "m" (c)); \
_ASM ("addpd %0, %%xmm0 \n\t" \
      "addpd %1, %%xmm1 \n\t" \
      "addpd %2, %%xmm2 \n\t" \
      "addpd %3, %%xmm3 \n\t" \
      : \
      : \
      "m" (*((r))), \
      "m" (*((r)+1)), \
      "m" (*((r)+2)), \
      "m" (*((r)+3)));	      \
_ASM ("addpd %0, %%xmm4 \n\t"		\
      "addpd %1, %%xmm5 \n\t" \
      "addpd %2, %%xmm6 \n\t" \
      "addpd %3, %%xmm7 \n\t" \
      : \
      : \
      "m" (*((r)+4)), \
      "m" (*((r)+5)), \
      "m" (*((r)+6)), \
      "m" (*((r)+7))); \
_ASM ("movapd %%xmm0, %0 \n\t" \
      "movapd %%xmm1, %1 \n\t" \
      "movapd %%xmm2, %2 \n\t" \
      "movapd %%xmm3, %3 \n\t" \
      : \
      "=m" (*((r))), \
      "=m" (*((r)+1)), \
      "=m" (*((r)+2)), \
      "=m" (*((r)+3)));	       \
_ASM ("movapd %%xmm4, %0 \n\t"	\
      "movapd %%xmm5, %1 \n\t" \
      "movapd %%xmm6, %2 \n\t" \
      "movapd %%xmm7, %3 \n\t" \
      : \
      "=m" (*((r)+4)), \
      "=m" (*((r)+5)), \
      "=m" (*((r)+6)), \
      "=m" (*((r)+7))); }

#define  _sse_double_multiply_16(r,c,s) { \
_ASM ("movapd %0, %%xmm0 \n\t" \
      "movapd %1, %%xmm1 \n\t" \
      "movapd %2, %%xmm2 \n\t" \
      "movapd %3, %%xmm3 \n\t" \
      : \
      : \
      "m" (*((s))), \
      "m" (*((s)+1)), \
      "m" (*((s)+2)), \
      "m" (*((s)+3)));	       \
_ASM ("movapd %0, %%xmm4 \n\t"	\
      "movapd %1, %%xmm5 \n\t" \
      "movapd %2, %%xmm6 \n\t" \
      "movapd %3, %%xmm7 \n\t" \
      : \
      : \
      "m" (*((s)+4)), \
      "m" (*((s)+5)), \
      "m" (*((s)+6)), \
      "m" (*((s)+7))); \
_ASM ("mulpd  %0, %%xmm0 \n\t" \
      "mulpd  %0, %%xmm1 \n\t" \
      "mulpd  %0, %%xmm2 \n\t" \
      "mulpd  %0, %%xmm3 \n\t" \
      "mulpd  %0, %%xmm4 \n\t" \
      "mulpd  %0, %%xmm5 \n\t" \
      "mulpd  %0, %%xmm6 \n\t" \
      "mulpd  %0, %%xmm7 \n\t" \
      : \
      : \
      "m" (c)); \
_ASM ("movapd %%xmm0, %0 \n\t" \
      "movapd %%xmm1, %1 \n\t" \
      "movapd %%xmm2, %2 \n\t" \
      "movapd %%xmm3, %3 \n\t" \
      : \
      "=m" (*((r))), \
      "=m" (*((r)+1)), \
      "=m" (*((r)+2)), \
      "=m" (*((r)+3)));	       \
_ASM ("movapd %%xmm4, %0 \n\t"	\
      "movapd %%xmm5, %1 \n\t" \
      "movapd %%xmm6, %2 \n\t" \
      "movapd %%xmm7, %3 \n\t" \
      : \
      "=m" (*((r)+4)), \
      "=m" (*((r)+5)), \
      "=m" (*((r)+6)), \
      "=m" (*((r)+7))); }


static void _sse_check_alignment(void* var, unsigned int base) {
  unsigned int af1=(unsigned int) var;
  if (af1!=(af1&~base)) {
    error("_sse_check_alignment()\nVariable not aligned properly");
  }
}

