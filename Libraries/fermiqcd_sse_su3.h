/////////////////////////////////////////////////////////////////
/// @file fermiqcd_sse_su3.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Stuff for SSE/SSE2 compile with -DSSE2
///
/// Distributed under GPL2 License
/// 
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////

#if defined(USE_DOUBLE_PRECISION) && defined(SSE2)

inline void _sse_mulABC_set_331(mdp_complex* a, mdp_complex* b, mdp_complex* c) 
    {                                             
      _sse_double_load_123(b[0],b[1],b[2]);       
      _sse_double_su3_multiply(*((_sse_su3*) a)); 
      _sse_double_store_up_123(c[0],c[1],c[2]);      
    }
inline void _sse_mulABC_add_331(mdp_complex* a, mdp_complex* b, mdp_complex* c) 
    {                                             
      _sse_double_load_123(b[0],b[1],b[2]);       
      _sse_double_su3_multiply(*((_sse_su3*) a)); 
      _sse_double_load_123(c[0],c[1],c[2]);       
      _sse_double_vector_add();                   
      _sse_double_store_123(c[0],c[1],c[2]);      
    }
inline void _sse_mulABC_sub_331(mdp_complex* a, mdp_complex* b, mdp_complex* c) 
    {                                             
      _sse_double_load_123(b[0],b[1],b[2]);       
      _sse_double_su3_multiply(*((_sse_su3*) a)); 
      _sse_double_load_123(c[0],c[1],c[2]);       
      _sse_double_vector_sub();                   
      _sse_double_store_123(c[0],c[1],c[2]);      
    }

inline void _sse_mulAHBC_set_331(mdp_complex* a, mdp_complex* b, mdp_complex* c) 
    {                                             
      _sse_double_load_123(b[0],b[1],b[2]);       
      _sse_double_su3_inverse_multiply(*((_sse_su3*) a)); 
      _sse_double_store_up_123(c[0],c[1],c[2]);      
    }
inline void _sse_mulAHBC_add_331(mdp_complex* a, mdp_complex* b, mdp_complex* c) 
    {                                             
      _sse_double_load_123(b[0],b[1],b[2]);       
      _sse_double_su3_inverse_multiply(*((_sse_su3*) a)); 
      _sse_double_load_123(c[0],c[1],c[2]);       
      _sse_double_vector_add();                   
      _sse_double_store_123(c[0],c[1],c[2]);      
    }
inline void _sse_mulAHBC_sub_331(mdp_complex* a, mdp_complex* b, mdp_complex* c) 
    {                                             
      _sse_double_load_123(b[0],b[1],b[2]);       
      _sse_double_su3_inverse_multiply(*((_sse_su3*) a)); 
      _sse_double_load_123(c[0],c[1],c[2]);       
      _sse_double_vector_sub();                   
      _sse_double_store_123(c[0],c[1],c[2]);      
    }


inline void _sse_mulABC_set_333(mdp_complex* a, mdp_complex* b, mdp_complex* c) 
    {                                             
      _sse_double_load_123(b[0],b[3],b[6]);       
      _sse_double_su3_multiply(*((_sse_su3*) a)); 
      _sse_double_store_up_123(c[0],c[3],c[6]);      
      _sse_double_load_123(b[1],b[4],b[7]);       
      _sse_double_su3_multiply(*((_sse_su3*) a)); 
      _sse_double_store_up_123(c[1],c[4],c[7]);      
      _sse_double_load_123(b[2],b[5],b[8]);       
      _sse_double_su3_multiply(*((_sse_su3*) a)); 
      _sse_double_store_up_123(c[2],c[5],c[8]);      
    }

inline void _sse_mulABC_add_333(mdp_complex* a, mdp_complex* b, mdp_complex* c) 
    {                                             
      _sse_double_load_123(b[0],b[3],b[6]);       
      _sse_double_su3_multiply(*((_sse_su3*) a)); 
      _sse_double_load_123(c[0],c[3],c[6]);       
      _sse_double_vector_add();                   
      _sse_double_store_123(c[0],c[3],c[6]);      
      _sse_double_load_123(b[1],b[4],b[7]);       
      _sse_double_su3_multiply(*((_sse_su3*) a)); 
      _sse_double_load_123(c[1],c[4],c[7]);       
      _sse_double_vector_add();                   
      _sse_double_store_123(c[1],c[4],c[7]);      
      _sse_double_load_123(b[2],b[5],b[8]);       
      _sse_double_su3_multiply(*((_sse_su3*) a)); 
      _sse_double_load_123(c[2],c[5],c[8]);       
      _sse_double_vector_add();                   
      _sse_double_store_123(c[2],c[5],c[8]);      
    }

inline void _sse_mulABC_sub_333(mdp_complex* a, mdp_complex* b, mdp_complex* c) 
    {                                             
      _sse_double_load_123(b[0],b[3],b[6]);       
      _sse_double_su3_multiply(*((_sse_su3*) a)); 
      _sse_double_load_123(c[0],c[3],c[6]);       
      _sse_double_vector_sub();                   
      _sse_double_store_123(c[0],c[3],c[6]);      
      _sse_double_load_123(b[1],b[4],b[7]);       
      _sse_double_su3_multiply(*((_sse_su3*) a)); 
      _sse_double_load_123(c[1],c[4],c[7]);       
      _sse_double_vector_sub();                   
      _sse_double_store_123(c[1],c[4],c[7]);      
      _sse_double_load_123(b[2],b[5],b[8]);       
      _sse_double_su3_multiply(*((_sse_su3*) a)); 
      _sse_double_load_123(c[2],c[5],c[8]);       
      _sse_double_vector_sub();                   
      _sse_double_store_123(c[2],c[5],c[8]);      
    }

inline void _sse_mulAHBC_set_333(mdp_complex* a, mdp_complex* b, mdp_complex* c) 
    {                                             
      _sse_double_load_123(b[0],b[3],b[6]);       
      _sse_double_su3_inverse_multiply(*((_sse_su3*) a)); 
      _sse_double_store_up_123(c[0],c[3],c[6]);      
      _sse_double_load_123(b[1],b[4],b[7]);       
      _sse_double_su3_inverse_multiply(*((_sse_su3*) a)); 
      _sse_double_store_up_123(c[1],c[4],c[7]);      
      _sse_double_load_123(b[2],b[5],b[8]);       
      _sse_double_su3_inverse_multiply(*((_sse_su3*) a)); 
      _sse_double_store_up_123(c[2],c[5],c[8]);      
    }

inline void _sse_mulAHBC_add_333(mdp_complex* a, mdp_complex* b, mdp_complex* c) 
    {                                             
      _sse_double_load_123(b[0],b[3],b[6]);       
      _sse_double_su3_inverse_multiply(*((_sse_su3*) a)); 
      _sse_double_load_123(c[0],c[3],c[6]);       
      _sse_double_vector_add();                   
      _sse_double_store_123(c[0],c[3],c[6]);      
      _sse_double_load_123(b[1],b[4],b[7]);       
      _sse_double_su3_inverse_multiply(*((_sse_su3*) a)); 
      _sse_double_load_123(c[1],c[4],c[7]);       
      _sse_double_vector_add();                   
      _sse_double_store_123(c[1],c[4],c[7]);      
      _sse_double_load_123(b[2],b[5],b[8]);       
      _sse_double_su3_inverse_multiply(*((_sse_su3*) a)); 
      _sse_double_load_123(c[2],c[5],c[8]);       
      _sse_double_vector_add();                   
      _sse_double_store_123(c[2],c[5],c[8]);      
    }

inline void _sse_mulAHBC_sub_333(mdp_complex* a, mdp_complex* b, mdp_complex* c) 
    {                                             
      _sse_double_load_123(b[0],b[3],b[6]);       
      _sse_double_su3_inverse_multiply(*((_sse_su3*) a)); 
      _sse_double_load_123(c[0],c[3],c[6]);       
      _sse_double_vector_sub();                   
      _sse_double_store_123(c[0],c[3],c[6]);      
      _sse_double_load_123(b[1],b[4],b[7]);       
      _sse_double_su3_inverse_multiply(*((_sse_su3*) a)); 
      _sse_double_load_123(c[1],c[4],c[7]);       
      _sse_double_vector_sub();                   
      _sse_double_store_123(c[1],c[4],c[7]);      
      _sse_double_load_123(b[2],b[5],b[8]);       
      _sse_double_su3_inverse_multiply(*((_sse_su3*) a)); 
      _sse_double_load_123(c[2],c[5],c[8]);       
      _sse_double_vector_sub();                   
      _sse_double_store_123(c[2],c[5],c[8]);      
    }


inline void _sse_mulABHC_set_333(mdp_complex* a, mdp_complex* b, mdp_complex* c) 
    {                          
      static _sse_su3_vector v ALIGN16;
      v.c1=conj(b[0]); v.c2=conj(b[1]); v.c3=conj(b[2]);
      _sse_double_load(v);       
      _sse_double_su3_multiply(*((_sse_su3*) a)); 
      _sse_double_store_up_123(c[0],c[3],c[6]);      

      v.c1=conj(b[3]); v.c2=conj(b[4]); v.c3=conj(b[5]);
      _sse_double_load(v);       
      _sse_double_su3_multiply(*((_sse_su3*) a)); 
      _sse_double_store_up_123(c[1],c[4],c[7]);      

      v.c1=conj(b[6]); v.c2=conj(b[7]); v.c3=conj(b[8]);
      _sse_double_load(v);       
      _sse_double_su3_multiply(*((_sse_su3*) a)); 
      _sse_double_store_up_123(c[2],c[5],c[8]);      
    }

inline void _sse_mulABHC_add_333(mdp_complex* a, mdp_complex* b, mdp_complex* c) 
    {  
      static _sse_su3_vector v ALIGN16;
      v.c1=conj(b[0]); v.c2=conj(b[1]); v.c3=conj(b[2]);
      _sse_double_load(v);       
      _sse_double_su3_multiply(*((_sse_su3*) a)); 
      _sse_double_load_123(c[0],c[3],c[6]);       
      _sse_double_vector_add();                   
      _sse_double_store_123(c[0],c[3],c[6]);      

      v.c1=conj(b[3]); v.c2=conj(b[4]); v.c3=conj(b[5]);
      _sse_double_load(v);
      _sse_double_su3_multiply(*((_sse_su3*) a)); 
      _sse_double_load_123(c[1],c[4],c[7]);       
      _sse_double_vector_add();                   
      _sse_double_store_123(c[1],c[4],c[7]);      

      v.c1=conj(b[6]); v.c2=conj(b[7]); v.c3=conj(b[8]);
      _sse_double_load(v);       
      _sse_double_su3_multiply(*((_sse_su3*) a)); 
      _sse_double_load_123(c[2],c[5],c[8]);       
      _sse_double_vector_add();                   
      _sse_double_store_123(c[2],c[5],c[8]);      
    }

inline void _sse_mulABHC_sub_333(mdp_complex* a, mdp_complex* b, mdp_complex* c) 
    {                             
      static _sse_su3_vector v ALIGN16;
      v.c1=conj(b[0]); v.c2=conj(b[1]); v.c3=conj(b[2]);
      _sse_double_load(v);       
      _sse_double_su3_multiply(*((_sse_su3*) a)); 
      _sse_double_load_123(c[0],c[3],c[6]);       
      _sse_double_vector_sub();                   
      _sse_double_store_123(c[0],c[3],c[6]);      

      v.c1=conj(b[3]); v.c2=conj(b[4]); v.c3=conj(b[5]);
      _sse_double_load(v);
      _sse_double_su3_multiply(*((_sse_su3*) a)); 
      _sse_double_load_123(c[1],c[4],c[7]);       
      _sse_double_vector_sub();                   
      _sse_double_store_123(c[1],c[4],c[7]);      

      v.c1=conj(b[6]); v.c2=conj(b[7]); v.c3=conj(b[8]);
      _sse_double_load(v);       
      _sse_double_su3_multiply(*((_sse_su3*) a)); 
      _sse_double_load_123(c[2],c[5],c[8]);       
      _sse_double_vector_sub();                   
      _sse_double_store_123(c[2],c[5],c[8]);      
    }

inline void _sse_mulAbC_set_31(mdp_complex* a, mdp_complex b, mdp_complex* c) {
  
  static _sse_double real, imag ALIGN16;
  real.c1=real.c2=b.real();
  imag.c1=imag.c2=b.imag();
  _sse_double_load_123(a[0],a[1],a[2]);       
  _sse_double_vector_mul_complex(real,imag);
  _sse_double_store_up_123(c[0],c[1],c[2]);   
}

inline void _sse_mulAbC_add_31(mdp_complex* a, mdp_complex b, mdp_complex* c) 
    {                                             
      static _sse_double real, imag ALIGN16;
      real.c1=real.c2=b.real();
      imag.c1=imag.c2=b.imag();
      _sse_double_load_123(a[0],a[1],a[2]);       
      _sse_double_vector_mul_complex(real, imag);
      _sse_double_load_123(c[0],c[1],c[2]);          
      _sse_double_vector_add();                   
      _sse_double_store_123(c[0],c[1],c[2]);       
    }
inline void _sse_mulAbC_sub_31(mdp_complex* a, mdp_complex b, mdp_complex* c) 
    {                                            
      static _sse_double real, imag ALIGN16;
      real.c1=real.c2=b.real();
      imag.c1=imag.c2=b.imag();          
      _sse_double_load_123(a[0],a[1],a[2]);       
      _sse_double_vector_mul_complex(real, imag);                  
      _sse_double_load_123(c[0],c[1],c[2]);          
      _sse_double_vector_sub();                   
      _sse_double_store_123(c[0],c[1],c[2]);       
    }

inline void _sse_sumAC_set_31(mdp_complex* a, mdp_complex* c) 
    {                                      
      _sse_double_load_123(a[0],a[1],a[2]); 
      _sse_double_store_123(c[0],c[1],c[2]); 
    }
inline void _sse_sumAC_add_31(mdp_complex* a, mdp_complex* c) 
    {                                          
      _sse_double_load_123(c[0],c[1],c[2]);    
      _sse_double_load_up_123(a[0],a[1],a[2]); 
      _sse_double_vector_add();                
      _sse_double_store_123(c[0],c[1],c[2]);   
    }
inline void _sse_sumAC_sub_31(mdp_complex* a, mdp_complex* c)
    {                                          
      _sse_double_load_123(c[0],c[1],c[2]);    
      _sse_double_load_up_123(a[0],a[1],a[2]); 
      _sse_double_vector_sub();                
      _sse_double_store_123(c[0],c[1],c[2]);   
    }

#endif
