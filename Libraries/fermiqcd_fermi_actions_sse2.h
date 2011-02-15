/////////////////////////////////////////////////////////////////
/// @file fermiqcd_fermi_actions_sse2.h
/// @version 2009-12-21
/// @author Martin Luescher and Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Basic actions for Wilson Fermions optimized in assembler
///
//////////////////////////////////////////////////////////////////


#if defined(SSE2)

/// @brief Wilson/Clover action SSE/SSE2
///
/// Only on Pentium 4 or compatible. Compile with -DSSE2
/// Example:
/// @verbatim
/// gauge_field U(lattice,nc);
/// fermi_field psi(lattice,nc);
/// fermi_field chi(lattice,nc);
/// coefficients coeff;
/// coeff["kappa_s"]=0.11;
/// coeff["kappa_t"]=0.11;
/// coeff["r_s"]=1.0;
/// coeff["r_t"]=1.0;
/// coeff["c_{sw}"]=1.0;
/// coeff["c_E"]=1.0;
/// coeff["c_B"]=1.0;
/// default_fermi_action=FermiCloverActionSSE2::mul_Q;
/// if(coeff["c_{sw}"]!=0) compute_em_field(U);
/// mul_Q(chi,psi,U,coeff);
/// @endverbatim
/// Attention: uses always FERMILAB convention of Gamma matrices
class FermiCloverActionSSE2 {
 public:
  static void mul_Q(fermi_field &chi_out, 
		    fermi_field &psi_in, 
		    gauge_field &U_in, 
		    coefficients &coeff,
		    int parity=EVENODD) { 

    register int   ndim=psi_in.lattice().ndim;
    register int   nspin=psi_in.nspin;
    register int   nc=psi_in.nc;
    register mdp_real kappa_t=0;
    register mdp_real kappa_s=0;
    register mdp_real r_t;
    register mdp_real r_s;
    register mdp_real cSW;
    register mdp_real c_E;
    register mdp_real c_B;
    register int sign;

    if(coeff.has_key("kappa")) kappa_s=kappa_t=coeff["kappa"];
    if(coeff.has_key("kappa_t")) kappa_t=coeff["kappa_t"];
    if(coeff.has_key("kappa_s")) kappa_s=coeff["kappa_s"];
    if(kappa_t==0) error("FermiCloverActionSSE2\nkappa_t=0 or undeclared");
    if(kappa_s==0) error("FermiCloverActionSSE2\nkappa_s=0 or undeclared");
    if(coeff.has_key("r_t")) r_t=coeff["r_t"];       else r_t=1;
    if(coeff.has_key("r_s")) r_s=coeff["r_s"];       else r_s=1;
    if(coeff.has_key("c_{sw}")) cSW=coeff["c_{sw}"]; else cSW=0;
    if(coeff.has_key("c_E")) c_E=coeff["c_E"];       else c_E=1;
    if(coeff.has_key("c_B")) c_B=coeff["c_B"];       else c_B=1;
    if(coeff.has_key("sign")) sign=(int) coeff["sign"];   else sign=+1;
    if(parity!=EVENODD) error("FermiCloverActionSSE2\nparity must be EVENODD here");

#if !defined(USE_DOUBLE_PRECISION)

    _sse_spinor *chi=(_sse_spinor*) chi_out.physical_address();
    _sse_check_alignment((void*) chi, 0xf);
    _sse_spinor *psi=(_sse_spinor*) psi_in.physical_address();
    _sse_check_alignment((void*) psi, 0xf);
    _sse_su3    *U=(_sse_su3*) U_in.physical_address();
    _sse_check_alignment((void*) U, 0xf);
    _sse_su3    *uem=0;
    if(cSW!=0.0) {
      uem=(_sse_su3*) U_in.em.physical_address();
      _sse_check_alignment((void*) uem, 0xf);
    }
    mdp_int   **iup=U_in.lattice().up;
    mdp_int   **idw=U_in.lattice().dw;
    mdp_int   start=U_in.lattice().start_index(ME,0); // even
    mdp_int   stop =U_in.lattice().stop_index(ME,1);  // odd
    
    _sse_float fact1 ALIGN16;
    _sse_float fact2 ALIGN16;
    _sse_float fact3 ALIGN16;
    _sse_float fact4 ALIGN16;
    _sse_float fact5 ALIGN16;
    _sse_float fact6 ALIGN16;
    _sse_vector r12_1  ALIGN16;
    _sse_vector r34_1  ALIGN16;
    _sse_vector r12_2  ALIGN16;
    _sse_vector r34_2  ALIGN16;
    _sse_vector r0  ALIGN16;
    mdp_int ix1,iy1,iy2,iz1;
    float rho;
    _sse_su3 *up1,*um1,*um2;   
    _sse_spinor *s1,*sp1,*sp2,*sm1,*sm2,*sn1;
    
    if(sign!=1) exit(1);

    if((stop-start)%2 !=0) 
      error("FermiCloverActionSSE2\nProblem with parallelization: odd # of sites on process!");
    
    if(r_t!=1.0) 
      error("FermiCloverActionSSE2\nr_t!=1 not compatible with SSE2\n");

    _sse_check_alignment((void*) &fact1, 0xf);
    _sse_check_alignment((void*) &fact2, 0xf);
    _sse_check_alignment((void*) &fact3, 0xf);
    _sse_check_alignment((void*) &fact4, 0xf);
    _sse_check_alignment((void*) &fact5, 0xf);
    _sse_check_alignment((void*) &fact6, 0xf);
    _sse_check_alignment((void*) &r12_1, 0xf);
    _sse_check_alignment((void*) &r34_1, 0xf);
    _sse_check_alignment((void*) &r12_2, 0xf);
    _sse_check_alignment((void*) &r34_2, 0xf);
    
    /* empty vector */
    
    r0.c1.c1=r0.c1.c2=r0.c1.c3=r0.c1.c4=0;
    r0.c2.c1=r0.c2.c2=r0.c2.c3=r0.c2.c4=0;
    r0.c3.c1=r0.c3.c2=r0.c3.c3=r0.c3.c4=0;
    
    rho=-1.0f/kappa_s;
    
    /* coefficient of (1 +- Gamma[i]) */
    
    fact1.c1=rho;
    fact1.c2=rho;
    fact1.c3=rho;
    fact1.c4=rho;
    
    fact2.c1=-1.0f*kappa_s;
    fact2.c2=fact2.c1;
    fact2.c3=fact2.c1;
    fact2.c4=fact2.c1;
    
    /* coefficient of rho/2*(1 +- Gamma[0]) */
    
    fact3.c1=(1.0f+r_t)*kappa_t/kappa_s;
    fact3.c2=fact3.c1;
    fact3.c3=fact3.c1;
    fact3.c4=fact3.c1;
    
    fact4.c1=-1.0f;
    fact4.c2=-1.0f;
    fact4.c3=-1.0f;
    fact4.c4=-1.0f;
    
    /* coefficient of Sigma[i][j]*U.em(x,i,j) */
    
    fact5.c1=1.0f*kappa_s*cSW*c_B;
    fact5.c2=fact5.c1;
    fact5.c3=fact5.c1;
    fact5.c4=fact5.c1;

    /* coefficient of fact5*Sigma[0][i]*U.em(x,0,i) */
    
    fact6.c1=1.0f*c_E/c_B;
    fact6.c2=fact6.c1;
    fact6.c3=fact6.c1;
    fact6.c4=fact6.c1;
    
    sp1=(_sse_spinor*) &psi[iup[start][0]];
    sp2=(_sse_spinor*) &psi[iup[start+1][0]];
    up1=(_sse_su3*) U+4*start;
    
    /************************ loop over all lattice sites ***************/
    
    
    for (ix1=start; ix1<stop; ix1+=2) {   


      s1=psi+ix1;
      _sse_float_prefetch_spinor(s1);


      /******************************* direction +0 ***********************/
      
      iy1=idw[ix1][0];
      iy2=idw[ix1+1][0];
      sm1=psi+iy1;
      sm2=psi+iy2;
      _sse_float_prefetch_spinor(sm1);
      _sse_float_prefetch_spinor(sm2);
      
      _sse_float_pair_load((*sp1).c3,(*sp1).c4);
      _sse_float_vector_mul(fact3);

      _sse_float_su3_multiply((*up1));

      _sse_float_pair_load((*s1).c1,(*s1).c2);
      _sse_float_vector_mul(fact1);      
      _sse_float_vector_store(r12_1);
      
      _sse_float_pair_load((*s1).c3,(*s1).c4);
      _sse_float_vector_mul(fact1);      
      _sse_float_vector_add();
      _sse_float_vector_store(r34_1);
      
      um1=U+iy1*4;
      _sse_float_prefetch_su3(um1);
      um2=U+iy2*4;
      _sse_float_prefetch_su3(um2);
      
      _sse_float_pair_load((*(sp2)).c3,(*(sp2)).c4);
      _sse_float_vector_mul(fact3);
      
      _sse_float_su3_multiply((*(up1+4)));
      
      _sse_float_pair_load((*(s1+1)).c1,(*(s1+1)).c2);
      _sse_float_vector_mul(fact1);      
      _sse_float_vector_store(r12_2);
    
    _sse_float_pair_load((*(s1+1)).c3,(*(s1+1)).c4);
    _sse_float_vector_mul(fact1);      
    _sse_float_vector_add();
    _sse_float_vector_store(r34_2);

    
    /******************************* direction -0 ***********************/

    sp1=psi+iup[ix1][1];
    sp2=psi+iup[ix1+1][1];
    _sse_float_prefetch_spinor(sp1);
    _sse_float_prefetch_spinor(sp2);
        
    _sse_float_pair_load((*sm1).c1,(*sm1).c2);
    _sse_float_vector_mul(fact3);      
    
    _sse_float_su3_inverse_multiply((*um1));
    
    _sse_float_vector_load(r12_1);
    _sse_float_vector_add();
    _sse_float_vector_store(r12_1);
    
    up1++;
    _sse_float_prefetch_su3(up1);      

    _sse_float_pair_load((*(sm2)).c1,(*(sm2)).c2);
    _sse_float_vector_mul(fact3);      
    
    _sse_float_su3_inverse_multiply((*(um2)));
    
    _sse_float_vector_load(r12_2);
    _sse_float_vector_add();
    _sse_float_vector_store(r12_2);
      
    /******************************* direction +1 ***********************/
    
    iy1=idw[ix1][1];
    iy2=idw[ix1+1][1];
    sm1=psi+iy1;
    sm2=psi+iy2;
    _sse_float_prefetch_spinor(sm1);
    _sse_float_prefetch_spinor(sm2);
    
    _sse_float_pair_load((*sp1).c1,(*sp1).c2);
    _sse_float_pair_load_up((*sp1).c4,(*sp1).c3);
    _sse_float_vector_sub();
    
    _sse_float_su3_multiply((*up1));
    
    _sse_float_vector_load(r12_1);
    _sse_float_vector_add();
    _sse_float_vector_store(r12_1);
    
    _sse_float_vector_load(r34_1);
    _sse_float_vector_xch();
    _sse_float_vector_sub();
    _sse_float_vector_store(r34_1);
    
    um1=U+iy1*4+1;
    _sse_float_prefetch_su3(um1);      
    um2=U+iy2*4+1;
    _sse_float_prefetch_su3(um2);      
    
    _sse_float_pair_load((*(sp2)).c1,(*(sp2)).c2);
    _sse_float_pair_load_up((*(sp2)).c4,(*(sp2)).c3);
    _sse_float_vector_sub();
    
    _sse_float_su3_multiply((*(up1+4)));
    
    _sse_float_vector_load(r12_2);
    _sse_float_vector_add();
    _sse_float_vector_store(r12_2);
    
    _sse_float_vector_load(r34_2);
    _sse_float_vector_xch();
    _sse_float_vector_sub();
    _sse_float_vector_store(r34_2);
    
    /******************************* direction -1 ***********************/
    
    sp1=psi+iup[ix1][2];
    sp2=psi+iup[ix1+1][2];
    _sse_float_prefetch_spinor(sp1);
    _sse_float_prefetch_spinor(sp2);
    
    _sse_float_pair_load((*sm1).c1,(*sm1).c2);
    _sse_float_pair_load_up((*sm1).c4,(*sm1).c3);
    _sse_float_vector_add();
    
    _sse_float_su3_inverse_multiply((*um1));
    
    _sse_float_vector_load(r12_1);
    _sse_float_vector_add();
    _sse_float_vector_store(r12_1);      

    _sse_float_vector_load(r34_1);
    _sse_float_vector_xch();
    _sse_float_vector_add();
    _sse_float_vector_store(r34_1);      
    
    up1++;
    _sse_float_prefetch_su3(up1);      
    
    _sse_float_pair_load((*(sm2)).c1,(*(sm2)).c2);
    _sse_float_pair_load_up((*(sm2)).c4,(*(sm2)).c3);
    _sse_float_vector_add();
    
    _sse_float_su3_inverse_multiply((*(um2)));
    
    _sse_float_vector_load(r12_2);
    _sse_float_vector_add();
    _sse_float_vector_store(r12_2);      
    
    _sse_float_vector_load(r34_2);
    _sse_float_vector_xch();
    _sse_float_vector_add();
    _sse_float_vector_store(r34_2); 
    
    /******************************* direction +2 ***********************/
    
    iy1=idw[ix1][2];
    iy2=idw[ix1+1][2];
    sm1=psi+iy1;
    sm2=psi+iy2;
    _sse_float_prefetch_spinor(sm1);
    _sse_float_prefetch_spinor(sm2);

    _sse_float_pair_load((*sp1).c1,(*sp1).c2);
    _sse_float_pair_load_up((*sp1).c4,(*sp1).c3);
    _sse_float_vector_i_addsub();
    
    _sse_float_su3_multiply((*up1));
    
    _sse_float_vector_load(r12_1);
    _sse_float_vector_add();
    _sse_float_vector_store(r12_1);       
    
    _sse_float_vector_load(r34_1);
    _sse_float_vector_xch();
    _sse_float_vector_i_addsub();
    _sse_float_vector_store(r34_1);       
    

    um1=U+iy1*4+2;
    um2=U+iy2*4+2;     
    _sse_float_prefetch_su3(um1);
    _sse_float_prefetch_su3(um2);
    
    _sse_float_pair_load((*(sp2)).c1,(*(sp2)).c2);
    _sse_float_pair_load_up((*(sp2)).c4,(*(sp2)).c3);
    _sse_float_vector_i_addsub();
    
    _sse_float_su3_multiply((*(up1+4)));
    
    _sse_float_vector_load(r12_2);
    _sse_float_vector_add();
    _sse_float_vector_store(r12_2);       
    
    _sse_float_vector_load(r34_2);
    _sse_float_vector_xch();
    _sse_float_vector_i_addsub();
    _sse_float_vector_store(r34_2); 
    
    /******************************* direction -2 ***********************/
    
    sp1=psi+iup[ix1][3];
    sp2=psi+iup[ix1+1][3];
    _sse_float_prefetch_spinor(sp1);
    _sse_float_prefetch_spinor(sp2);

    _sse_float_pair_load((*sm1).c1,(*sm1).c2);
    _sse_float_pair_load_up((*sm1).c4,(*sm1).c3);
    _sse_float_vector_i_subadd();      
    
    _sse_float_su3_inverse_multiply((*um1));
    
    _sse_float_vector_load(r12_1);
    _sse_float_vector_add();
    _sse_float_vector_store(r12_1);
    
    _sse_float_vector_load(r34_1);
    _sse_float_vector_xch();
    _sse_float_vector_i_subadd();
    _sse_float_vector_store(r34_1);
    
    up1++;
    _sse_float_prefetch_su3(up1);
    
    _sse_float_pair_load((*(sm2)).c1,(*(sm2)).c2);
    _sse_float_pair_load_up((*(sm2)).c4,(*(sm2)).c3);
    _sse_float_vector_i_subadd();      
    
    _sse_float_su3_inverse_multiply((*(um2)));
    
    _sse_float_vector_load(r12_2);
    _sse_float_vector_add();
    _sse_float_vector_store(r12_2);
    
    _sse_float_vector_load(r34_2);
    _sse_float_vector_xch();
    _sse_float_vector_i_subadd();
    _sse_float_vector_store(r34_2);
    
    /******************************* direction +3 ***********************/

    iy1=idw[ix1][3];
    iy2=idw[ix1+1][3];
    sm1=psi+iy1;
    sm2=psi+iy2;
    _sse_float_prefetch_spinor(sm1);
    _sse_float_prefetch_spinor(sm2);

    _sse_float_pair_load((*sp1).c1,(*sp1).c2);
    _sse_float_pair_load_up((*sp1).c3,(*sp1).c4);
    _sse_float_vector_subadd();
    
    _sse_float_su3_multiply((*up1));
      
    _sse_float_vector_load(r12_1);
    _sse_float_vector_add();
    _sse_float_vector_store(r12_1);
    
    _sse_float_vector_load(r34_1);
    _sse_float_vector_subadd();
    _sse_float_vector_store(r34_1);      
    
    um1=U+iy1*4+3;
    _sse_float_prefetch_su3(um1);
    um2=U+iy2*4+3;     
    _sse_float_prefetch_su3(um2);

    _sse_float_pair_load((*sp2).c1,(*sp2).c2);
    _sse_float_pair_load_up((*sp2).c3,(*sp2).c4);
    _sse_float_vector_subadd();
    
    _sse_float_su3_multiply((*(up1+4)));
    
    _sse_float_vector_load(r12_2);
    _sse_float_vector_add();
    _sse_float_vector_store(r12_2);
    
    _sse_float_vector_load(r34_2);
    _sse_float_vector_subadd();
    _sse_float_vector_store(r34_2); 
    
/******************************* direction -3 ***********************/

    sn1=(_sse_spinor*) &chi[ix1];      
    _sse_float_prefetch_spinor(sn1);
    
    iz1=ix1+2;
    if (iz1<stop) {
      sp1=(_sse_spinor*) &psi[iup[iz1][0]];
      sp2=(_sse_spinor*) &psi[iup[iz1+1][0]];
      _sse_float_prefetch_spinor(sp1);
      _sse_float_prefetch_spinor(sp2);  
    }

    _sse_float_pair_load((*sm1).c1,(*sm1).c2);
    _sse_float_pair_load_up((*sm1).c3,(*sm1).c4);
    _sse_float_vector_addsub();
    
    _sse_float_su3_inverse_multiply((*um1));
        
    _sse_float_vector_load(r12_1);
    _sse_float_vector_add();
    _sse_float_vector_mul(fact2);
    _sse_float_pair_store((*sn1).c1,(*sn1).c2);
    
    _sse_float_vector_load(r34_1);
    _sse_float_vector_addsub();
    _sse_float_vector_mul(fact2);
    _sse_float_pair_store((*sn1).c3,(*sn1).c4);      
    
    up1=U+iz1*4;
    _sse_float_prefetch_su3(up1);

    
    _sse_float_pair_load((*sm2).c1,(*sm2).c2);
    _sse_float_pair_load_up((*sm2).c3,(*sm2).c4);
    _sse_float_vector_addsub();
    
    _sse_float_su3_inverse_multiply((*um2));
    
    _sse_float_vector_load(r12_2);
    _sse_float_vector_add();
    _sse_float_vector_mul(fact2);
    _sse_float_pair_store((*(sn1+1)).c1,(*(sn1+1)).c2);
    
    _sse_float_vector_load(r34_2);
    _sse_float_vector_addsub();
    _sse_float_vector_mul(fact2);
    _sse_float_pair_store((*(sn1+1)).c3,(*(sn1+1)).c4); 
    
    /******************************** end of loop ***********************/
  }

  if(cSW==0) return;

  /* 
     everything here must be in agreement with gauge_field::ordered_index()
  */
    
  /*********** loop over all lattice sites for clover term *************/
    
  um1=uem+6*start;
  
  for (ix1=start; ix1<stop; ix1+=2) {   
    s1=psi+ix1;
    _sse_float_prefetch_spinor(s1);

    //************************* mu=0, nu=1 ***********************
      
      _sse_float_prefetch_su3(um1+1);
      _sse_float_prefetch_su3(um1+7);
      
      _sse_float_pair_load((*s1).c4,(*s1).c3);
      _sse_float_vector_mul(fact4);
      _sse_float_su3_multiply((*um1));
      // set this to zero 
      _sse_float_vector_load(r0);
      _sse_float_vector_add();
      _sse_float_vector_store(r12_1);

      _sse_float_pair_load((*s1).c2,(*s1).c1);
      _sse_float_su3_multiply((*um1));    
      _sse_float_vector_load(r0);
      _sse_float_vector_add();
      _sse_float_vector_store(r34_1);
      
      _sse_float_pair_load((*(s1+1)).c4,(*(s1+1)).c3);
      _sse_float_vector_mul(fact4);
      _sse_float_su3_multiply((*(um1+6)));    
      // set this to zero 
      _sse_float_vector_load(r0);
      _sse_float_vector_add();
      _sse_float_vector_store(r12_2);

      _sse_float_pair_load((*(s1+1)).c2,(*(s1+1)).c1);
      _sse_float_su3_multiply((*(um1+6)));    
      // set this to zero 
      _sse_float_vector_load(r0);
      _sse_float_vector_add();
      _sse_float_vector_store(r34_2);
      
      um1++;
      
      //************************* mu=0, nu=2 ***********************
      
      _sse_float_prefetch_su3(um1+1);
      _sse_float_prefetch_su3(um1+7);
      
      _sse_float_pair_load((*s1).c4,(*s1).c3);
      _sse_float_su3_multiply((*um1));
      _sse_float_vector_load(r12_1);
      _sse_float_vector_i_addsub();
      _sse_float_vector_store(r12_1);
      
      _sse_float_pair_load((*s1).c2,(*s1).c1);
      _sse_float_su3_multiply((*um1));
      _sse_float_vector_load(r34_1);
      _sse_float_vector_i_subadd();
      _sse_float_vector_store(r34_1);
      
      _sse_float_pair_load((*(s1+1)).c4,(*(s1+1)).c3);
      _sse_float_su3_multiply((*(um1+6)));
      _sse_float_vector_load(r12_2);
      _sse_float_vector_i_addsub();
      _sse_float_vector_store(r12_2);
      
      _sse_float_pair_load((*(s1+1)).c2,(*(s1+1)).c1);
      _sse_float_su3_multiply((*(um1+6)));
      _sse_float_vector_load(r34_2);
      _sse_float_vector_i_subadd();
      _sse_float_vector_store(r34_2);
      
      um1++;
      
      //************************* mu=0, nu=3 ***********************
      
      _sse_float_prefetch_su3(um1+1);
      _sse_float_prefetch_su3(um1+7);
      
      _sse_float_pair_load((*s1).c3,(*s1).c4);
      _sse_float_su3_multiply((*um1));
      _sse_float_vector_load(r12_1);
      _sse_float_vector_subadd();
      _sse_float_vector_mul(fact6);
      _sse_float_vector_store(r12_1);
      
      _sse_float_pair_load((*s1).c1,(*s1).c2);
      _sse_float_su3_multiply((*um1));
      _sse_float_vector_load(r34_1);
      _sse_float_vector_addsub();
      _sse_float_vector_mul(fact6);
      _sse_float_vector_store(r34_1);
      
      _sse_float_pair_load((*(s1+1)).c3,(*(s1+1)).c4);
      _sse_float_su3_multiply((*(um1+6)));
      _sse_float_vector_load(r12_2);
      _sse_float_vector_subadd();
      _sse_float_vector_mul(fact6);
      _sse_float_vector_store(r12_2);
      
      _sse_float_pair_load((*(s1+1)).c1,(*(s1+1)).c2);
      _sse_float_su3_multiply((*(um1+6)));
      _sse_float_vector_load(r34_2);
      _sse_float_vector_addsub();
      _sse_float_vector_mul(fact6);
      _sse_float_vector_store(r34_2);
      
      um1++;
      
      //************************* mu=1, nu=2 ***********************
      
      _sse_float_prefetch_su3(um1+1);
      _sse_float_prefetch_su3(um1+7);
      
      _sse_float_pair_load((*s1).c1,(*s1).c2);
      _sse_float_su3_multiply((*um1));
      _sse_float_vector_load(r12_1);
      _sse_float_vector_i_subadd();
      _sse_float_vector_store(r12_1);
      
      _sse_float_pair_load((*s1).c3,(*s1).c4);
      _sse_float_su3_multiply((*um1));
      _sse_float_vector_load(r34_1);
      _sse_float_vector_i_subadd();
      _sse_float_vector_store(r34_1);
      
      _sse_float_pair_load((*(s1+1)).c1,(*(s1+1)).c2);
      _sse_float_su3_multiply((*(um1+6)));
      _sse_float_vector_load(r12_2);
      _sse_float_vector_i_subadd();
      _sse_float_vector_store(r12_2);
      
      _sse_float_pair_load((*(s1+1)).c3,(*(s1+1)).c4);
      _sse_float_su3_multiply((*(um1+6)));
      _sse_float_vector_load(r34_2);
      _sse_float_vector_i_subadd();
      _sse_float_vector_store(r34_2);
      
      um1++;
      
      //************************* mu=1, nu=3 ***********************
      
      _sse_float_prefetch_su3(um1+1);
      _sse_float_prefetch_su3(um1+7);
      
      _sse_float_pair_load((*s1).c2,(*s1).c1);
      _sse_float_su3_multiply((*um1));
      _sse_float_vector_load(r12_1);
      _sse_float_vector_addsub();
      _sse_float_vector_store(r12_1);
      
      _sse_float_pair_load((*s1).c4,(*s1).c3);
      _sse_float_su3_multiply((*um1));
      _sse_float_vector_load(r34_1);
      _sse_float_vector_addsub();
      _sse_float_vector_store(r34_1);
      
      _sse_float_pair_load((*(s1+1)).c2,(*(s1+1)).c1);
      _sse_float_su3_multiply((*(um1+6)));
      _sse_float_vector_load(r12_2);
      _sse_float_vector_addsub();
      _sse_float_vector_store(r12_2);
      
      _sse_float_pair_load((*(s1+1)).c4,(*(s1+1)).c3);
      _sse_float_su3_multiply((*(um1+6)));
      _sse_float_vector_load(r34_2);
      _sse_float_vector_addsub();
      _sse_float_vector_store(r34_2);
      
      um1++;
      
      //************************* mu=2, nu=3 ***********************
      
      sn1=(_sse_spinor*) &chi[ix1];      
      _sse_float_prefetch_spinor(sn1);
      
      _sse_float_pair_load((*s1).c2,(*s1).c1);
      _sse_float_su3_multiply((*um1));
      _sse_float_vector_load(r12_1);
      _sse_float_vector_i_sub();
      _sse_float_vector_store(r12_1);
      
      _sse_float_pair_load((*s1).c4,(*s1).c3);
      _sse_float_su3_multiply((*um1));
      _sse_float_vector_load(r34_1);
      _sse_float_vector_i_sub();
      _sse_float_vector_store(r34_1);
      
      _sse_float_pair_load((*(s1+1)).c2,(*(s1+1)).c1);
      _sse_float_su3_multiply((*(um1+6)));
      _sse_float_vector_load(r12_2);
      _sse_float_vector_i_sub();
      _sse_float_vector_store(r12_2);
      
      _sse_float_pair_load((*(s1+1)).c4,(*(s1+1)).c3);
      _sse_float_su3_multiply((*(um1+6)));
      _sse_float_vector_load(r34_2);
      _sse_float_vector_i_sub();
      _sse_float_vector_store(r34_2);
      
      um1+=7;
      if(ix1<stop-1) {
	_sse_float_prefetch_su3(um1);
	_sse_float_prefetch_su3(um1+6);
      }


      _sse_float_pair_load_up((*sn1).c1,(*sn1).c2);
      _sse_float_vector_load(r12_1);
      _sse_float_vector_mul(fact5);
      _sse_float_vector_add();
      _sse_float_pair_store((*sn1).c1,(*sn1).c2);

      _sse_float_pair_load_up((*sn1).c3,(*sn1).c4);
      _sse_float_vector_load(r34_1);
      _sse_float_vector_mul(fact5);
      _sse_float_vector_add();
      _sse_float_pair_store((*sn1).c3,(*sn1).c4);

      _sse_float_pair_load_up((*(sn1+1)).c1,(*(sn1+1)).c2);
      _sse_float_vector_load(r12_2);
      _sse_float_vector_mul(fact5);
      _sse_float_vector_add();
      _sse_float_pair_store((*(sn1+1)).c1,(*(sn1+1)).c2);

      _sse_float_pair_load_up((*(sn1+1)).c3,(*(sn1+1)).c4);
      _sse_float_vector_load(r34_2);
      _sse_float_vector_mul(fact5);
      _sse_float_vector_add();
      _sse_float_pair_store((*(sn1+1)).c3,(*(sn1+1)).c4);


      //************** end of loop ***********************
  }


#else

  _sse_spinor *chi=(_sse_spinor*) chi_out.physical_address();
  _sse_spinor *psi=(_sse_spinor*) psi_in.physical_address();
  _sse_su3    *U=(_sse_su3*) U_in.physical_address();
  _sse_su3    *uem=(_sse_su3*) U_in.em.physical_address();
  mdp_int   **iup=U_in.lattice().up;
  mdp_int   **idw=U_in.lattice().dw;
  mdp_int   start=U_in.lattice().start_index(ME,0); // even
  mdp_int   stop =U_in.lattice().stop_index(ME,1);  // odd
  
  _sse_double fact1 ALIGN16;
  _sse_double fact2 ALIGN16;
  _sse_double fact3 ALIGN16;
  _sse_double fact4 ALIGN16;
  _sse_double fact5 ALIGN16;
  _sse_double fact6 ALIGN16;
  _sse_spinor rs ALIGN16;
  _sse_spinor r0 ALIGN16;
  mdp_int ix,iy,iz;
  double rho;
  _sse_su3 *up, *um;   
  _sse_spinor *s,*sp,*sm,*sn;

  if(sign!=1) exit(1);
  if((stop-start)%2 !=0)
    error("FermiCloverActionSSE2\nProblem with parallelization: odd # of sites on process!");

  if(r_t!=1.0) 
    error("FermiCloverActionSSE2\nr_t!=1 not compatible with SSE2\n");

  _sse_check_alignment((void*) &fact1, 0xf);
  _sse_check_alignment((void*) &fact2, 0xf);
  _sse_check_alignment((void*) &fact3, 0xf);
  _sse_check_alignment((void*) &fact4, 0xf);
  _sse_check_alignment((void*) &fact5, 0xf);
  _sse_check_alignment((void*) &fact6, 0xf);
  _sse_check_alignment((void*) &rs, 0xf);
  _sse_check_alignment((void*) &r0, 0xf);

  r0.c1.c1.real()=r0.c1.c2.real()=r0.c1.c3.real()=0;
  r0.c2.c1.real()=r0.c2.c2.real()=r0.c2.c3.real()=0;
  r0.c3.c1.real()=r0.c3.c2.real()=r0.c3.c3.real()=0;
  r0.c4.c1.real()=r0.c4.c2.real()=r0.c4.c3.real()=0;

  rho=-1.0/kappa_s;

  fact1.c1=rho;
  fact1.c2=rho;
  
  fact2.c1=-1.0*kappa_s;
  fact2.c2=fact2.c1;

  fact3.c1=(1.0+r_t)*kappa_t/kappa_s;
  fact3.c2=fact3.c1;

  fact4.c1=-1.0;
  fact4.c2=-1.0;

  fact5.c1=1.0*kappa_s*cSW*c_E;
  fact5.c2=fact5.c1;

  fact6.c1=1.0*c_E/c_B;
  fact6.c2=fact6.c1;

  sp=(_sse_spinor*) &psi[iup[start][0]];
  up=(_sse_su3*) U+4*start;

   /************************ loop over all lattice sites ***************/

  for (ix=start; ix<stop; ix++) {   
    s=psi+ix;
    _sse_double_prefetch_spinor(s);
    
    /******************************* direction +0 ***********************/

    iy=idw[ix][0];
    sm=psi+iy;
    _sse_double_prefetch_spinor(sm);

    _sse_double_load((*s).c1);
    _sse_double_vector_mul(fact1);      
    _sse_double_store(rs.c1);
    _sse_double_load((*s).c2);
    _sse_double_vector_mul(fact1);      
    _sse_double_store(rs.c2);

    um=U+iy*4;
    _sse_double_prefetch_su3(um);

    _sse_double_load((*sp).c3);
    _sse_double_vector_mul(fact3);
    _sse_double_su3_multiply((*up));
    _sse_double_load((*s).c3);
    _sse_double_vector_mul(fact1);      
    _sse_double_vector_add();
    _sse_double_store(rs.c3);

    _sse_double_load((*sp).c4);
    _sse_double_vector_mul(fact3);
    _sse_double_su3_multiply((*up));
    _sse_double_load((*s).c4);
    _sse_double_vector_mul(fact1);      
    _sse_double_vector_add();
    _sse_double_store(rs.c4);
   
    
    /******************************* direction -0 ***********************/

    sp=psi+iup[ix][1];
    _sse_double_prefetch_spinor(sp);
    up++;
    _sse_double_prefetch_su3(up);      
        
    _sse_double_load((*sm).c1);
    _sse_double_vector_mul(fact3);          
    _sse_double_su3_inverse_multiply((*um));
    _sse_double_load(rs.c1);
    _sse_double_vector_add();
    _sse_double_store(rs.c1);

    _sse_double_load((*sm).c2);
    _sse_double_vector_mul(fact3);          
    _sse_double_su3_inverse_multiply((*um));
    _sse_double_load(rs.c2);
    _sse_double_vector_add();
    _sse_double_store(rs.c2);
      
    /******************************* direction +1 ***********************/
        
    iy=idw[ix][1];
    sm=psi+iy;
    _sse_double_prefetch_spinor(sm);
    um=U+iy*4+1;
    _sse_double_prefetch_su3(um);      

    _sse_double_load((*sp).c1);
    _sse_double_load_up((*sp).c4);
    _sse_double_vector_sub();
    _sse_double_su3_multiply((*up));
    _sse_double_load(rs.c1);
    _sse_double_vector_add();
    _sse_double_store(rs.c1);
    _sse_double_load(rs.c4);
    _sse_double_vector_sub();
    _sse_double_store(rs.c4);

    _sse_double_load((*sp).c2);
    _sse_double_load_up((*sp).c3);
    _sse_double_vector_sub();
    _sse_double_su3_multiply((*up));
    _sse_double_load(rs.c2);
    _sse_double_vector_add();
    _sse_double_store(rs.c2);
    _sse_double_load(rs.c3);
    _sse_double_vector_sub();
    _sse_double_store(rs.c3);
    
    
    /******************************* direction -1 ***********************/
    
    sp=psi+iup[ix][2];
    _sse_double_prefetch_spinor(sp);
    up++;
    _sse_double_prefetch_su3(up);      

    _sse_double_load((*sm).c1);
    _sse_double_load_up((*sm).c4);
    _sse_double_vector_add();
    _sse_double_su3_inverse_multiply((*um));
    _sse_double_load(rs.c1);
    _sse_double_vector_add();
    _sse_double_store(rs.c1);      
    _sse_double_load(rs.c4);
    _sse_double_vector_add();
    _sse_double_store(rs.c4);      
  
    _sse_double_load((*sm).c2);
    _sse_double_load_up((*sm).c3);
    _sse_double_vector_add();
    _sse_double_su3_inverse_multiply((*um));
    _sse_double_load(rs.c2);
    _sse_double_vector_add();
    _sse_double_store(rs.c2);      
    _sse_double_load(rs.c3);
    _sse_double_vector_add();
    _sse_double_store(rs.c3);      

    /******************************* direction +2 ***********************/
    
    iy=idw[ix][2];
    sm=psi+iy;
    _sse_double_prefetch_spinor(sm);
    um=U+iy*4+2;
    _sse_double_prefetch_su3(um);

    _sse_double_load((*sp).c1);
    _sse_double_load_up((*sp).c4);
    _sse_double_vector_i_mul(); _sse_double_vector_add();
    _sse_double_su3_multiply((*up));
    _sse_double_load(rs.c1);
    _sse_double_vector_add();
    _sse_double_store(rs.c1);       
    _sse_double_load(rs.c4);
    _sse_double_vector_i_mul(); _sse_double_vector_sub();
    _sse_double_store(rs.c4);       
    
    _sse_double_load((*sp).c2);
    _sse_double_load_up((*sp).c3);
    _sse_double_vector_i_mul(); _sse_double_vector_sub();
    _sse_double_su3_multiply((*up));
    _sse_double_load(rs.c2);
    _sse_double_vector_add();
    _sse_double_store(rs.c2);       
    _sse_double_load(rs.c3);
    _sse_double_vector_i_mul(); _sse_double_vector_add();
    _sse_double_store(rs.c3); 

    
    /******************************* direction -2 ***********************/
    
    sp=psi+iup[ix][3];
    _sse_double_prefetch_spinor(sp);
    up++;
    _sse_double_prefetch_su3(up);

    _sse_double_load((*sm).c1);
    _sse_double_load_up((*sm).c4);
    _sse_double_vector_i_mul(); _sse_double_vector_sub();      
    _sse_double_su3_inverse_multiply((*um));
    _sse_double_load(rs.c1);
    _sse_double_vector_add();
    _sse_double_store(rs.c1);
    _sse_double_load(rs.c4);
    _sse_double_vector_i_mul(); _sse_double_vector_add();
    _sse_double_store(rs.c4);
    
    _sse_double_load((*sm).c2);
    _sse_double_load_up((*sm).c3);
    _sse_double_vector_i_mul(); _sse_double_vector_add();      
    _sse_double_su3_inverse_multiply((*um));
    _sse_double_load(rs.c2);
    _sse_double_vector_add();
    _sse_double_store(rs.c2);
    _sse_double_load(rs.c3);
    _sse_double_vector_i_mul(); _sse_double_vector_sub();
    _sse_double_store(rs.c3);
    
    /******************************* direction +3 ***********************/

    iy=idw[ix][3];
    sm=psi+iy;
    _sse_double_prefetch_spinor(sm);
    um=U+iy*4+3;
    _sse_double_prefetch_su3(um);

    _sse_double_load((*sp).c1);
    _sse_double_load_up((*sp).c3);
    _sse_double_vector_sub();
    _sse_double_su3_multiply((*up));
    _sse_double_load(rs.c1);
    _sse_double_vector_add();
    _sse_double_store(rs.c1);
    _sse_double_load(rs.c3);
    _sse_double_vector_sub();
    _sse_double_store(rs.c3);      
    
    _sse_double_load((*sp).c2);
    _sse_double_load_up((*sp).c4);
    _sse_double_vector_add();
    _sse_double_su3_multiply((*up));
    _sse_double_load(rs.c2);
    _sse_double_vector_add();
    _sse_double_store(rs.c2);
    _sse_double_load(rs.c4);
    _sse_double_vector_add();
    _sse_double_store(rs.c4);
    
    /******************************* direction -3 ***********************/

    sn=chi+ix;      
    _sse_double_prefetch_spinor(sn);
    
    iz=ix+1;
    if (iz<stop) {
      sp=psi+iup[iz][0];
      _sse_double_prefetch_spinor(sp);
      up=U+iz*4;
      _sse_double_prefetch_su3(up);
    }

    _sse_double_load((*sm).c1);
    _sse_double_load_up((*sm).c3);
    _sse_double_vector_add();
    _sse_double_su3_inverse_multiply((*um));
    _sse_double_load(rs.c1);
    _sse_double_vector_add();
    _sse_double_vector_mul(fact2);
    _sse_double_store((*sn).c1);
    _sse_double_load(rs.c3);
    _sse_double_vector_add();
    _sse_double_vector_mul(fact2);
    _sse_double_store((*sn).c3);
    
    _sse_double_load((*sm).c2);
    _sse_double_load_up((*sm).c4);
    _sse_double_vector_sub(); 
    _sse_double_su3_inverse_multiply((*um));
    _sse_double_load(rs.c2);
    _sse_double_vector_add();
    _sse_double_vector_mul(fact2);
    _sse_double_store((*sn).c2);
    _sse_double_load(rs.c4);
    _sse_double_vector_sub();
    _sse_double_vector_mul(fact2);
    _sse_double_store((*sn).c4);      
    
    /******************************** end of loop ***********************/    
  }

  if(cSW==0) return;

  /* 
     everything here must be in agreement with gauge_field::ordered_index()
  */
    
  /*********** loop over all lattice sites for clover term *************/
    
  um=uem+6*start;

  for (ix=start; ix<stop; ix++) {   
    s=psi+ix;
    _sse_double_prefetch_spinor(s);
    
    /************************** mu=0, nu=1 ***********************/

      
      _sse_double_prefetch_su3(um+1);
      
      _sse_double_load((*s).c4);
      _sse_double_vector_mul(fact4);
      _sse_double_su3_multiply((*um));
      /* set this to zero */
      _sse_double_load(r0.c1);
      _sse_double_vector_add();
      _sse_double_store(rs.c1);

      _sse_double_load((*s).c3);
      _sse_double_vector_mul(fact4);
      _sse_double_su3_multiply((*um));
      /* set this to zero */
      _sse_double_load(r0.c2);
      _sse_double_vector_add();
      _sse_double_store(rs.c2);

      _sse_double_load((*s).c2);
      _sse_double_su3_multiply((*um));    
      _sse_double_load(r0.c3);
      _sse_double_vector_add();
      _sse_double_store(rs.c3);
      
      _sse_double_load((*s).c1);
      _sse_double_su3_multiply((*um));    
      _sse_double_load(r0.c4);
      _sse_double_vector_add();
      _sse_double_store(rs.c4);
    
      um++;
      
      /************************** mu=0, nu=2 ***********************/
      
      _sse_double_prefetch_su3(um+1);
      
      _sse_double_load((*s).c4);
      _sse_double_su3_multiply((*um));
      _sse_double_load(rs.c1);
      _sse_double_vector_i_mul(); _sse_double_vector_add();
      _sse_double_store(rs.c1);
      
      _sse_double_load((*s).c3);
      _sse_double_su3_multiply((*um));
      _sse_double_load(rs.c2);
      _sse_double_vector_i_mul(); _sse_double_vector_sub();
      _sse_double_store(rs.c2);
      
      _sse_double_load((*s).c2);
      _sse_double_su3_multiply((*um));
      _sse_double_load(rs.c3);
      _sse_double_vector_i_mul(); _sse_double_vector_sub();
      _sse_double_store(rs.c3);
      
      _sse_double_load((*s).c1);
      _sse_double_su3_multiply((*um));
      _sse_double_load(rs.c4);
      _sse_double_vector_i_mul(); _sse_double_vector_add();
      _sse_double_store(rs.c4); 
      
      um++;
      
      /************************** mu=0, nu=3 ***********************/
      
      _sse_double_prefetch_su3(um+1);
      
      _sse_double_load((*s).c3);
      _sse_double_su3_multiply((*um));
      _sse_double_load(rs.c1);
      _sse_double_vector_sub();
      _sse_double_vector_mul(fact6);
      _sse_double_store(rs.c1);
      
      _sse_double_load((*s).c4);
      _sse_double_su3_multiply((*um));
      _sse_double_load(rs.c2);
      _sse_double_vector_add();
      _sse_double_vector_mul(fact6);
      _sse_double_store(rs.c2);
      
      _sse_double_load((*s).c1);
      _sse_double_su3_multiply((*um));
      _sse_double_load(rs.c3);
      _sse_double_vector_add();
      _sse_double_vector_mul(fact6);
      _sse_double_store(rs.c3);
          
      _sse_double_load((*s).c2);
      _sse_double_su3_multiply((*um));
      _sse_double_load(rs.c4);
      _sse_double_vector_sub();
      _sse_double_vector_mul(fact6);
      _sse_double_store(rs.c4);
          
      um++;
      
      /************************** mu=1, nu=2 ***********************/
      
      _sse_double_prefetch_su3(um+1);
      
      _sse_double_load((*s).c1);
      _sse_double_su3_multiply((*um));
      _sse_double_load(rs.c1);
      _sse_double_vector_i_mul(); _sse_double_vector_sub();
      _sse_double_store(rs.c1);
      
      _sse_double_load((*s).c2);
      _sse_double_su3_multiply((*um));
      _sse_double_load(rs.c2);
      _sse_double_vector_i_mul(); _sse_double_vector_add();
      _sse_double_store(rs.c2);
      
      _sse_double_load((*s).c3);
      _sse_double_su3_multiply((*um));
      _sse_double_load(rs.c3);
      _sse_double_vector_i_mul(); _sse_double_vector_sub();
      _sse_double_store(rs.c3);

      _sse_double_load((*s).c4);
      _sse_double_su3_multiply((*um));
      _sse_double_load(rs.c4);
      _sse_double_vector_i_mul(); _sse_double_vector_add();
      _sse_double_store(rs.c4);
      
      um++;
      
      /************************** mu=1, nu=3 ***********************/
      
      _sse_double_prefetch_su3(um+1);
      
      _sse_double_load((*s).c2);
      _sse_double_su3_multiply((*um));
      _sse_double_load(rs.c1);
      _sse_double_vector_add();
      _sse_double_store(rs.c1);
      
      _sse_double_load((*s).c1);
      _sse_double_su3_multiply((*um));
      _sse_double_load(rs.c2);
      _sse_double_vector_sub();
      _sse_double_store(rs.c2);
      
      _sse_double_load((*s).c4);
      _sse_double_su3_multiply((*um));
      _sse_double_load(rs.c3);
      _sse_double_vector_add();
      _sse_double_store(rs.c3);
      
      _sse_double_load((*s).c3);
      _sse_double_su3_multiply((*um));
      _sse_double_load(rs.c4);
      _sse_double_vector_sub();
      _sse_double_store(rs.c4);
      
      um++;
      
      /************************** mu=2, nu=3 ***********************/

      sn=(_sse_spinor*) &chi[ix];      
      _sse_double_prefetch_spinor(sn);
      
      _sse_double_load((*s).c2);
      _sse_double_su3_multiply((*um));
      _sse_double_load(rs.c1);
      _sse_double_vector_i_mul(); _sse_double_vector_sub();
      _sse_double_store(rs.c1);
      
      _sse_double_load((*s).c1);
      _sse_double_su3_multiply((*um));
      _sse_double_load(rs.c2);
      _sse_double_vector_i_mul(); _sse_double_vector_sub();
      _sse_double_store(rs.c2);
      
      _sse_double_load((*s).c4);
      _sse_double_su3_multiply((*um));
      _sse_double_load(rs.c3);
      _sse_double_vector_i_mul(); _sse_double_vector_sub();
      _sse_double_store(rs.c3);
      
      _sse_double_load((*s).c3);
      _sse_double_su3_multiply((*um));
      _sse_double_load(rs.c4);
      _sse_double_vector_i_mul(); _sse_double_vector_sub();
      _sse_double_store(rs.c4);

      um++;
      if(ix<stop) {
	_sse_double_prefetch_su3(um);
      }

      _sse_double_load_up((*sn).c1);
      _sse_double_load(rs.c1);
      _sse_double_vector_mul(fact5);
      _sse_double_vector_add();
      _sse_double_store((*sn).c1);

      _sse_double_load_up((*sn).c2);
      _sse_double_load(rs.c2);
      _sse_double_vector_mul(fact5);
      _sse_double_vector_add();
      _sse_double_store((*sn).c2);

      _sse_double_load_up((*sn).c3);
      _sse_double_load(rs.c3);
      _sse_double_vector_mul(fact5);
      _sse_double_vector_add();
      _sse_double_store((*sn).c3);

      _sse_double_load_up((*sn).c4);
      _sse_double_load(rs.c4);
      _sse_double_vector_mul(fact5);
      _sse_double_vector_add();
      _sse_double_store((*sn).c4);
      
      /*************** end of loop ***********************/
  }
#endif // if defined(USE_DOUBLE_PRECISION)
 
  }
};

#endif // if defined(SSE2)


