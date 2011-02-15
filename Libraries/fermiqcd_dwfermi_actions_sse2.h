/////////////////////////////////////////////////////////////////
/// @file fermiqcd_fermi_actions_sse2.h
/// @version 2009-12-21
/// @author Martin Luescher and Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Basic actions for Wilson Fermions optimized in assembler
///
//////////////////////////////////////////////////////////////////

/// WARNING!!!! THIS WAS NOT TESTED !!!!

#if defined(SSE2)

/// @brief Domain Wilson action SSE/SSE2
///
/// Only on Pentium 4 or compatible. Compile with -DSSE2
/// Example:
/// @verbatim
/// gauge_field U(lattice,nc);
/// dwfermi_field psi(lattice,nc);
/// dwfermi_field chi(lattice,nc);
/// coefficients coeff;
/// coeff["m_5"]=0.11;
/// coeff["m_f"]=0.11;
/// default_dwfermi_action=DWFermiActionSSE2::mul_Q;
/// mul_Q(chi,psi,U,coeff);
/// @endverbatim
/// Attention: uses always FERMILAB convention of Gamma matrices
class DwFermiActionSSE2 {
 public:
  static void mul_Q(dwfermi_field &chi_out, 
		    dwfermi_field &psi_in, 
		    gauge_field &U_in, 
		    coefficients &coeff) { 
    
    if(psi_in.nspin!=4) error("fermiqcd_dwfermi_algorithms/dwfermi_mul_Q_ONE: nspin!=4");
    if(psi_in.nc!=U_in.nc) error("fermiqcd_dwfermi_algorithms/dwfermi_mul_Q_ONE: gauge and spinor have different nc");
    
    register int    ndim=psi_in.lattice().ndim;
    register int    nspin=psi_in.nspin;
    register int    nc=psi_in.nc;
    register int    L5=psi_in.L5;
    register mdp_real m_5,m_f,sign;
    if(coeff.has_key("m_5")) m_5=coeff["m_5"];
    else error("coefficients m_5 undeclared");
    if(coeff.has_key("m_f")) m_f=coeff["m_f"];
    else error("coefficients m_f undeclared");
    if(coeff.has_key("sign")) sign=coeff["sign"];
    else sign=1;
    
    // check the sign and the 6.0 here
    register mdp_real kappa5=0.5/(m_5-6.0);
    register mdp_real kappaf=-m_f*kappa5;
    
#if !defined(USE_DOUBLE_PRECISION)

    _sse_spinor *chi=(_sse_spinor*) chi_out.physical_address();
    _sse_check_alignment((void*) chi, 0xf);
    _sse_spinor *psi=(_sse_spinor*) psi_in.physical_address();
    _sse_check_alignment((void*) psi, 0xf);
    _sse_su3    *U=(_sse_su3*) U_in.physical_address();
    _sse_check_alignment((void*) U, 0xf);
    mdp_int   **iup=U_in.lattice().up;
    mdp_int   **idw=U_in.lattice().dw;
    mdp_int   start=U_in.lattice().start_index(ME,0); // even
    mdp_int   stop =U_in.lattice().stop_index(ME,1);  // odd
    
    _sse_float fact1 ALIGN16;
    _sse_float fact2 ALIGN16;
    _sse_float fact3 ALIGN16;
    _sse_float fact4 ALIGN16;
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
      error("DWFermiActionSSE2\nProblem with parallelization: odd # of sites on process!");
    
    _sse_check_alignment((void*) &fact1, 0xf);
    _sse_check_alignment((void*) &fact2, 0xf);
    _sse_check_alignment((void*) &fact3, 0xf);
    _sse_check_alignment((void*) &fact4, 0xf);
    _sse_check_alignment((void*) &r12_1, 0xf);
    _sse_check_alignment((void*) &r34_1, 0xf);
    _sse_check_alignment((void*) &r12_2, 0xf);
    _sse_check_alignment((void*) &r34_2, 0xf);
    
    /* empty vector */
    
    r0.c1.c1=r0.c1.c2=r0.c1.c3=r0.c1.c4=0;
    r0.c2.c1=r0.c2.c2=r0.c2.c3=r0.c2.c4=0;
    r0.c3.c1=r0.c3.c2=r0.c3.c3=r0.c3.c4=0;
    
    rho=-1.0f/kappa5;
    
    /* coefficient of (1 +- Gamma[i]) */
    
    fact1.c1=rho;
    fact1.c2=rho;
    fact1.c3=rho;
    fact1.c4=rho;
    
    fact2.c1=-1.0f*kappa5;
    fact2.c2=fact2.c1;
    fact2.c3=fact2.c1;
    fact2.c4=fact2.c1;
    
    /* coefficient of rho/2*(1 +- Gamma[0]) */
    
    fact3.c1=2.0f;
    fact3.c2=fact3.c1;
    fact3.c3=fact3.c1;
    fact3.c4=fact3.c1;
    
    fact4.c1=-1.0f;
    fact4.c2=-1.0f;
    fact4.c3=-1.0f;
    fact4.c4=-1.0f;
    
    
    
    /************************ loop over all lattice sites ***************/
    
    
    
    for (ix1=start; ix1<stop; ix1++) {
      
      up1=(_sse_su3*) U+4*ix1;
      
      for(int ell=0; ell<L5; ell++) {
	
	s1=psi+ix1*L5+ell;
	_sse_float_prefetch_spinor(s1);
	
	/******************************* direction +0 ***********************/
	
	sm1=psi+(iy1=idw[ix1][0])*L5+ell;
	sp1=psi+iup[ix1][0]*L5+ell;
	_sse_float_prefetch_spinor(sm1);
	
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
	
	
	/******************************* direction -0 ***********************/
	
	sp1=psi+iup[ix1][1]*L5+ell;
	_sse_float_prefetch_spinor(sp1);
	
	_sse_float_pair_load((*sm1).c1,(*sm1).c2);
	_sse_float_vector_mul(fact3);          
	_sse_float_su3_inverse_multiply((*um1));
	
	_sse_float_vector_load(r12_1);
	_sse_float_vector_add();
	_sse_float_vector_store(r12_1);
	
	up1++;
	_sse_float_prefetch_su3(up1);      
      
	/******************************* direction +1 ***********************/
	
	sm1=psi+(iy1=idw[ix1][1])*L5+ell;
	_sse_float_prefetch_spinor(sm1);
    
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
        
	/******************************* direction -1 ***********************/
	
	sp1=psi+iup[ix1][2]*L5;
	_sse_float_prefetch_spinor(sp1);
	
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
        
	/******************************* direction +2 ***********************/
    
	sm1=psi+(iy1=idw[ix1][2])*L5+ell;
	_sse_float_prefetch_spinor(sm1);
	
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
	_sse_float_prefetch_su3(um1);
        
	/******************************* direction -2 ***********************/
    
	sp1=psi+iup[ix1][3]*L5+ell;
	_sse_float_prefetch_spinor(sp1);
	
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
        
	/******************************* direction +3 ***********************/

	sm1=psi+(iy1=idw[ix1][3])*L5+ell;
	_sse_float_prefetch_spinor(sm1);
	
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
    
	/******************************* direction -3 ***********************/
	
	sn1=chi+ix1*L5+ell;      
	_sse_float_prefetch_spinor(sn1);
	
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
	
	
	/******************************** end of loop ***********************/
      }
    }      
#else

      error("NOT IMPLEMENTED SORRY");

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

#endif // if defined(USE_DOUBLE_PRECISION)
 
  }
};

#endif // if defined(SSE2)


