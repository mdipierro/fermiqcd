/////////////////////////////////////////////////////////////////
/// @file fermiqcd_staggered_actions_sse2.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Stuff for SSE/SSE2 compile with -DSSE2
///
/// Distributed under GPL2 License
/// 
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////

#if defined(SSE2)

/// @brief Staggered/Asqtad action optimized for Pentium 4
///
/// Compile with -DSSE2
///
/// Example:
/// @verbatim
/// gauge_field U(lattice,nc);
/// staggered_field psi(lattice,nc);
/// staggered_field chi(lattice,nc);
/// coefficients coeff;
/// coeff["mass"]=2.0;
/// default_staggered_action=StaggeredAsqtadActionSSE2::mul_Q;
/// mul_Q(chi,psi,U,coeff);
/// @endverbatim
class StaggeredAsqtadActionSSE2 {
 public:
  static void mul_Q(staggered_field &chi_out,
		    staggered_field &psi_in,
		    gauge_field &U_in,
		    coefficients &coeff, 
		    int parity=EVENODD) {
    
    _sse_su3_vector* chi=(_sse_su3_vector*) chi_out.physical_address();
    _sse_su3_vector* psi=(_sse_su3_vector*) psi_in.physical_address();
    
    _sse_su3*  U=(_sse_su3*) U_in.physical_address();
    _sse_su3*  LL=(_sse_su3*) U_in.long_links.physical_address();
    mdp_int** iup=U_in.lattice().up;
    mdp_int** idw=U_in.lattice().dw;
    
    mdp_int   start=U_in.lattice().start_index(ME,parity); 
    mdp_int   stop =U_in.lattice().stop_index(ME,parity);  
    
    mdp_real two_mass;
    int sign;
    if(coeff.has_key("mass")) two_mass=2.0*coeff["mass"];
    else error("coefficient mass undefined");
    if(coeff.has_key("sign")) sign=(int) coeff["sign"];
    else sign=1;

    site x(psi_in.lattice());
    
#if defined(USE_DOUBLE_PRECISION) 

  static _sse_double fact1 ALIGN16;
  static _sse_su3_vector r ALIGN16;
  _sse_su3 *up, *um;   
  _sse_su3_vector *s,*sp,*sm,*sn;
  mdp_int ix,iy,iz;
  int sign0, sign1, sign2, sign3;

  fact1.c1=1.0f*two_mass;
  fact1.c2=1.0f*two_mass;

  sp=psi+iup[start][0];
  up=U+4*start;

  /************************ loop over all lattice sites ***************/
  
  for(ix=start; ix<stop; ix++) {

    s=psi+ix;
    
    x.idx=ix;
    if(sign>0) {
      sign0=(int) (psi_in.eta(x,0)+1);
      sign1=(int) (psi_in.eta(x,1)+1);
      sign2=(int) (psi_in.eta(x,2)+1);
      sign3=(int) (psi_in.eta(x,3)+1);
    } else {
      sign0=(int) (1.0-psi_in.eta(x,0));
      sign1=(int) (1.0-psi_in.eta(x,1));
      sign2=(int) (1.0-psi_in.eta(x,2));
      sign3=(int) (1.0-psi_in.eta(x,3));
    }

    /**************** mu =0 *****************/

    iy=idw[ix][0];
    sm=psi+iy;
    _sse_float_prefetch_spinor(sm);
    um=(_sse_su3*) U+iy*4;
    _sse_float_prefetch_su3(um);

    _sse_double_load(*sp);
    _sse_double_su3_multiply(*up);
    _sse_double_load(*s);
    _sse_double_vector_mul(fact1);
    if(sign0) _sse_double_vector_add();
    else       _sse_double_vector_sub();
    _sse_double_store(r);
    
    up++;
    sp=psi+iup[ix][1];
    _sse_float_prefetch_spinor(sp);
    
    _sse_double_load(*sm);
    _sse_double_su3_inverse_multiply(*um);
    _sse_double_load(r);
    if(sign0) _sse_double_vector_sub();
    else       _sse_double_vector_add();
    _sse_double_store(r);


    /**************** mu =1 *****************/

    iy=idw[ix][1];
    sm=psi+iy;
    _sse_float_prefetch_spinor(sm);
    um=U+iy*4+1;
    _sse_float_prefetch_su3(um);

    _sse_double_load(*sp);
    _sse_double_su3_multiply(*up);
    _sse_double_load(r);
    if(sign1) _sse_double_vector_add();
    else       _sse_double_vector_sub();
    _sse_double_store(r);

    up++;
    sp=psi+iup[ix][2];
    _sse_float_prefetch_spinor(sp);

    _sse_double_load(*sm);
    _sse_double_su3_inverse_multiply(*um);
    _sse_double_load(r);
    if(sign1) _sse_double_vector_sub();
    else       _sse_double_vector_add();
    _sse_double_store(r);

    /**************** mu = 2 *****************/

    iy=idw[ix][2];
    sm=psi+iy;
    _sse_float_prefetch_spinor(sm);
    um=U+iy*4+2;
    _sse_float_prefetch_su3(um);

    _sse_double_load(*sp);
    _sse_double_su3_multiply(*up);
    _sse_double_load(r);
    if(sign2) _sse_double_vector_add();
    else       _sse_double_vector_sub();
    _sse_double_store(r);

    up++;
    sp=psi+iup[ix][3];
    _sse_float_prefetch_spinor(sp);

    _sse_double_load(*sm);
    _sse_double_su3_inverse_multiply(*um);
    _sse_double_load(r);
    if(sign2) _sse_double_vector_sub();
    else       _sse_double_vector_add();
    _sse_double_store(r);
 
    /**************** mu =3 *****************/

    sn=chi+ix;
    _sse_float_prefetch_spinor(sn);

    iy=idw[ix][3];
    sm=psi+iy;
    _sse_float_prefetch_spinor(sm);
    um=U+iy*4+3;
    _sse_float_prefetch_su3(um);

    _sse_double_load(*sp);
    _sse_double_su3_multiply(*up);
    _sse_double_load(r);
    if(sign3) _sse_double_vector_add();
    else       _sse_double_vector_sub();
    _sse_double_store(r);

    iz=ix+1;
    if(iz<stop) {
      sp=psi+iup[iz][0];
      _sse_float_prefetch_spinor(sp);
      up++;
    }
    _sse_double_load(*sm);
    _sse_double_su3_inverse_multiply(*um);
    _sse_double_load(r);
    if(sign3) _sse_double_vector_sub();
    else       _sse_double_vector_add();
    _sse_double_store(*sn);

  }
  
  if(LL==0) return;

  /************************ loop over all lattice sites for Naik term ***************/

  sp=psi+iup[iup[iup[start][0]][0]][0];
  up=LL+4*start;
  
  for(ix=start; ix<stop; ix++) {
    
    sn=chi+ix;
    _sse_float_prefetch_spinor(sn);
    
   x.idx=ix;
    if(sign>0) {
      sign0=(int) (psi_in.eta(x,0)+1);
      sign1=(int) (psi_in.eta(x,1)+1);
      sign2=(int) (psi_in.eta(x,2)+1);
      sign3=(int) (psi_in.eta(x,3)+1);
    } else {
      sign0=(int) (1.0-psi_in.eta(x,0));
      sign1=(int) (1.0-psi_in.eta(x,1));
      sign2=(int) (1.0-psi_in.eta(x,2));
      sign3=(int) (1.0-psi_in.eta(x,3));
    }

    /**************** mu =0 *****************/

    iy=idw[idw[idw[ix][0]][0]][0];
    
    sm=psi+iy;
    _sse_float_prefetch_spinor(sm);
    um=(_sse_su3*) LL+iy*4;
    _sse_float_prefetch_su3(um);
    
    _sse_double_load(*sp);
    _sse_double_su3_multiply(*up);
    _sse_double_load(*sn);
    if(sign0) _sse_double_vector_add();
    else      _sse_double_vector_sub();
    _sse_double_store(r);
    
    up++;
      sp=psi+iup[iup[iup[ix][1]][1]][1];  
    _sse_float_prefetch_spinor(sp);
    
    _sse_double_load(*sm);
    _sse_double_su3_inverse_multiply(*um);
    _sse_double_load(r);
    if(sign0) _sse_double_vector_sub();
    else       _sse_double_vector_add();
    _sse_double_store(r);


    /**************** mu =1 *****************/

    iy=idw[idw[idw[ix][1]][1]][1];
    
    sm=psi+iy;
    _sse_float_prefetch_spinor(sm);
    um=LL+iy*4+1;
    _sse_float_prefetch_su3(um);

    _sse_double_load(*sp);
    _sse_double_su3_multiply(*up);
    _sse_double_load(r);
    if(sign1) _sse_double_vector_add();
    else       _sse_double_vector_sub();
    _sse_double_store(r);

    up++;
    sp=psi+iup[iup[iup[ix][2]][2]][2];
    
    _sse_float_prefetch_spinor(sp);

    _sse_double_load(*sm);
    _sse_double_su3_inverse_multiply(*um);
    _sse_double_load(r);
    if(sign1) _sse_double_vector_sub();
    else       _sse_double_vector_add();
    _sse_double_store(r);

    /**************** mu = 2 *****************/

    iy=idw[idw[idw[ix][2]][2]][2];
   
    sm=psi+iy;
    _sse_float_prefetch_spinor(sm);
    um=LL+iy*4+2;
    _sse_float_prefetch_su3(um);

    _sse_double_load(*sp);
    _sse_double_su3_multiply(*up);
    _sse_double_load(r);
    if(sign2) _sse_double_vector_add();
    else       _sse_double_vector_sub();
    _sse_double_store(r);

    up++;
    sp=psi+iup[iup[iup[ix][3]][3]][3];   
    _sse_float_prefetch_spinor(sp);

    _sse_double_load(*sm);
    _sse_double_su3_inverse_multiply(*um);
    _sse_double_load(r);
    if(sign2) _sse_double_vector_sub();
    else       _sse_double_vector_add();
    _sse_double_store(r);
 
    /**************** mu =3 *****************/

    _sse_float_prefetch_spinor(sn);

    iy=idw[idw[idw[ix][3]][3]][3];
  
    sm=psi+iy;
    _sse_float_prefetch_spinor(sm);
    um=LL+iy*4+3;
    _sse_float_prefetch_su3(um);

    _sse_double_load(*sp);
    _sse_double_su3_multiply(*up);
    _sse_double_load(r);
    if(sign3) _sse_double_vector_add();
    else       _sse_double_vector_sub();
    _sse_double_store(r);

    iz=ix+1;
    if(iz<stop) {
      sp=psi+iup[iup[iup[iz][0]][0]][0];
      _sse_float_prefetch_spinor(sp);
      up++;
    }
    _sse_double_load(*sm);
    _sse_double_su3_inverse_multiply(*um);
    _sse_double_load(r);
    if(sign3) _sse_double_vector_sub();
    else       _sse_double_vector_add();
    _sse_double_store(*sn);

  }

#else

  static _sse_float fact1 ALIGN16;
  static _sse_vector r   ALIGN16;
  _sse_su3 *up, *um;   
  _sse_su3_vector *s,*sp,*sm,*sn;
  _sse_su3_vector dump;
  mdp_int ix,iy,iz;
  int *peta;
  int sign0, sign1, sign2, sign3;

  _sse_check_alignment(&fact1, 0xf);
  _sse_check_alignment(&r, 0xf);

  fact1.c1=1.0f*two_mass;
  fact1.c2=1.0f*two_mass;
  fact1.c3=1.0f*two_mass;
  fact1.c4=1.0f*two_mass;

  dump.c1.real()=0.0f;
  dump.c1.imag()=0.0f;
  dump.c2.real()=0.0f;
  dump.c2.imag()=0.0f;
  dump.c3.real()=0.0f;
  dump.c3.imag()=0.0f;

  sp=psi+iup[start][0];
  up=U+4*start;

  /************************ loop over all lattice sites ***************/
  
  for(ix=start; ix<stop; ix++) {
    
    s=psi+ix;

    x.idx=ix;
	     
    if(sign>0) {
      sign0=(int) (psi_in.eta(x,0)+1);
      sign1=(int) (psi_in.eta(x,1)+1);
      sign2=(int) (psi_in.eta(x,2)+1);
      sign3=(int) (psi_in.eta(x,3)+1);
    } else {
      sign0=(int) (1.0-psi_in.eta(x,0));
      sign1=(int) (1.0-psi_in.eta(x,1));
      sign2=(int) (1.0-psi_in.eta(x,2));
      sign3=(int) (1.0-psi_in.eta(x,3));
    }

    /**************** mu =0 *****************/

    iy=idw[ix][0];
    sm=psi+iy;
    _sse_float_prefetch_spinor(sm);
    um=(_sse_su3*) U+iy*4;
    _sse_float_prefetch_su3(um);
    
    _sse_float_pair_load(*sp,dump);
    _sse_float_su3_multiply(*up);
    _sse_float_pair_load(*s,dump);
    _sse_float_vector_mul(fact1);
    if(sign0) _sse_float_vector_add();
    else       _sse_float_vector_sub();
    _sse_float_vector_store(r);
    
    up++;
    sp=psi+iup[ix][1];
    _sse_float_prefetch_spinor(sp);
    
    _sse_float_pair_load(*sm,dump);
    _sse_float_su3_inverse_multiply(*um);
    _sse_float_vector_load(r);
    if(sign0) _sse_float_vector_sub();
    else       _sse_float_vector_add();
    _sse_float_vector_store(r);


    /**************** mu =1 *****************/

    iy=idw[ix][1];
    sm=psi+iy;
    _sse_float_prefetch_spinor(sm);
    um=U+iy*4+1;
    _sse_float_prefetch_su3(um);

    _sse_float_pair_load(*sp,dump);
    _sse_float_su3_multiply(*up);
    _sse_float_vector_load(r);
    if(sign1) _sse_float_vector_add();
    else       _sse_float_vector_sub();
    _sse_float_vector_store(r);

    up++;
    sp=psi+iup[ix][2];
    _sse_float_prefetch_spinor(sp);

    _sse_float_pair_load(*sm,dump);
    _sse_float_su3_inverse_multiply(*um);
    _sse_float_vector_load(r);
    if(sign1) _sse_float_vector_sub();
    else       _sse_float_vector_add();
    _sse_float_vector_store(r);

    /**************** mu = 2 *****************/

    iy=idw[ix][2];
    sm=psi+iy;
    _sse_float_prefetch_spinor(sm);
    um=U+iy*4+2;
    _sse_float_prefetch_su3(um);

    _sse_float_pair_load(*sp,dump);
    _sse_float_su3_multiply(*up);
    _sse_float_vector_load(r);
    if(sign2) _sse_float_vector_add();
    else       _sse_float_vector_sub();
    _sse_float_vector_store(r);

    up++;
    sp=psi+iup[ix][3];
    _sse_float_prefetch_spinor(sp);

    _sse_float_pair_load(*sm,dump);
    _sse_float_su3_inverse_multiply(*um);
    _sse_float_vector_load(r);
    if(sign2) _sse_float_vector_sub();
    else       _sse_float_vector_add();
    _sse_float_vector_store(r);
 
    /**************** mu =3 *****************/

    sn=chi+ix;
    _sse_float_prefetch_spinor(sn);

    iy=idw[ix][3];
    sm=psi+iy;
    _sse_float_prefetch_spinor(sm);
    um=U+iy*4+3;
    _sse_float_prefetch_su3(um);


    _sse_float_pair_load(*sp,dump);
    _sse_float_su3_multiply(*up);
    _sse_float_vector_load(r);
    if(sign3) _sse_float_vector_add();
    else       _sse_float_vector_sub();
    _sse_float_vector_store(r);

    iz=ix+1;
    if(iz<stop) {
      sp=psi+iup[iz][0];
      _sse_float_prefetch_spinor(sp);
      up++;
    }

    _sse_float_pair_load(*sm,dump);
    _sse_float_su3_inverse_multiply(*um);
    _sse_float_vector_load(r);
    if(sign3) _sse_float_vector_sub();
    else       _sse_float_vector_add();
    _sse_float_pair_store(*sn,dump);

  }
  
  if(LL==0) return;

  /************************ loop over all lattice sites for Naik term ***************/

  sp=psi+iup[iup[iup[start][0]][0]][0];
  up=LL+4*start;
  
  for(ix=start; ix<stop; ix++) {
    
    sn=chi+ix;
    _sse_float_prefetch_spinor(sn);
    
    x.idx=ix;
    if(sign>0) {
      sign0=(int) (psi_in.eta(x,0)+1);
      sign1=(int) (psi_in.eta(x,1)+1);
      sign2=(int) (psi_in.eta(x,2)+1);
      sign3=(int) (psi_in.eta(x,3)+1);
    } else {
      sign0=(int) (1.0-psi_in.eta(x,0));
      sign1=(int) (1.0-psi_in.eta(x,1));
      sign2=(int) (1.0-psi_in.eta(x,2));
      sign3=(int) (1.0-psi_in.eta(x,3));
    }   

    /**************** mu =0 *****************/

    iy=idw[idw[idw[ix][0]][0]][0];

    sm=psi+iy;
    _sse_float_prefetch_spinor(sm);
    um=(_sse_su3*) LL+iy*4;
    _sse_float_prefetch_su3(um);
    
    _sse_float_pair_load(*sp,dump);
    _sse_float_su3_multiply(*up);
    _sse_float_pair_load(*sn,dump);
    if(sign0) _sse_float_vector_add();
    else       _sse_float_vector_sub();
    _sse_float_vector_store(r);
    
    up++; 
    sp=psi+iup[iup[iup[ix][1]][1]][1];
    _sse_float_prefetch_spinor(sp);
    
    _sse_float_pair_load(*sm,dump);
    _sse_float_su3_inverse_multiply(*um);
    _sse_float_vector_load(r);
    if(sign0) _sse_float_vector_sub();
    else       _sse_float_vector_add();
    _sse_float_vector_store(r);


    /**************** mu =1 *****************/

    iy=idw[idw[idw[ix][1]][1]][1];

    sm=psi+iy;
    _sse_float_prefetch_spinor(sm);
    um=LL+iy*4+1;
    _sse_float_prefetch_su3(um);

    _sse_float_pair_load(*sp,dump);
    _sse_float_su3_multiply(*up);
    _sse_float_vector_load(r);
    if(sign1) _sse_float_vector_add();
    else       _sse_float_vector_sub();
    _sse_float_vector_store(r);

    up++;
    sp=psi+iup[iup[iup[ix][2]][2]][2];
    _sse_float_prefetch_spinor(sp);

    _sse_float_pair_load(*sm,dump);
    _sse_float_su3_inverse_multiply(*um);
    _sse_float_vector_load(r);
    if(sign1) _sse_float_vector_sub();
    else       _sse_float_vector_add();
    _sse_float_vector_store(r);

    /**************** mu = 2 *****************/

    iy=idw[idw[idw[ix][2]][2]][2];
    sm=psi+iy;
    _sse_float_prefetch_spinor(sm);
    um=LL+iy*4+2;
    _sse_float_prefetch_su3(um);

    _sse_float_pair_load(*sp,dump);
    _sse_float_su3_multiply(*up);
    _sse_float_vector_load(r);
    if(sign2) _sse_float_vector_add();
    else       _sse_float_vector_sub();
    _sse_float_vector_store(r);

    up++;
    sp=psi+iup[iup[iup[ix][3]][3]][3];
    _sse_float_prefetch_spinor(sp);

    _sse_float_pair_load(*sm,dump);
    _sse_float_su3_inverse_multiply(*um);
    _sse_float_vector_load(r);
    if(sign2) _sse_float_vector_sub();
    else       _sse_float_vector_add();
    _sse_float_vector_store(r);
 
    /**************** mu =3 *****************/

    _sse_float_prefetch_spinor(sn);

    iy=idw[idw[idw[ix][3]][3]][3];

    sm=psi+iy;
    _sse_float_prefetch_spinor(sm);
    um=LL+iy*4+3;
    _sse_float_prefetch_su3(um);

    _sse_float_pair_load(*sp,dump);
    _sse_float_su3_multiply(*up);
    _sse_float_vector_load(r);
    if(sign3) _sse_float_vector_add();
    else       _sse_float_vector_sub();
    _sse_float_vector_store(r);

    iz=ix+1;
    if(iz<stop) {
      sp=psi+iup[iup[iup[iz][0]][0]][0];
      _sse_float_prefetch_spinor(sp);
      up++;
    }
    _sse_float_pair_load(*sm,dump);
    _sse_float_su3_inverse_multiply(*um);
    _sse_float_vector_load(r);
    if(sign3) _sse_float_vector_sub();
    else       _sse_float_vector_add();
    _sse_float_pair_store(*sn,dump);

  }
#endif
  }
};

#endif // id fefined(SSE2)


