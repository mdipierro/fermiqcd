/////////////////////////////////////////////////////////////////
/// @file fermiqcd_gauge_routines.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Various gauge multiplication routines (some call SSE/SSE2 macors)
/// 
/// Distributed under GPL2 License
/// 
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////

// ///////////////////////////
// If use class mdp_matrix
// ///////////////////////////
inline mdp_matrix staple(gauge_field &U, register site x, 
		     int mu, int s1, int nu) {
  mdp_matrix tmp(U.nc,U.nc);
    if(s1==+1) {
      tmp=U(x,nu)*U(x+nu,mu)*hermitian(U(x+mu,nu));
    } else {
      site y(U.lattice());
      y=x-nu;
      tmp=hermitian(U(y,nu))*U(y,mu)*U(y+mu,nu);
    }
    return tmp;
}
inline mdp_matrix staple(gauge_field &U, register site x, int mu) {
  mdp_matrix tmp(U.nc,U.nc);
  site y(U.lattice());
  int nu;
  tmp=0;
    for(nu=0; nu<U.ndim; nu++) if(nu!=mu) {
      tmp+=U(x,nu)*U(x+nu,mu)*hermitian(U(x+mu,nu));
      y=x-nu;
      tmp+=hermitian(U(y,nu))*U(y,mu)*U(y+mu,nu);
    }
    return tmp;
  }

inline mdp_matrix staple_H(gauge_field &U, register site x, 
		       int mu, int s1, int nu) {
  mdp_matrix tmp(U.nc,U.nc);
  if(s1==+1) {
    tmp=U(x+mu,nu)*hermitian(U(x+nu,mu))*hermitian(U(x,nu));
  } else {
      site y(U.lattice());
      y=x-nu;
      tmp=hermitian(U(y+mu,nu))*hermitian(U(y,mu))*U(y,nu);
  }
  return tmp;
}
inline mdp_matrix staple_H(gauge_field &U, register site x, int mu) {
  mdp_matrix tmp(U.nc,U.nc);
  site y(U.lattice());
  int nu;
  tmp=0;
  for(nu=0; nu<U.ndim; nu++) if(nu!=mu) {
    tmp+=U(x+mu,nu)*hermitian(U(x+nu,mu))*hermitian(U(x,nu));
    y=x-nu;
    tmp+=hermitian(U(y+mu,nu))*hermitian(U(y,mu))*U(y,nu);
  }
  return tmp;
}

inline mdp_matrix staple_0i_H(gauge_field &U, register site x, int mu) {
  mdp_matrix tmp(U.nc,U.nc);
  site y(U.lattice());
  int nu;
  if(mu==0) {
    tmp=0;
    for(nu=1; nu<U.ndim; nu++) {
      tmp+=U(x+mu,nu)*hermitian(U(x+nu,mu))*hermitian(U(x,nu));
      y=x-nu;
      tmp+=hermitian(U(y+mu,nu))*hermitian(U(y,mu))*U(y,nu);
    }
  } else {
    nu=0;
    tmp=U(x+mu,nu)*hermitian(U(x+nu,mu))*hermitian(U(x,nu));
    y=x-nu;
    tmp+=hermitian(U(y+mu,nu))*hermitian(U(y,mu))*U(y,nu);
    }
  
  return tmp;
}  

inline mdp_matrix staple_ij_H(gauge_field &U, register site x, int mu) {
  mdp_matrix tmp(U.nc,U.nc);
  site y(U.lattice());
  int nu;
  tmp=0;
  if(mu!=0)
    for(nu=1; nu<U.ndim; nu++) if(nu!=mu) {
      tmp+=U(x+mu,nu)*hermitian(U(x+nu,mu))*hermitian(U(x,nu));
	y=x-nu;
	tmp+=hermitian(U(y+mu,nu))*hermitian(U(y,mu))*U(y,nu);
    }
  return tmp;
}

// ///////////////////////////
// plaquette
// ////////////////////////// 
inline mdp_matrix plaquette(gauge_field &U, site x, int mu, int nu) {
  mdp_matrix tmp(U.nc,U.nc);
  tmp=U(x,mu)*U(x+mu,nu)*hermitian(U(x+nu,mu))*hermitian(U(x,nu));
  return tmp;
}

/* obsolete stuff 

void fast_mul_AB_to_C(mdp_complex *a, mdp_complex *b, mdp_complex *c, 
		      int &ni) {
  register int i,j,k, ink, inj;
  for(i=0; i<ni; i++) {
    ink=i*ni; inj=i*ni;
    for(k=0; k<ni; k++) {
      c[ink+k]=a[inj]*b[k];
      for(j=1; j<ni; j++) c[ink+k]+=a[inj+j]*b[j*ni+k];
    }
  }
}
void fast_mul_ABH_to_C(mdp_complex *a, mdp_complex *b, mdp_complex *c, 
		       int &ni) {
  register int i,j,k, ink, inj;
  for(i=0; i<ni; i++) {
    ink=i*ni; inj=i*ni;
    for(k=0; k<ni; k++) {
      c[ink+k]=a[inj]*conj(b[k*ni]);
      for(j=1; j<ni; j++) c[ink+k]+=a[inj+j]*conj(b[k*ni+j]);
    }
  }
}
void fast_mul_AHB_to_C(mdp_complex *a, mdp_complex *b, mdp_complex *c, 
		       int &ni) {
  register int i,j,k, ink;
  for(i=0; i<ni; i++) {
    ink=i*ni;
    for(k=0; k<ni; k++) {
      c[ink+k]=conj(a[i])*b[k];
      for(j=1; j<ni; j++) c[ink+k]+=conj(a[j*ni+i])*b[j*ni+k];
    }
  }
}
void fast_mul_AHBH_to_C(mdp_complex *a, mdp_complex *b, mdp_complex *c, 
		       int &ni) {
  register int i,j,k, ink, inj;
  for(i=0; i<ni; i++) {
    ink=i*ni; inj=i*ni;
    for(k=0; k<ni; k++) {
      c[ink+k]=conj(a[i]*b[k*ni]);
      for(j=1; j<ni; j++) c[ink+k]+=conj(a[j*ni+i]*b[k*ni+j]);
    }
  }
}


void fast_mul_AB_addto_C(mdp_complex *a, mdp_complex *b, mdp_complex *c, 
		      int &ni) {
  register int i,j,k, ink, inj;
  for(i=0; i<ni; i++) {
    ink=i*ni; inj=i*ni;
    for(k=0; k<ni; k++) {
      c[ink+k]+=a[inj]*b[k];
      for(j=1; j<ni; j++) c[ink+k]+=a[inj+j]*b[j*ni+k];
    }
  }
}
void fast_mul_ABH_addto_C(mdp_complex *a, mdp_complex *b, mdp_complex *c, 
		       int &ni) {
  register int i,j,k, ink, inj;
  for(i=0; i<ni; i++) {
    ink=i*ni; inj=i*ni;
    for(k=0; k<ni; k++) {
      c[ink+k]+=a[inj]*conj(b[k*ni]);
      for(j=1; j<ni; j++) c[ink+k]+=a[inj+j]*conj(b[k*ni+j]);
    }
  }
}
void fast_mul_AHB_addto_C(mdp_complex *a, mdp_complex *b, mdp_complex *c, 
		       int &ni) {
  register int i,j,k, ink;
  for(i=0; i<ni; i++) {
    ink=i*ni;
    for(k=0; k<ni; k++) {
      c[ink+k]+=conj(a[i])*b[k];
      for(j=1; j<ni; j++) c[ink+k]+=conj(a[j*ni+i])*b[j*ni+k];
    }
  }
}
void fast_add_AB_to_C(mdp_complex *a, mdp_complex *b, mdp_complex *c, 
		      int &ni) {
  for(register int i=0; i<ni*ni; i++) c[i]=a[i]+b[i];
}
void fast_add_AB_addto_C(mdp_complex *a, mdp_complex *b, mdp_complex *c, 
			 int &ni) {
  for(register int i=0; i<ni*ni; i++) c[i]+=a[i]+b[i];
}
void fast_add_aAbB_to_C(mdp_complex c1, mdp_complex *a1, mdp_complex c2, mdp_complex *a2,
			mdp_complex *c, int &ni) {
  for(register int i=0; i<ni*ni; i++) c[i]=c1*a1[i]+c2*a2[i];
}
void fast_add_aAbB_addto_C(mdp_complex c1, mdp_complex *a1, mdp_complex c2, mdp_complex *a2,
			mdp_complex *c, int &ni) {
  for(register int i=0; i<ni*ni; i++) c[i]+=c1*a1[i]+c2*a2[i];
}
void fast_add_aAbB_to_C(mdp_complex c1, mdp_complex *a1, mdp_complex c2, mdp_complex *a2,
			mdp_complex c3, mdp_complex *a3, mdp_complex *c, int &ni) {
  for(register int i=0; i<ni*ni; i++) c[i]=c1*a1[i]+c2*a2[i]+c3*a3[i];
}
void fast_add_aAbB_addto_C(mdp_complex c1, mdp_complex *a1, mdp_complex c2, mdp_complex *a2,
			mdp_complex c3, mdp_complex *a3, mdp_complex *c, int &ni) {
  for(register int i=0; i<ni*ni; i++) c[i]+=c1*a1[i]+c2*a2[i]+c3*a3[i];
}
void fast_scalar_mul_AB_to_C(mdp_complex a, mdp_complex *b, mdp_complex *c, 
			     int &ni) {
  for(register int i=0; i<ni*ni; i++) c[i]=a*b[i];
}
void fast_scalar_mul_AB_addto_C(mdp_complex a, mdp_complex *b, mdp_complex *c, 
				int &ni) {
  for(register int i=0; i<ni*ni; i++) c[i]+=a*b[i];
}




void fast_mul_AB_subto_C(mdp_complex *a, mdp_complex *b, mdp_complex *c, 
		      int &ni) {
  register int i,j,k, ink, inj;
  for(i=0; i<ni; i++) {
    ink=i*ni; inj=i*ni;
    for(k=0; k<ni; k++) {
      c[ink+k]+=a[inj]*b[k];
      for(j=1; j<ni; j++) c[ink+k]-=a[inj+j]*b[j*ni+k];
    }
  }
}
void fast_mul_ABH_subto_C(mdp_complex *a, mdp_complex *b, mdp_complex *c, 
		       int &ni) {
  register int i,j,k, ink, inj;
  for(i=0; i<ni; i++) {
    ink=i*ni; inj=i*ni;
    for(k=0; k<ni; k++) {
      c[ink+k]+=a[inj]*conj(b[k*ni]);
      for(j=1; j<ni; j++) c[ink+k]-=a[inj+j]*conj(b[k*ni+j]);
    }
  }
}
void fast_mul_AHB_subto_C(mdp_complex *a, mdp_complex *b, mdp_complex *c, 
		       int &ni) {
  register int i,j,k, ink;
  for(i=0; i<ni; i++) {
    ink=i*ni;
    for(k=0; k<ni; k++) {
      c[ink+k]+=conj(a[i])*b[k];
      for(j=1; j<ni; j++) c[ink+k]-=conj(a[j*ni+i])*b[j*ni+k];
    }
  }
}
void fast_add_AB_subto_C(mdp_complex *a, mdp_complex *b, mdp_complex *c, 
			 int &ni) {
  for(register int i=0; i<ni*ni; i++) c[i]-=a[i]+b[i];
}
void fast_add_aAbB_subto_C(mdp_complex c1, mdp_complex *a1, mdp_complex c2, mdp_complex *a2,
			mdp_complex *c, int &ni) {
  for(register int i=0; i<ni*ni; i++) c[i]-=c1*a1[i]+c2*a2[i];
}
void fast_add_aAbB_subto_C(mdp_complex c1, mdp_complex *a1, mdp_complex c2, mdp_complex *a2,
			mdp_complex c3, mdp_complex *a3, mdp_complex *c, int &ni) {
  for(register int i=0; i<ni*ni; i++) c[i]-=c1*a1[i]+c2*a2[i]+c3*a3[i];
}
void fast_scalar_mul_AB_subto_C(mdp_complex a, mdp_complex *b, mdp_complex *c, 
				int &ni) {
  for(register int i=0; i<ni*ni; i++) c[i]-=a*b[i];
}

// ///////////////////////////
// else use fast_mul functions
// gains a factor 1.45 in time
  // ///////////////////////////

inline mdp_matrix staple(gauge_field &U, register site x, 
		     int mu, int s1, int nu) {
  int nc=U.nc;
  mdp_matrix b1(nc,nc);
  site y(U.lattice());
  mdp_matrix tmp(nc,nc);
  if(s1==+1) {
    fast_mul_AB_to_C(&U(x,nu,0,0),&U(x+nu,mu,0,0),&b1(0,0),nc);
      fast_mul_ABH_to_C(&b1(0,0),&U(x+mu,nu,0,0),&tmp(0,0),nc);
  } else {
    y=x-nu;
    fast_mul_AB_to_C(&U(y,mu,0,0),&U(y+mu,nu,0,0),&b1(0,0),nc);
    fast_mul_AHB_to_C(&U(y,nu,0,0),&b1(0,0),&tmp(0,0),nc);
  }
    return tmp;
}
inline mdp_matrix staple(gauge_field &U, register site x, int mu) {
  int nc=U.nc;
    mdp_matrix b1(nc,nc);
    mdp_matrix tmp(nc,nc);
    site y(U.lattice());
    int nu;
    tmp=0;
    for(nu=0; nu<U.ndim; nu++) if(nu!=mu) {
      fast_mul_AB_to_C(&U(x,nu,0,0),&U(x+nu,mu,0,0),&b1(0,0),nc);
      fast_mul_ABH_addto_C(&b1(0,0),&U(x+mu,nu,0,0),&tmp(0,0),nc);
      y=x-nu;
      fast_mul_AB_to_C(&U(y,mu,0,0),&U(y+mu,nu,0,0),&b1(0,0),nc);
      fast_mul_AHB_addto_C(&U(y,nu,0,0),&b1(0,0),&tmp(0,0),nc);
    }
    return tmp;
}

inline mdp_matrix staple_H(gauge_field &U, register site x, 
			      int mu, int s1, int nu) {
  int nc=U.nc;
  mdp_matrix b1(nc,nc);
  mdp_matrix tmp(nc,nc);
  site y(U.lattice());
  if(s1==+1) {
    fast_mul_ABH_to_C(&U(x+mu,nu,0,0),&U(x+nu,mu,0,0),&b1(0,0),nc);
    fast_mul_ABH_to_C(&b1(0,0),&U(x,nu,0,0),&tmp(0,0),nc);
  } else {
    y=x-nu;
    fast_mul_AHB_to_C(&U(y,mu,0,0),&U(y,nu,0,0),&b1(0,0),nc);
    fast_mul_AHB_to_C(&U(y+mu,nu,0,0),&b1(0,0),&tmp(0,0),nc);
  }
  return tmp;
}
inline mdp_matrix staple_H(gauge_field &U, register site x, int mu) {
  int nc=U.nc;
  mdp_matrix b1(nc,nc);
  mdp_matrix tmp(nc,nc);
  site y(U.lattice());
  int nu;
  tmp=0;
  for(nu=0; nu<U.ndim; nu++) if(nu!=mu) {
    fast_mul_ABH_to_C(&U(x+mu,nu,0,0),&U(x+nu,mu,0,0),&b1(0,0),nc);
    fast_mul_ABH_addto_C(&b1(0,0),&U(x,nu,0,0),&tmp(0,0),nc);
    y=x-nu;
    fast_mul_AHB_to_C(&U(y,mu,0,0),&U(y,nu,0,0),&b1(0,0),nc);
    fast_mul_AHB_addto_C(&U(y+mu,nu,0,0),&b1(0,0),&tmp(0,0),nc);
  }
  return tmp;
}

inline mdp_matrix staple_0i_H(gauge_field &U, register site x, int mu) {
  int nc=U.nc;
  mdp_matrix b1(nc,nc);
  mdp_matrix tmp(U.nc,U.nc);
  site y(U.lattice());
  int nu;
  if(mu==0) {
    tmp=0;
    for(nu=1; nu<U.ndim; nu++) {
      fast_mul_ABH_to_C(&U(x+mu,nu,0,0),&U(x+nu,mu,0,0),&b1(0,0),nc);
      fast_mul_ABH_addto_C(&b1(0,0),&U(x,nu,0,0),&tmp(0,0),nc);
      y=x-nu;
      fast_mul_AHB_to_C(&U(y,mu,0,0),&U(y,nu,0,0),&b1(0,0),nc);
      fast_mul_AHB_addto_C(&U(y+mu,nu,0,0),&b1(0,0),&tmp(0,0),nc);
    }
  } else {
    nu=0;
    fast_mul_ABH_to_C(&U(x+mu,nu,0,0),&U(x+nu,mu,0,0),&b1(0,0),nc);
    fast_mul_ABH_to_C(&b1(0,0),&U(x,nu,0,0),&tmp(0,0),nc);
    y=x-nu;
    fast_mul_AHB_to_C(&U(y,mu,0,0),&U(y,nu,0,0),&b1(0,0),nc);
    fast_mul_AHB_addto_C(&U(y+mu,nu,0,0),&b1(0,0),&tmp(0,0),nc);
  }
  return tmp;
}  

inline mdp_matrix staple_ij_H(gauge_field &U, register site x, int mu) {
  int nc=U.nc;
  mdp_matrix b1(nc,nc);
  mdp_matrix tmp(U.nc,U.nc);
  site y(U.lattice());
  int nu;
  tmp=0;
  if(mu!=0)
    for(nu=1; nu<U.ndim; nu++) if(nu!=mu) {
      fast_mul_ABH_to_C(&U(x+mu,nu,0,0),&U(x+nu,mu,0,0),&b1(0,0),nc);
      fast_mul_ABH_addto_C(&b1(0,0),&U(x,nu,0,0),&tmp(0,0),nc);
      y=x-nu;
      fast_mul_AHB_to_C(&U(y,mu,0,0),&U(y,nu,0,0),&b1(0,0),nc);
      fast_mul_AHB_addto_C(&U(y+mu,nu,0,0),&b1(0,0),&tmp(0,0),nc);
    }
  return tmp;
}

// ///////////////////////////
// plaquette
// ////////////////////////// 
inline mdp_matrix plaquette(gauge_field &U, site x, int mu, int nu) {
  int nc=U.nc;
  mdp_matrix b1(nc,nc);
  mdp_matrix b2(nc,nc);
  mdp_matrix tmp(U.nc,U.nc);
  fast_mul_AB_to_C(&U(x,mu,0,0),&U(x+mu,nu,0,0),&b1(0,0),nc);
  fast_mul_ABH_to_C(&b1(0,0),&U(x+nu,mu,0,0),&b2(0,0),nc);
  fast_mul_ABH_to_C(&b2(0,0),&U(x,nu,0,0),&tmp(0,0),nc);
  return tmp;
}
*/

