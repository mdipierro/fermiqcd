/////////////////////////////////////////////////////////////////
/// @file fermiqcd_sdwf_field.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// WORK IN PROGRESS
///
/// Distributed under GPL2 License
/// 
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////

/// @brief field for domain wall staggered fermions 
class sdwf_field: public mdp_complex_field {
 public:
  int nc, ndim, nspin, L5;

  sdwf_field(mdp_lattice &a, int L5_, int nc_, int nspin_=4) {
    // attention here that nspin_ is ignored!
    L5=L5_;
    nc=nc_;
    nspin=nspin_;
    ndim=a.ndim;
    allocate_field(a,L5*nc);
  }
  sdwf_field(sdwf_field &chi) {
    // attention here that nspin_ is ignored!
    nc=chi.nc;
    nspin=chi.nspin;
    L5=chi.L5;
    ndim=chi.lattice().ndim;
    allocate_field(chi.lattice(),L5*nc);
  }
  mdp_matrix operator() (site x, int x5) {
    return mdp_matrix(address(x,x5*nc),nc,1);
  }
  inline mdp_complex& operator() (site x, int x5, int i) { 
    return *(address(x,x5*nc+i));
  }
  inline const mdp_complex& operator() (site x, int x5, int i) const { 
    return *(address(x,x5*nc+i));
  }
  void operator= (mdp_complex a) {
    for(mdp_int i=0; i<size; i++) m[i]=a;
  }  
  inline mdp_real component(site x, int mu) { 
    return x(mu) % 2;
  }
  inline mdp_real eta(site x, int mu) { 
    int tmp;
    int i;
    int i_max=(mu+ndim-1) % ndim;
    tmp=0;
    for(i=1; i<=i_max; i++) tmp+=x(i);
    return mdp_mod2sign(tmp);
  }
  /*
  inline mdp_real zeta(site x, int mu) { 
    
  }
  */
  inline mdp_real eps (site x) { 
    int tmp;
    int i;
    tmp=x(0);
    for(i=1; i<ndim; i++) tmp+=x(i);
    return mdp_mod2sign(tmp);
  }
  inline mdp_real type (site x) { 
    mdp_real tmp;
    int i;
    tmp=x(0) % 2;
    for(i=1; i<ndim; i++) tmp+=(x(i) % 2)*pow(2.0,i);
    return tmp;
  }
  inline site chiral_shift(site x) {
    int i;
    for(i=0; i<ndim; i++) if(x(i)%2==1) x=x-i; else x=x+i;
    return x;
  }
  inline mdp_real chiral_phase(site x) { // (Gamma5 (x) 1)
    int tmp=ndim/2;
    int i;
    for(i=1; i<ndim; i+=2) tmp+=x(i);
    return (mdp_real) mdp_mod2sign(tmp);
  }
  inline mdp_real chiral_phase2(site x) { // (Gamma5 (x) Gamma5)
    int tmp=0;
    int i;
    for(i=0; i<ndim; i++) tmp+=x(i);
    return (mdp_real) mdp_mod2sign(tmp);
  }
};

