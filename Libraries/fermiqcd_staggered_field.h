/////////////////////////////////////////////////////////////////
/// @file fermiqcd_staggered_field.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Stuff for SSE/SSE2 compile with -DSSE2
///
/// Distributed under GPL2 License
/// 
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////


/// @brief staggered fermionic field
///
/// Example:
/// @verbatim
/// staggered_field psi(lattice,nc);
/// mdp_site x(lattice);
/// forallsites(x)
///   for(int i=0; i<nc; i++)
///     psi(x,i)=0.0+0.0*I;
/// @endverbatim
class staggered_field: public mdp_complex_field {
 public:
  int nc, ndim, nspin;

  staggered_field(mdp_lattice &a, int nc_, int nspin_=4) {
    // attention here that nspin_ is ignored!
    nc=nc_;
    nspin=nspin_;
    ndim=a.ndim;
    allocate_field(a,nc);
  }
  staggered_field(const staggered_field &chi) : mdp_complex_field(chi) {
    nc=chi.nc;
    ndim=chi.ndim;
    nspin=chi.nspin;
  }
  void operator=(const staggered_field &chi) {
    nc=chi.nc;
    ndim=chi.ndim;
    nspin=chi.nspin;
    mdp_complex_field::operator=(chi);
  }
  inline mdp_matrix operator() (site x) {
    return mdp_matrix(address(x),nc,1);
  }
  inline mdp_complex& operator() (site x, int i) { 
    return *(address(x,i));
  }
  inline const mdp_complex& operator() (site x, int i) const { 
    return *(address(x,i));
  }
  void operator= (mdp_complex a) {
    for(mdp_int i=0; i<size; i++) m[i]=a;
  }
   
  inline mdp_real component(site x, int mu) { 
    return x(mu) % 2;
  }
  inline mdp_real eta(site x, int mu) { 
#ifdef USE_GOLTERMAN
    int i, tmp, i_max=(mu+ndim-1) % ndim;
    tmp=0;
    for(i=1; i<=i_max; i++) tmp+=x(i);
#else
    int i, tmp;
    tmp=0;
    for(i=0; i<mu; i++) tmp+=x(i);
#endif  
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


};

