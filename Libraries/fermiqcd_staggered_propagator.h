/////////////////////////////////////////////////////////////////
/// @file fermiqcd_staggered_propagator.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Various stuff for staggered fermions
///
/// Distributed under GPL2 License
/// 
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////

/// @brief staggared quark propagator
/// 
/// On a (2n) dimensional lattice this makes 3*(2^n) sources at the 
/// vertices of the hypercube at the origin of the lattice
/// and inverts the Staggered/Asqtad action on them.
///
/// Example:
/// @verbatim
/// mdp_gauge U(lattice,nc);
/// staggered_propagator S(lattice,nc);
/// mdp_site x(lattice);
/// mdp_site y(lattice);
/// coefficients coeff;
/// coeff["mass"]=1.0;
/// generate(S,U,coeff);
/// for(int i=0; i<(int) pow(2,lattice.ndim); i++) {
///    x=binary2versor(a);
///    cout << "source at:" << x << "\nprop:\n"; 
///    forallsites(y) cout << S(x,a) << endl;
/// }
/// @endverbatim
class staggered_propagator: public mdp_field<mdp_complex> {
 public:
  int nc;
  staggered_propagator(mdp_lattice &mylattice, int nc_) {
    nc=nc_;
    int ndim=mylattice.ndim;
    allocate_field(mylattice, ndim*ndim*nc*nc);
  }
  inline mdp_matrix operator() (site x, int a) {
    mdp_matrix tmp(address(x,a*nc*nc),nc,nc);
    return tmp;
  }
  inline mdp_complex &operator() (site x, int a, int i, int j) {
    return *(address(x)+a*nc*nc+i*nc+j);
  }
  friend void generate(staggered_propagator &S, gauge_field &U,
		       coefficients &coeff,
		       mdp_real absolute_precision=fermi_inversion_precision,
		       mdp_real relative_precision=0,
		       int max_steps=2000,
		       void (*smf)(staggered_field&,gauge_field&)=0, 
		       int comp=0) {
    staggered_field psi(S.lattice(),S.nc);
    staggered_field chi(S.lattice(),S.nc);
    site        x(S.lattice());
    int ndim=S.lattice().ndim;
    int nc=S.nc;
    int i,j,mu,a;
    
    double time=mpi.time();

    if(ME==0 && shutup==false) {
      printf("BEGIN Generating ordinary propagator\n");
      fflush(stdout);
    }

    for(a=0; a<(0x1 << ndim); a++) 
       for(j=0; j<nc; j++) {
	forallsitesandcopies(x) 
	  for(i=0; i<nc; i++) 
	    psi(x,i)=0;
	
	x=binary2versor(a);
	if(ME==0 && shutup==false) {
	  printf("(source at (");
	  for(mu=0; mu<ndim; mu++) printf("%i ", x(mu));
	  printf("), Color: %i\n", j);
	  fflush(stdout);
	}
	if(x.is_here()) psi(x,j)=1;

	/*
	  If a smearing function is passed (smf)
	  the source is smeared before the inversion
	  the sink must be smeared using smear_sink.
	*/      
	if(smf!=0) (*smf)(psi,U);
	mul_invQ(chi,psi,U,coeff,absolute_precision,relative_precision,max_steps);
	
	forallsites(x)
	  for(i=0; i<nc; i++)
	    S(x,a,i,j)=chi(x,i);
       }
    if(ME==0 && shutup==false) {
      printf("END Generating ordinary propagator. Time: %f (sec)\n",
	     mpi.time()-time);
      fflush(stdout);
    }
  }
};
