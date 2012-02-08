/////////////////////////////////////////////////////////////////
/// @file fermiqcd_fermi_propagator.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains class fermi_propagator
///
/// Distributed under GPL2 License
/// 
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////

/// @brief a Wilson/Clover quark propagator (all 12 components)
///
/// Example of how to make a pion:
/// @verbatim
/// gauge_field U(lattice,nc);
/// U.load("myfield");
/// fermi_propagator S(lattice,nc);
/// coefficients quark;
/// quark["kappa"]=1.12;
/// generate(S,U,quark);
/// vector<float> sum(U.lattice.size(TIME));
/// forallsites(x) 
///   for(int alpha=0; alpha<4; alpha++)
///     for(int beta=0; beta<4; beta++)
///        sum(x(0))+=real(trace(S(x,alpha,beta)*
///                   hermitian(S(x,beta,alpha))));
/// @endverbatim
/// Note that S(x,alpha,beta,i,j) is 
/// \f$ \left<0|\bar q^i_\alpha(x), q^j_\beta(0)|\right> \f$ 
class fermi_propagator: public mdp_complex_field {
public:
  int nspin, nc;
  int gamma5_a[10][10]; 
  int gamma5_b[10][10]; 
  mdp_complex gamma5_val[10][10];
  fermi_propagator() {
    reset_field();
  };
  void set_gamma5(mdp_matrix gamma5) {
    mdp_matrix a(nspin,nspin);
    for(int i=0; i<nspin; i++)
      for(int j=0; j<nspin; j++)
	a(i,j)=i*nspin+j+1;
    a = gamma5 * a * gamma5;
    for(int i=0; i<nspin; i++)
      for(int j=0; j<nspin; j++) {
	int k = int(abs(a(i,j)))-1;
	gamma5_a[i][k]=k/nspin;
	gamma5_b[i][k]=k%nspin;
	gamma5_val[i][j] = a(i,j)/abs(a(i,j));
      }
  }
  fermi_propagator(mdp_lattice &mylattice, int nc_, int nspin_=4) {
    reset_field();
    nspin=nspin_;
    nc=nc_;
    allocate_field(mylattice, nspin*nspin*nc*nc);
  };
  void allocate_fermi_propagator(mdp_lattice &mylattice, 
				 int nc_, int nspin_=4) {
    deallocate_field();
    nspin=nspin_;
    nc=nc_;
    allocate_field(mylattice, nspin*nspin*nc*nc);
  };
  
  inline mdp_matrix operator() (site x, int a, int b) {
    mdp_matrix tmp(address(x,(a*nspin+b)*nc*nc), nc, nc);
    return tmp;
  };
  inline mdp_matrix operator() (int zero, site x, int a, int b) {
    mdp_matrix tmp(address(x,(a*nspin+b)*nc*nc), nc, nc);
    return tmp;
  };
  inline mdp_matrix operator() (site x, int zero, int a, int b) {
    // this needs testing!
    // S(x,zero) = gamma5 * S(zero,x)^h * gamma5
    mdp_matrix tmp(nc,nc);
    int c = gamma5_a[a][b];
    int d = gamma5_b[a][b];
    mdp_complex v = gamma5_val[a][b];
    for(int i=0; i<nc; i++)
      for(int j=0; j<nc; j++)
	tmp(i,j)=v*conj((*this)(x,d,c,j,i));
    return tmp;
  };
  inline mdp_complex &operator() (site x, int a, int b, int i, int j) {
    return *(address(x,((a*nspin+b)*nc+i)*nc+j));
  };
  const mdp_complex operator() (site x, int a, int b, int i, int j) const {
    return *(address(x,((a*nspin+b)*nc+i)*nc+j));
  };
  /// makes the quark propagator
  ///
  /// @param S the output propagator
  /// @param U the input gauge configuration
  /// @param coeff the parameters to be passed to the action
  /// @param absolute_precision the target absolute precision for inversion
  /// @param relative_precision the target relative precision for invcersion
  /// @param max_steps the max number of steps in inversion
  /// @param smf pointer to smearing function (smear sources)
  /// @param smear_coeff parameters for smearing
  friend void generate(fermi_propagator &S, gauge_field &U,
		       coefficients &coeff,
		       mdp_real absolute_precision=fermi_inversion_precision,
		       mdp_real relative_precision=0,
		       int max_steps=2000,
		       void (*smf)(fermi_field&,
				   gauge_field&,
				   coefficients&)=0,
		       coefficients smear_coeff=coefficients(),
		       int comp=0) {
    fermi_field psi(S.lattice(),S.nc,S.nspin);
    fermi_field chi(S.lattice(),S.nc,S.nspin);
    site        x(S.lattice());
    int i,j,a,b;
    double time=mpi.time();
    inversion_stats stats;
    begin_function("generate");
    mdp << "BEGIN Generating ordinary propagator\n";
 
    for(b=0; b<psi.nspin; b++)
      for(j=0; j<psi.nc; j++) {

	mdp << "Source: spin=" << b << ", color=" << j << '\n';

	forallsitesandcopies(x)
	  for(a=0; a<psi.nspin; a++)
	    for(i=0; i<psi.nc; i++) {
	      if((x.is_equal(0)) && (a==b) && (i==j)) psi(x,a,i)=1;
	      else                                    psi(x,a,i)=0;
	    }
	/*
	  If a smearing function is passed (smf)
	  the source is smeared before the inversion
	  the sink must be smeared using smear_sink.
	*/
	if(smf!=0) (*smf)(psi,U,smear_coeff);
	stats=mul_invQ(chi,psi,U,coeff,
		       absolute_precision,relative_precision,max_steps);

	forallsites(x)
	  for(a=0; a<psi.nspin; a++)
	    for(i=0; i<psi.nc; i++) {
	      S(x,a,b,i,j)=chi(x,a,i);
	  }
	cout << "Statistics: residue=" << stats.residue
	     << ", steps=" << stats.steps 
	     << ", time=" << stats.time << '\n';	  
      }
    mdp << "END Generating ordinary propagator. ";
    mdp << "time=" << mpi.time()-time << '\n';
    end_function("generate");
  }
};

// //////////////////////////////////////////////////////////////////////
// function to print a propagator by components. It prompts for a site //
// similar to the CANOPY one. I used this to test the converter        //
/////////////////////////////////////////////////////////////////////////

void print_propagator(fermi_propagator &S) {
  begin_function("print_propagator");
  int x0,x1,x2,x3;
  site x(S.lattice());
  int spin_source, spin_sink;
  int color_source, color_sink;
  mdp_complex tmp;
  int do_exit=false;
  int nc=S.nc;
  do {
    mdp << "\nCheck point!\n";
    mdp << "Here you called the function to print the propagator\n";
    mdp << "Enter the coordinates (x0,x1,x2,x3 or 'quit' to end): ";
    if(ME==0) {
      string stringa;
      cin >> stringa;
      if(stringa=="quit") do_exit=true;
      else sscanf(stringa.c_str(),"%i,%i,%i,%i", &x0,&x1,&x2,&x3);
    }
    mpi.broadcast(do_exit,0);
    if(do_exit==true) {
      mdp << '\n';
      break;
    };
    mpi.broadcast(x0,0);
    mpi.broadcast(x1,0);
    mpi.broadcast(x2,0);
    mpi.broadcast(x3,0);
    if(on_which_process(S.lattice(),x0,x1,x2,x3)==ME) {
      x.set(x0,x1,x2,x3);
      for(color_source=0; color_source<nc; color_source++)
	for(spin_source=0; spin_source<4; spin_source++) {
	  mdp << "Source: spin=" << spin_source 
	      << ", color=" << color_source << '\n';
	  for(spin_sink=0; spin_sink<4; spin_sink++) {
	    mdp << "[ ";
	    for(color_sink=0; color_sink<nc; color_sink++) {
	      tmp=S(x,spin_sink, spin_source, color_sink, color_source);
	      mdp << tmp;
	      if(color_sink<nc-1) mdp << ",\t";
	    };
	    mdp << " ]\n";
	  };
	}; 
      fflush(stdout); 
    };
  } while(1);
  begin_function("print_propagator");
};
