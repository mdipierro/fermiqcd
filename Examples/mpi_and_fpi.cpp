#include "fermiqcd.h"       // include FermiQCD libraries
int main(int argc, char** argv) { 
  mdp.open_wormholes(argc,argv);    // START
  define_base_matrices("FERMILAB"); // set gamma matrices
  int L[]={32,8,8,8};       // lattice volume
  int n=3;                  // SU(n) gauge group
  int N=100;                // number of gauge configurations
  mdp_lattice lattice(4,L); // make a 4D lattice
  gauge_field U(lattice,n); // make a gauge field U
  coefficients gauge;       // set physical gauge parameters
  gauge["beta"]=6.0;        // beta=6/g^2 sets lattice spacing
  coefficients quark;       // set physical quark parameters
  quark["kappa"]=0.12345;   // kappa sets the quark mass
  mdp_array<float,1> C(L[0]); // make a 1D array of 32 float elements
  fermi_field psi(lattice,n); // make a fermionic field
  fermi_field chi(lattice,n); // make another fermionic field
  mdp_site x(lattice);      // variable to loop on the lattice
  set_hot(U);               // make a random gauge configuration  
  for(int t=0; t<L[0]; t++) C(t)=0;          // zero the output
  for(int k=0; k<N; k++) {                   // loop over MCMC
    WilsonGaugeAction::heatbath(U,gauge,10); // do 10 MCMC steps
    for(int alpha=0; alpha<4; alpha++)       // source spin index
      for(int i=0; i<n; i++) {	             // source color index
	forallsites(x) {                     // build source
	  psi(x)=0; if(x(0)==0) psi(x,alpha,i)=1;
	}
	mul_invQ(chi,psi,U,quark);           // invert propagator
	forallsites(x) {                     // contract sink
	  for(int beta=0; beta<4; beta++)
	    for(int j=0; j<n; j++) {	
	      // change following line for other mesons
	      C(x(0))+=real(conj(chi(x,beta,j))*chi(x,beta,j));
	    }
	}
      }						
  }
  mpi.add(C.address(),C.size());             // add parallel results
  for(int t=0; t<L[0]; t++)                  // print output 
    mdp << "C(" << t << ")=" << C(t) << endl; 
  mdp.close_wormholes();                     // STOP
  return 0;
}
