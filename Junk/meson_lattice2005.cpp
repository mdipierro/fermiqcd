#include "fermiqcd.h"
int main(int argc, char **argv) {
  mdp.open_wormholes(argc, argv);   // START
  define_base_matrices("FERMIQCD"); // set Gamma convention
  int n=3, nconfig=100;
  int L[] = {16,8,8,8};
  mdp_lattice lattice(4,L);         // declare lattice
  gauge_field U(lattice, n);        // declare fields
  fermi_propagator S(lattice,n);    // declare propagator
  mdp_site x(lattice);              // declare a site var
  coefficients gauge; gauge["beta"]=6.0; // set parameters
  coefficients quark; quark["kappa"]=0.1234; quark["c_{SW}"]=0.0;
  default_fermi_action=FermiCloverActionFast::mul_Q;
  mdp_array<float,1> Cpi(L[TIME]);    // declare and zero Cpi
  for(int t=0; t<L[TIME]; t++) Cpi(t)=0;
  set_hot(U);
  for(int k=0; k < nconfig; k++) { 
    WilsonGaugeAction::heatbath(U,gauge); // do heatbath
    U.save(string("gauge")+tostring(k));  // save config
    if(quark["c_{SW}"]!=0) compute_em_field(U);
    generate(S,U,quark,1e-20,1e-12);      // make propagator
    forallsites(x)                        // contract pion
      for(int a=0; a<4; a++)              // source spin
	for(int b=0; b<4; b++)            // sink spin
	  Cpi[x(TIME)]+=real(trace(S(x,a,b)*hermitian(S(x,a,b))));
    mpi.add(Cpi.address(),Cpi.size());    // parallel add
    for(int t=0; t<L[TIME]; t++)    
      mdp << t << " " << Cpi(t) << endl;  // print output
  } 
  mdp.close_wormholes();            // STOP
  return 0;
} 
