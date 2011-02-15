#include "fermiqcd.h"        
 
int main(int argc, char **argv) {
  mdp.open_wormholes(argc, argv);      
 
  int nc=4, Nconfig=5;
  int mybox[] = {2,2,2,2,2};         
  generic_lattice mylattice(5, mybox); 
  gauge_field     U(mylattice, nc);    
  coefficients coeff;
  coeff["beta"]=2.0;         

  set_hot(U);                          
  
  for(int config=0; config < Nconfig; config++) {
      WilsonGaugeAction::heatbath(U,coeff,1);
      mdp << "config=" << config 
          << "plaquette=" 
          << average_plaquette(U) << endl;
  }
  mdp.close_wormholes();     
  return 0;
}

