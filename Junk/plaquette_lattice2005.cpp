#include "fermiqcd.h"
int main(int argc, char **argv) {
  mdp.open_wormholes(argc, argv); // START
  int n=4, nconfig=100;           
  int L[] = {8,8,8,8,8};          
  mdp_lattice lattice(5,L);       // declare lattice
  gauge_field U(lattice, n);      // declare fields
  coefficients gauge; gauge["beta"]=6.0; // set parameters
  set_cold(U);
  for(int k=0; k < nconfig; k++) { 
     WilsonGaugeAction::heatbath(U,gauge,10); // do heathbath
     U.save(string("gauge")+tostring(k));     // save gauge config
     mdp << average_plaquette(U) << endl;     // print plaquette
  } 
  mdp.close_wormholes();         // STOP
  return 0;
} 
