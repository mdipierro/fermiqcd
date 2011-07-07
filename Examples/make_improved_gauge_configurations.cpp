#include "fermiqcd.h"

int main(int argc, char** argv) {
  mdp.open_wormholes(argc,argv);    
  int nc=3;
  int box[]={8,4,4,4};
  mdp_lattice lattice(4,box);
  gauge_field U(lattice,nc);
  coefficients gauge;
  gauge["beta"]=6.0;
  gauge["zeta"]=1.0; // has to be one
  gauge["u_s"]=1.0;
  gauge["u_t"]=1.0;
  set_cold(U);
  for(int k=0; k<100; k++) {
    ImprovedGaugeAction::heatbath(U,gauge,10,"MILC");
    mdp << "plaquette:" << average_plaquette(U) << endl;
    U.save(string("gauge.")+tostring(k));  
  }
  mdp.close_wormholes();
  return 0;
}
