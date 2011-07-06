#include "fermiqcd.h"

int main(int argc, char** argv) {
  mdp.open_wormholes(argc,argv); // mpirun 

  int L[]={4,10,10,10};
  mdp_lattice lattice(4,L);
  gauge_field U(lattice,3);
  coefficients gauge;
  gauge["beta"]=5.0;
  // mdp_field<float> Q(lattice);
    

  set_cold(U);
  
  for(int k=0; k<10; k++) {
    cout << k << endl;
    WilsonGaugeAction::heatbath(U,gauge,1);
  }
  topological_charge_vtk(U,"top_charge_simple.vtk");
    

  mdp.close_wormholes();
  return 0;
}

