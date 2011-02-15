/*
    python qcd.py code +hot TxXxYxZ=4x4x4x4 +loop n=10 { +heatbath +plaquette } +loop n=1 { +ape_smear steps=1 +plaquette +topological_charge_vtk }

    +hot warning: assuming default argument nc=3
    +heatbath warning: assuming default argument beta=5.0
    +heatbath warning: assuming default argument steps=1
    +ape_smear warning: assuming default argument alpha=0.7
    +ape_smear warning: assuming default argument cooling_steps=10
    +topological_charge_vtk warning: assuming default argument prefix=topological_charge
    +topological_charge_vtk warning: assuming default argument t=-1
*/
#include "fermiqcd.h"

int main(int argc, char** argv) {
   mdp.open_wormholes(argc,argv);
   string filename;
   coefficients coeff;
   int L[]={4,4,4,4};
   mdp_lattice spacetime(4,L);
   int nc=3;
   gauge_field U(spacetime,nc);
   set_hot(U);

   for(int i0=0; i0<10; i0++) {
      coeff["beta"]=5.0;
      WilsonGaugeAction::heatbath(U,coeff,1);
      mdp << "average_plaquette=" << average_plaquette(U) << endl;
   }

   for(int i0=0; i0<10; i0++) {
      ApeSmearing::smear(U,0.7,1,10);
      U.save("test*");
      mdp << "average_plaquette=" << average_plaquette(U) << endl;
      {float tc=topological_charge_vtk(U,"topological_charge*",-1);
      mdp << "total topological charge=" << tc << endl; }
   }


   mdp.close_wormholes();
   return 0;
}
