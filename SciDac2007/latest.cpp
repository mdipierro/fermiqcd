/*
    python qcd.py code +hot +heatbath steps=100 +loop n=100 { +ape_smear steps=1 +topological_charge }

    +hot warning: assuming default argument TxXxYxZ=16x4x4x4
    +hot warning: assuming default argument nc=3
    +heatbath warning: assuming default argument beta=5.0
    +ape_smear warning: assuming default argument alpha=0.7
    +ape_smear warning: assuming default argument cooling_steps=10
    +topological_charge warning: assuming default argument filename=topological_charge_*.vtk
    +topological_charge warning: assuming default argument t=-1
*/
#include "fermiqcd.h"

int main(int argc, char** argv) {
   mdp.open_wormholes(argc,argv);
   string filename;
   coefficients coeff;
   int L[]={16,4,4,4};
   mdp_lattice spacetime(4,L);
   int nc=3;
   gauge_field U(spacetime,nc);
   set_hot(U);

   coeff["beta"]=5.0;
   WilsonGaugeAction::heatbath(U,coeff,100);
   for(int i0=0; i0<100; i0++) {

      ApeSmearing::smear(U,0.7,1,10);
      {float tc=topological_charge_vtk(U,"topological_charge_*.vtk",-1);
      mdp << "topological_charge=" << tc << endl; }
   }


   mdp.close_wormholes();
   return 0;
}
