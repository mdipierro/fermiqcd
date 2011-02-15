#include "fermiqcd.h"

int main(int argc, char** argv) {    
  mdp.open_wormholes(argc, argv);   
  
  int box[]={4,4,4,4}, nc=3;
  generic_lattice lattice(4,box);
  gauge_field U(lattice,nc);
  coefficients coeff;
  coeff["beta"]=6.0;
  set_hot(U);
  mdp << "average plaquette = " << average_plaquette(U) << endl;
  mdp.close_wormholes();             
  
  return 0;                         
}
                                           
