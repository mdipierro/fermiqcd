#include "fermiqcd.h"

int main(int argc, char** argv) {
  mdp.open_wormholes(argc,argv);    
  int nc=3;
  int box[]={8,4,4,4};
  mdp_lattice lattice(4,box);
  gauge_field U(lattice,nc);
  set_cold(U);
  U.save("gauge.cold");  
  mdp.close_wormholes();
  return 0;
}
