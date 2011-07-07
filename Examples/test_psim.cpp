#include "fermiqcd.h"

int main(int argc, char** argv) {
  mdp.open_wormholes(argc, argv);
  int L[]={256,32,32,32};

  mdp_lattice lattice(4,L,default_partitioning0,torus_topology,0,1,false);
  mdp << "here\n";

  gauge_field U(lattice,3);
  mdp << "there\n";

  U.update();
  mdp << "there2\n";

  gauge_field V(lattice,3);
  mdp << "there3\n";

  mdp.close_wormholes();
  return 0;
}
