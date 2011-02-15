#include "fermiqcd.h"

int main(int argc, char** argv) {
  mdp.open_wormholes(argc,argv);
  int L[]={32,24,24,24};
  mdp_lattice lattice(4,L,default_partitioning0,torus_topology,1,0,false);
  gauge_field U(lattice,3);
  U.load(argv[1]);
  mdp_site x(lattice),y(lattice);

  x.set(0,0,0,0);
  y.set(31,23,23,23);
  for(int mu=0; mu<4; mu++) U(y,mu)=U(x,mu);
  cout << average_plaquette(U) << endl;
  
  mdp.close_wormholes();
  return 0;
}

