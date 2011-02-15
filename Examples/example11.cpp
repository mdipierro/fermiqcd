// Program: example11.C
#include "mdp.h" 

int main(int argc, char **argv) {
  mpi.open_wormholes(argc,argv);
  int mybox[]={8,8};
  generic_lattice mylattice(2,mybox);
  mpi.close_wormholes();
}
