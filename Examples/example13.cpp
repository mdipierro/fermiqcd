// Program: example13.C
#include "mdp.h" 

int main(int argc, char **argv) {
  mpi.open_wormholes(argc,argv);
  int mybox[]={10,10};
  mdp_lattice mylattice(2,mybox);
  mdp_site x(mylattice);
  int mu=0;
  if(ME==0) {
    x.set(0,0);
    do {
      cout << "x=(" << x(0) << "," << x(1) << ")\n";
      if(x.is_in_boundary()) error("I found the boundary");
      x=x+mu;
    } while(true);
  }
  mpi.close_wormholes();
}
