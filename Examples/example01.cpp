// Program: example01.cpp
#include "mdp.h"

int main(int argc, char **argv) {
  mpi.open_wormholes(argc,argv);
  int               ndim=3; 
  int               mybox[]={4,4,4}; 
  mdp_lattice       mylattice(ndim, mybox); 
  mdp_matrix_field  F(mylattice,5,5);
  site              x(mylattice);

  forallsites(x)
    F(x)=inv(exp(mylattice.random(x).SU(5)));
  F.update();
  forallsites(x) {
    F(x)=1.0/7.0*(F(x)+F(x+0)+F(x-0)+F(x+1)+F(x-1)+F(x+2)+F(x-2)); 
    cout << "F(x) for x=(" << x(0) << "," << x(1) << ")\n";
    cout << F(x) << '\n';
  }
  mpi.close_wormholes();
  return 0;
}
