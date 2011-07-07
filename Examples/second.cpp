#include "fermiqcd.h"

int main(int argc, char** argv) {
  mdp.open_wormholes(argc,argv); // mpirun 

  int L[]={10,10};
  mdp_lattice lattice(2,L);
  mdp_matrix_field psi(lattice,3,3);
  mdp_site x(lattice);
  
  forallsites(x) psi(x)=lattice.random(x).SU(3);

  psi.update();

  forallsites(x) psi(x) = (psi(x+0)+psi(x-0)+psi(x+1)+psi(x-1))/4;

  psi.update();
  

  mdp << "hello world\n";

  mdp.close_wormholes();
  return 0;
}
