#include "fermiqcd.h"

int main(int argc, char** argv) {
  mdp.open_wormholes(argc,argv); // mpirun 

  int L[]={6,6,6,6};
  mdp_lattice lattice(4,L);
  gauge_field U(lattice,3);
  mdp_site x(lattice);
  
  forallsites(x) 
   for(int mu=0; mu<4; mu++)
     U(x,mu)=lattice.random(x).SU(3);

  U.update();
  for(k=0; k<1000; k++) {
    /*
    forallsites(x) 
      for(int mu=0; mu<4; mu++)
	U(x,mu) = (U(x+0,mu)+U(x-0,mu)+U(x+1,mu)+U(x-1,mu))/4;
    */
    U.save('somename');
  }

  U.update();
  

  mdp << "hello world\n";

  mdp.close_wormholes();
  return 0;
}
