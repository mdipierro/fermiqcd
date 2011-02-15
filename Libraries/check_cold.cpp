#include "fermiqcd.h"

int main(int argc, char** argv) {
  mdp.open_wormholes(argc,argv);    
  int nc=3;
  int nt,nx;
  sscanf(argv[1],"%ix%i",&nt,&nx);
  int box[]={nt,nx,nx,nx};
  mdp_lattice lattice(4,box,default_partitioning0,torus_topology,0,2,false);
  gauge_field U(lattice,nc);
  U.load(argv[2]);
  mdp_site x(lattice);
  x.set(0,0,0,0);
  cout << U(x,0) << endl;
  forallsites(x) for(int mu=0;mu<4; mu++)
      if(max(U(x,mu)-1)>0) {
	  cout << x << mu << "\n" << U(x,mu) << endl;
	  break;
      }
  mdp << "plaquette:" << average_plaquette(U) << endl;
  mdp.close_wormholes();
  return 0;
}
