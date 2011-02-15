#include "mdp.h"

int main(int argc, char** argv) {
  mdp.open_wormholes(argc,argv);
  int L[]={100};
  mdp_lattice line(1,L);
  mdp_field<int> spin(line);
  mdp_site x(line);
  int dE=0, H=L[0], dH=0;
  float kappa=2.0;
  forallsites(x) spin(x)=+1;
  while(1) {
    dH=0;
    for(int parity=0; parity<2; parity++) {
      forallsitesofparity(x,parity) {
	dE=2*spin(x)*(spin(x-0)+spin(x+0));
	if(exp(-kappa*dE)>mdp_random.plain()) 
	  { spin(x)*=-1; dH=dH+2*spin(x); }
      }
      spin.update();
    }
      mpi.add(dH);
      H=H+dH;
    mdp << "magnetization=" << H << endl;
  }
  mdp.close_wormholes();
  return 1;
}
