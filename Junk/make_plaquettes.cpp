#include "fermiqcd.h"

mdp_real average_plaquette1(gauge_field &U,int mu, int nu) {
  mdp_lattice& lattice=U.lattice();
  mdp_site x(lattice);
  mdp_real sum=0.0;
  forallsites(x) 
    sum+=real(trace(U(x,mu)*U(x+mu,nu)*
		    hermitian(U(x,nu)*U(x+nu,mu))));
  mdp.add(sum);
  return sum/lattice.size();
}

mdp_real average_plaquette2(gauge_field &U,int mu, int nu) {
  mdp_lattice& lattice=U.lattice();
  mdp_site x(lattice);
  mdp_complex sum=0.0;
  int path[4][2]={{+1,mu},{+1,nu},{-1,mu},{-1,nu}};
  sum=average_path(U,4,path);
  return real(sum);
}

int main(int argc, char** argv) {
  mdp.open_wormholes(argc,argv);    
  int nc=3;
  int box[]={8,4,4,4};
  mdp_lattice lattice(4,box);
  gauge_field U(lattice,nc);
  int mu=0, nu=1;
  U.load(argv[1]);
  mdp << "plaquette:" << average_plaquette(U,mu,nu) << endl;
  mdp << "plaquette:" << average_plaquette1(U,mu,nu) << endl;
  mdp << "plaquette:" << average_plaquette2(U,mu,nu) << endl;
  mdp.close_wormholes();
  return 0;
}
