#include "fermiqcd.h"

int main(int argc, char** argv) {
  define_base_matrices("FERMILAB");
  mdp.open_wormholes(argc,argv);    
  int nc=3;
  int L[]={10,40,10,10};
  mdp_lattice lattice(4,L,default_partitioning0,torus_topology,0,2,false);
  gauge_field U(lattice,nc);
  char filename[128];
  mdp_site x(lattice);
  set_cold(U);
  vector<mdp_real> p(4);
  for(int mu=0; mu<4; mu++) p[mu]=L[mu]/2;

  for(int k=0; k<40; k++) {
    p[1]=10.0+0.25*k;
    Instanton4D A1(3,0,1, 1.0,3,p);
    p[1]=30.0-0.25*k;
    Instanton4D A2(3,0,1, 1.0,3,p);
    
    forallsites(x)
      for(int mu=0; mu<U.ndim; mu++)
	U(x,mu)=exp(-I*(A1(x,mu)+A2(x,mu)));
    
    mdp << "k=" << k << endl;
    //if(k>0) ApeSmearing::smear(U,0.7,1,10);
    sprintf(filename,"%s.%i.vtk",argv[0],k);
    mdp << "top=" << topological_charge_vtk(U,filename) << endl;    
  }
  mdp.close_wormholes();
  return 0;
}
