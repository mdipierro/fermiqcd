#include "fermiqcd.h"

int main(int argc, char** argv) {
  mdp.open_wormholes(argc,argv); // mpirun 

  float beta;
  for(beta=5; beta<10; btea++)
  for(int T=1; T<11; T++) { 
    int L[]={T,10,10,10};
    mdp_lattice lattice(4,L);
    mdp_lattice lattice3d(3,L+1);
    gauge_field U(lattice,3);
    mdp_matrix_field V(lattice3d,3,3);
    mdp_field<mdp_real> q(lattice3d,2);
    
    mdp_site x(lattice);
    mdp_site x3d(lattice3d);
    coefficients gauge;
    gauge["beta"]=beta;
    
    int k,mu=0,nu=1;
    mdp_complex s=0;
    
    forallsites(x) 
      for(int mu=0; mu<4; mu++)
	U(x,mu)=lattice.random(x).SU(3);
    
    forallsites(x3d)
      V(x3d)=1;
    
    char filename[128];
    sprintf(filename,"gauge_1000_%ix%i_%f.mdp",L[0],L[1],gauge["beta"]);
    
    // if file filename does not exists {
    if(true) {
	U.update();
	for(k=0; k<100; k++)
	  WilsonGaugeAction::heatbath(U,gauge,10);
	U.save(filename);
    }
    U.load(filename);

    for(k=0; k<100; k++) {
      WilsonGaugeAction::heatbath(U,gauge,1);
      
      
      for(int t=0; t<L[0]; t++) {
	forallsites(x3d) {
	  x.set(t,x3d(0),x3d(1),x3d(2));
	  V(x3d)=V(x3d)*U(x,0);
	}
      }
      
      forallsites(x3d) {
	mdp_complex z=trace(V(x3d));
	q(x3d,0)=real(z);
	q(x3d,1)=imag(z);
      }
      /// FIX FILE NAMES depends on T, beta, and K
      q.save_vtk("real_polyakov_line.vtk",-1,0,0,false);
      q.save_vtk("imag_polyakov_line.vtk",-1,1,0,false);
      return 0;
    }
  }
  mdp.close_wormholes();
  return 0;
}

