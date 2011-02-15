#include "fermiqcd.h"
#include "mdp_all.h"
#include "dump.h"
#include "fermiqcd_topological_charge.hpp"

int main(int argc, char** argv) {
  mdp.open_wormholes(argc,argv);
  define_base_matrices("FERMILAB");
  int nc=3;
  int L[]={12,12,12,12};
  mdp_lattice lattice(4,L);
  mdp_lattice cube(3,L+1);
  gauge_field U(lattice,nc);
  mdp_complex_field F(cube,12);
  mdp_field<float> Q3(cube);
  mdp_site x(lattice), y(lattice);
  mdp_site x3(cube);
  mdp_matrix A,B,C; 
  mdp_complex d=0.0;
  char name[100];
  coefficients gauge; gauge["beta"]=5.0;
  /*
  set_cold(U);
  for(int k=0; k<100; k++) { 
    cout << k << endl;
    WilsonGaugeAction::heatbath(U,gauge,1);
  }
  */
  U.load("sample_gauge_10x10x10x10.mdp");
  WilsonGaugeAction::heatbath(U,gauge,100);
  for(int conf=0; conf<100; conf++) {
    WilsonGaugeAction::heatbath(U,gauge,10);
    ape_smearing(U,0.7,5,10);
    compute_em_field(U);

    x.set(0,L[1]/4,L[2]/2,L[3]/2);
    A=U(x,0); for(int k=1; k<L[0]; k++) { x=x+0; A=A*U(x,0); }   
    x.set(0,L[1]*3/4,L[2]/2,L[3]/2);
    B=U(x,0); for(int k=1; k<L[0]; k++) { x=x+0; B=B*U(x,0); }   
    cout << "1" << endl;
    C=A*hermitian(B);
    d+=L[0]*trace(C); // quark and antiquark
    forallsites(x) {
      x3.set(x(1),x(2),x(3));
      F(x3,0)+=trace(C*plaquette(U,x,0,1));
      F(x3,1)+=trace(C*plaquette(U,x,0,2));
      F(x3,2)+=trace(C*plaquette(U,x,0,3));
      F(x3,3)+=trace(C*plaquette(U,x,1,2));
      F(x3,4)+=trace(C*plaquette(U,x,1,3));
      F(x3,5)+=trace(C*plaquette(U,x,2,3));
      F(x3,6)+=trace(plaquette(U,x,0,1));
      F(x3,7)+=trace(plaquette(U,x,0,2));
      F(x3,8)+=trace(plaquette(U,x,0,3));
      F(x3,9)+=trace(plaquette(U,x,1,2));
      F(x3,10)+=trace(plaquette(U,x,1,3));
      F(x3,11)+=trace(plaquette(U,x,2,3));
    }
    /*
    for(int k=0; k<6; k++) {
      for(int mu=0; mu<12; mu++) {
	forallsitesofparity(x3,k%2) F(x3,mu)=0.8*F(x3,mu)+0.2/6*(F(x3+0,mu)+
								 F(x3-0,mu)+
								 F(x3+1,mu)+
								 F(x3-1,mu)+
								 F(x3+2,mu)+
								 F(x3-2,mu));
      }
    }
    */
    cout << "projecting to 3D...\n";
    forallsites(x3) 
      Q3(x3)=abs((F(x3,0)+F(x3,1)+F(x3,2)-F(x3,3)-F(x3,4)-F(x3,5))/d-
		 (F(x3,6)+F(x3,7)+F(x3,8)-F(x3,9)-F(x3,10)-F(x3,11))/L[0]/(conf+1));
    cout << "saving vtk file\n";
    sprintf(name,"qqbar_%i.vtk", conf);
    dump(Q3,name);
  }
  // compute correlatiton between three (or four) polyakov lines.
  // average over gauge configurations.
  
  mdp.close_wormholes();
  return 0;
}
