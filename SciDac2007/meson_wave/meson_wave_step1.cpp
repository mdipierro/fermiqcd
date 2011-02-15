// #define SSE2
#include "fermiqcd.h"
#include "dump.h"
#include "fermiqcd_topological_charge.hpp"
/* 
   this program compute the wave function of a heavy-light meson by 
   computing s(x3) = |<B(x) q(y) \bar q (y) \bar B(z)>|^2
   where B(x)= h(x) \gamma_5 \bar q(x)
   and x(0)=0, y(0)=T/4, z(0)=T/2, x(1..3)=z(...3)=(0,0,0), and y(1..3)=x3
*/

int main(int argc, char** argv) {
  mdp.open_wormholes(argc,argv);
  define_base_matrices("FERMILAB");
  int L[4]={20,10,10,10};
  mdp_lattice lattice(4,L);
  mdp_site x(lattice);
  int nc=2;
  gauge_field U(lattice,nc);
  fermi_field chi(lattice,nc);
  fermi_field psi(lattice,nc);
  fermi_field phi(lattice,nc);
  mdp_lattice space(3,L+1);
  mdp_site x0(space);
  mdp_site x3(space);
  mdp_matrix A(nc,nc);
  mdp_field<mdp_complex> s(space,16);
  coefficients gauge; gauge["beta"]=2.3;
  // coefficients gauge; gauge["beta"]=5.5;
  coefficients quark; quark["kappa"]=0.130;
  int L8=8;
  char filename[128]; 
  // default_fermi_action=FermiCloverActionSSE2::mul_Q;

  set_cold(U);
#ifndef TEST
  for(int k=0; k<2000; k++) {
    cout << k << endl;
    WilsonGaugeAction::heatbath(U,gauge,10);
    U.save("gauge_20x10x10x10.mdp");
  }
#endif

  forallsites(x3)
    for(int b=0; b<4; b++)
      for(int c=0; c<4; c++)
	for(int i=0; i<nc; i++)
	  s(x3,4*c+b)=0;

#ifdef TEST  
  for(int k=0; k<1; k++) {
    cout << k << endl;
    WilsonGaugeAction::heatbath(U,gauge,1);
    U.save("gauge_20x10x10x10.mdp");
#else
  for(int k=0; k<1000; k++) {
    cout << k << endl;
    WilsonGaugeAction::heatbath(U,gauge,50);
    U.save("gauge_20x10x10x10.mdp");
#endif
    A=1;
    x0.set(L[1]/2,L[2]/2,L[3]/2);
    for(int t=0; t<L8; t++) {
      x.set(t,x0(0),x0(1),x0(2));
      A=A*U(x,0); // DOES NOT WORK IN PARALLEL!!!
    }
    for(int a=0; a<2; a++) {
      for(int j=0; j<nc; j++) {
	// make source x=(0,x3)
	chi=0; 
	x.set(0,x0(0),x0(1),x0(2)); chi(x,a,j)=1;
	chi.update();
	mul_invQ(psi,chi,U,quark,1e-7,1e-5,200);
	// make other source at x=(L[0]/2,x3) by mul by heavy quark propagator
	chi=0;
	x.set(L8,x0(0),x0(1),x0(2)); 
	for(int k=0; k<nc; k++) chi(x,Gamma5_idx[a],k)=Gamma5_val[a]*A(j,k);

	chi.update();
	mul_invQ(phi,chi,U,quark,1e-7,1e-5,200);

	forallsites(x3) {
	  x.set(L8/2,x3(0),x3(1),x3(2));
	  for(int b=0; b<4; b++)
	    for(int c=0; c<4; c++)
	      for(int i=0; i<nc; i++)
		s(x3,4*c+b)+=phi(x,c,i)*conj(psi(x,b,i));
	}
      }
    }
    sprintf(filename,"wave.%.3i.mdp",k);
    s.save(filename);
  }
  mdp.close_wormholes();
  return 0;
}
