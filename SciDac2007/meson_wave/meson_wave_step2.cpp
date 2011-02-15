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
  mdp_lattice space(3,L+1);
  mdp_site x(space);
  mdp_site y(space);
  mdp_field<mdp_complex> s(space,16);
  mdp_field<mdp_complex> sum(space);
  mdp_field<float> Q(space);
  char filename[128]; 

  forallsites(x)
    sum(x)=0;

  mdp_complex z;
  for(int k=0; k<1000; k++) {
    cout << k << endl;
    sprintf(filename,"wave.%.3i.mdp",k);
    s.load(filename);
    forallsites(x) {
      z=0;
      for(int c=0; c<4; c++)
	for(int b=0; b<4; b++)
	  z+=s(x,4*c+b)*(Gamma5)(c,b);      
      sum(x)+=z;
      
      y.set(L[1]-x(0),x(1),x(2)); sum(y)+=z;      
      y.set(x(0),L[2]-x(1),x(2)); sum(y)+=z;
      y.set(x(0),x(1),L[3]-x(2)); sum(y)+=z;
      y.set(L[1]-x(0),L[2]-x(1),x(2)); sum(y)+=z;
      y.set(L[1]-x(0),x(1),L[3]-x(2)); sum(y)+=z;
      y.set(x(0),L[2]-x(1),L[3]-x(2)); sum(y)+=z;
      y.set(L[1]-x(0),L[2]-x(1),L[3]-x(2)); sum(y)+=z;
      
    }
    forallsites(x) Q(x)=pow(abs(sum(x)),2)/(k+1)+1e-12;
    sprintf(filename,"meson_wave_plus.%.3i.vtk",k);
    dump(Q,0,filename);
  }
  mdp.close_wormholes();
  return 0;
}
