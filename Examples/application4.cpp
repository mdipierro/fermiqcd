// program: application4.cpp
#include "mdp.h"
 
int main(int argc, char **argv) {
  mdp.open_wormholes(argc, argv);
  int box[4]={8,12,12,12};
  mdp_lattice space(4,box);
  mdp_matrix_field U(space,3,3);
  mdp_matrix_field V(space,3,3);
  mdp_site x(space);
  forallsites(x) {
    V(x)=space.random(x).SU(3);
    U(x)=0;
  };
  U.update();
  V.update();
  for(int i=0; i<10000; i++) {
    forallsites(x)
      U(x)=0.125*(cos(U(x)+V(x))+
                  U(x+0)+U(x-0)+
                  U(x+1)+U(x-1)+
                  U(x+2)+U(x-2)+
                  U(x+3)+U(x-3));
    /* uncomment to print convergence
    if(onwhichprocess(space,0,0,0,0)==ME) {
      x.set(0,0,0,0);
      cout << real(U(x,0,0)) << '\n';
    }
    */
    U.update();
  };
  V.save("V_field.dat");
  U.save("U_field.dat");
  mdp.close_wormholes();
  return 0;
}; 
