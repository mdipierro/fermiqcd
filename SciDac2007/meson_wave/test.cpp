#include "fermiqcd.h"

int main(int argc, char** argv) {
  mdp.open_wormholes(argc,argv);
  define_base_matrices("FERMILAB");
  cout << Gamma[0] << endl;
  cout << Gamma5 << endl;
  cout << (1+Gamma[0])/2*Gamma5 << endl;
  return 0;
}
