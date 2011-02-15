// Program: example12.C
#include "mdp.h" 

int main(int argc, char **argv) {
    mpi.open_wormholes(argc,argv);
    int mybox[]={10,10};
    generic_lattice mylattice(2,mybox);
    site x(mylattice);
    forallsites(x)
      cout << "Site=(" << x(0) << "," << x(1) 
	   << ") is stored by " << ME << '\n';
    mpi.close_wormholes();
}
