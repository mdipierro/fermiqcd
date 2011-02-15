// Program: example14.C
#include "mdp.h"

struct mystruct { 
    float value; /* or any other structure */
};

int main(int argc, char **argv) {
    mpi.open_wormholes(argc,argv);
    int mybox[]={10,10};
    mdp_lattice         mylattice(2,mybox);
    mdp_field<mystruct> myfield(mylattice);
    mdp_site            x(mylattice);
    forallsites(x)
      myfield(x).value=0;
    myfield.update();
    mpi.close_wormholes();
}
