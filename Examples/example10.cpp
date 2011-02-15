// Program: example13.C
#include "mdp.h" 

int main(int argc, char **argv) {
   mpi.open_wormholes(argc,argv);
   int i,j;
   for(i=0; i<5; i++) {
     if(ME==0) j=i;
     else      j=0;
     if(i%2==0) mpi.barrier();
     mpi.broadcast(j,0);
     cout << "I am process " << ME 
	  << ", i=" << i 
	  << ", j=" << j << '\n';
   }
  mpi.close_wormholes();
}
