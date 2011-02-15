// Program: example08.C
#include "mdp.h"

int main(int argc, char **argv) {
   mpi.open_wormholes(argc,argv);
   if(mpi.nproc()==2) { 
     if(ME==1) {
       int b;
       b=4*8;
       mpi.put(b,0);
     }
     if(ME==0) {
       int a,b;
       a=5*7;
       mpi.get(b,1);
       mpi << a+b << '\n';
     }
   } else {
     mpi << "Sorry, this only runs on 2 processes\n";
   }
   mpi.close_wormholes();
}
