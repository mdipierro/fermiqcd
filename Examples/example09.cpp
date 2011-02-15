// Program: example09.C
#include "mdp.h" 

int main(int argc, char **argv) {
   mpi.open_wormholes(argc,argv);
   int a;
   if(mpi.nproc()==2) {
     if(ME==1) {
       a=4*8;
       mpi.add(a);
     }
     if(ME==0) {
       a=5*7;
       mpi.add(a);
       cout << "a=" << a << '\n';
     }
   } else {
     mpi << "Sorry, this only runs on 2 processes\n";
   }
   mpi.close_wormholes();
}
