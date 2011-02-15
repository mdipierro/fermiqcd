//#define BLOCKSITE 4
//#define TWISTED_BOUNDARY
#include "fermiqcd.h"

void test_gauge() {
  int box[]={4,4,4,4}, nc=3;
  generic_lattice lattice(4,box);
  gauge_field U(lattice,nc);
  set_hot(U);
  cout << average_plaquette(U,0,1) << endl;
  int path[][2]={{+1,0},{+1,1},{-1,0},{-1,1}};
  cout << average_path(U,4,path) << endl;
}

int main(int argc, char** argv) {    
  mpi.open_wormholes(argc, argv);   
  define_base_matrices("FERMILAB");    
  test_gauge();
  mpi.close_wormholes();             
  return 0;                         
}
                                           
