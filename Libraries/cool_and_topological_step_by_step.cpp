//#define BLOCKSITE 4
//#define TWISTED_BOUNDARY
#include "fermiqcd.h"

void test_gauge(int nt, int nx, char* filename) {
  int box[]={nt,nx,nx,nx}, nc=3;
  generic_lattice lattice(4,box,
			  default_partitioning0,
			  torus_topology,
			  0, 1,false);
  gauge_field U(lattice,nc);
  char filename2[200];
  U.load(filename);
  // U.switch_endianess_4bytes();
  for(int k=0; k<=20; k+=5) {
      sprintf(filename2,"%s.topological_charge_%i.vtk",filename,k);
      float tc=topological_charge_vtk(U,filename2,0);
      mdp << "topological_charge=" << tc << endl;
      ApeSmearing::smear(U,0.7,5,10);
  }
}

int main(int argc, char** argv) {    
  mpi.open_wormholes(argc, argv);   
  define_base_matrices("FERMILAB");
  int nt, nx;
  sscanf(argv[1],"%ix%i",&nt,&nx);
  test_gauge(nt,nx,argv[2]);
  mpi.close_wormholes();             
  return 0;                         
}
                                           
