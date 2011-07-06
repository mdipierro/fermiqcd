//#define BLOCKSITE 4
//#define TWISTED_BOUNDARY
#include "fermiqcd.h"

void make_slices(int nt, int nx, char* filename) {
  int box[]={nt,nx,nx,nx}, nc=3;
  generic_lattice lattice(4,box,
			  default_partitioning0,
			  torus_topology,
			  0, 1,false);  
  gauge_field U(lattice,nc);
  mdp << "success in allocating vector\n";
  char filename2[128];
  U.load(filename);
  sprintf(filename2,"%s.topological_charge_20.vtk",filename);
  float tc=topological_charge_vtk(U,filename2,-1);
  mdp << "topological_charge =" << tc << endl;
}

int main(int argc, char** argv) {    
  mpi.open_wormholes(argc, argv);   
  define_base_matrices("FERMILAB");
  mdp_field_file_header header = get_info(argv[1]);
  assert(header.ndim==4);
  assert(header.box[2]==header.box[1]);
  assert(header.box[3]==header.box[1]);
  int nt = header.box[0];
  int nx = header.box[1];  
  make_slices(nt,nx,argv[1]);
  mpi.close_wormholes();             
  return 0;                         
}
                                           
