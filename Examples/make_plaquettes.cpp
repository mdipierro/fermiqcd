#include "fermiqcd.h"

int main(int argc, char** argv) {
  mdp.open_wormholes(argc,argv);    
  int nc=3;
  mdp_field_file_header header = get_info(argv[1]);
  assert(header.ndim==4);
  assert(header.box[2]==header.box[1]);
  assert(header.box[3]==header.box[1]);
  int nt = header.box[0];
  int nx = header.box[1];
  int box[]={nt,nx,nx,nx};
  mdp_lattice lattice(4,box,default_partitioning0,torus_topology,0,2,false);
  gauge_field U(lattice,nc);
  U.load(argv[1]);
  mdp_site x(lattice);
  x.set(0,0,0,0);
  cout << U(x,0) << endl;
  cout << U(x,0)*hermitian(U(x,0)) << endl;
  mdp << "plaquette:" << average_plaquette(U) << endl;
  mdp.close_wormholes();
  return 0;
}
