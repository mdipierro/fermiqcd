//#define BLOCKSITE 4
//#define TWISTED_BOUNDARY
#include "fermiqcd.h"

void meson_prop_from_slice(int nt, int nx, char* filename) {
  int box[]={nt,nx,nx,nx}, nc=3;
  mdp_lattice lattice(4,box,
		      default_partitioning0,
		      torus_topology,
		      0, 1,false);  
  mdp_lattice space(3,box+1,
		      default_partitioning0,
		      torus_topology,
		      0, 1,false);  
  mdp_site x(lattice);
  mdp_site y(space);
  fermi_propagator S(lattice,nc);
  mdp_field<float> Q(space);
  mdp << "success in allocating vector\n";
  char filename2[128];
  S.load(filename);
  Matrix A=Gamma5*Gamma[1];
  mdp_complex w=0;
  forallsites(x) {
    y.set(x(1),x(2),x(3));
    Q(y)=0;
    for(int a=0; a<S.nspin; a++)
      for(int b=0; b<S.nspin; b++)
	for(int c=0; c<S.nspin; c++)
	  for(int d=0; d<S.nspin; d++) {
	    w=A(b,c)*A(d,a);
	    if(w!=0) Q(y)+=real(w*trace(S(x,a,b)*hermitian(S(x,d,c))));
	  }
  }
  sprintf(filename2,"%s.meson.vtk",filename);
  Q.save_vtk(filename2);
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
  meson_prop_from_slice(nt,nx,argv[1]);
  mpi.close_wormholes();             
  return 0;                         
}
                                           
