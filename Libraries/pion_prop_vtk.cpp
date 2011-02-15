#include "fermiqcd.h"

int main(int argc, char** argv) {
  mdp.open_wormholes(argc,argv);

  int nc=3;
  mdp_field_file_header header = get_info(argv[1]);
  mdp_lattice spaceslice(4,header.box); // assumes 1xNxNxN
  mdp_lattice space(3,header.box+1);
  mdp_site x(spaceslice);
  mdp_site x3(space);
  fermi_propagator S(spaceslice,nc);
  mdp_field<mdp_real> prop(space);
  char output[512];
  sprintf(output,"%s.pion.vtk",argv[1]);

  forallsites(x) {
    x3.set(x(1),x(2),x(3));
    prop(x)=0;
    for(int a=0; a<S.nspin; a++) 
      for(int b=0; b<S.nspin; b++) 
	for(int i=0; i<S.nc; i++) 
	  for(int j=0; j<S.nc; j++) 
	    prop(x3) += pow(abs(S(x,a,b,i,j)),2);    
  }

  prop.save_vtk(output);

  mdp.close_wormholes();
  return 0;
}
