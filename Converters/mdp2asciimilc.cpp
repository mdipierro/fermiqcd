#include "fermiqcd.h"

int main(int argc, char** argv) {
  mdp.open_wormholes(argc,argv);
  if(argc<3) {
    mdp << "usage:\n\t" << argv[0] << " input output\n\n";
    exit(1);
  }

  if(!is_file(argv[1])) {
    mdp << "Unbale to open gauge configuration\n";
    exit(1);
  }

  FILE *outfp=fopen(argv[2],"w");    
  
  // read the lattice size from the gauge configuration itself

  mdp_field_file_header header=get_info(argv[1]);
  int* L=header.box;
  mdp << "Lattice size: " 
      << L[0] << "x" << L[1] << "x" 
      << L[2] << "x" << L[3] << endl;

  int n=3; // number of colors;
  mdp_lattice lattice(4,L);
  mdp_site x(lattice);
  if(header.bytes_per_site==4*9*2*sizeof(mdp_real)) {
    
    gauge_field U(lattice,n);  
    U.load(argv[1]);

    int x0,x1,x2,x3,i,j,mu;

    // WHAT IS THIS
    // printf("%i\n%s %s %s %s %s\n%i %i %i %i\n",a,s[0],s[1],s[2],s[3],s[4], b[0],b[1],b[2], b[3]);
    // ???
    for(x0=0; x0<L[0]; x0++) 
      for(x3=0; x3<L[3]; x3++) 
	for(x2=0; x2<L[2]; x2++) 
	  for(x1=0; x1<L[1]; x1++) {
	    x.set(x0,x1,x2,x3);
	    for(mu=1; mu<=4; mu++) {
	      for(i=0; i<3; i++)
		for(j=0; j<3; j++) {
		  fprintf(outfp, "%f %f\n", 
			 real(U(x,mu % 4,i,j)), imag(U(x,mu % 4,i,j)));
		}
	    }
	  }
  }
  fclose(outfp);

  mdp.close_wormholes();
  return 0;
}
	    











