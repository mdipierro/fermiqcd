#include "fermiqcd.h"

int main(int argc, char** argv) {
  mdp.open_wormholes(argc,argv);    
  define_base_matrices("FERMILAB");
  int nc=3;
  int box[]={8,4,4,4};
  mdp_lattice lattice(4,box);
  gauge_field U(lattice,nc);
  fermi_propagator S(lattice,nc);
  mdp_site x(lattice);
  mdp_matrix c2(box[0],1);
  for(int t=0; t<c2.size(); t++) c2(t)=0;
  coefficients quark;
  quark["kappa"]=1.1;
  quark["c_{sw}"]=0.3;
  U.load(argv[1]);  
  compute_em_field(U);
  generate(S,U,quark);
  forallsites(x) 
    for(int mu=1; mu<3; mu++) {
      mdp_matrix G1=Gamma[mu]*Gamma5;
      for(int a=0; a<4; a++)
	for(int b=0; b<4; b++)
	  for(int c=0; c<4; c++)
	    for(int d=0; d<4; d++)
	      c2(x(TIME))-=real(trace(G1(a,b)*S(x,b,c)*G1(c,d)*
				      hermitian(S(x,d,a))));
    }
  mdp.add(c2);
  for(int t=0; t<c2.size(); t++) 
    mdp << t << "\t" << c2(t) << endl;
  mdp.close_wormholes();
  return 0;
}
