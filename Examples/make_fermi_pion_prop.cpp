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
    for(int alpha=0; alpha<4; alpha++)
      for(int beta=0; beta<4; beta++)
	c2(x(TIME))+=real(trace(S(x,alpha,beta)*
				hermitian(S(x,beta,alpha))));
  mdp.add(c2);
  for(int t=0; t<c2.size(); t++) 
    mdp << t << "\t" << c2(t) << endl;
  mdp.close_wormholes();
  return 0;
}
