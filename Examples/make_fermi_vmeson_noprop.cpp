#include "fermiqcd.h"

int main(int argc, char** argv) {
  mdp.open_wormholes(argc,argv);    
  define_base_matrices("FERMILAB");
  int nc=3;
  int box[]={8,4,4,4};
  mdp_lattice lattice(4,box);
  gauge_field U(lattice,nc);
  fermi_field source(lattice,nc);
  fermi_field sink(lattice,nc);
  mdp_site x(lattice);
  mdp_matrix c2(box[0],1);
  for(int t=0; t<c2.size(); t++) c2(t)=0;
  coefficients quark;
  quark["kappa"]=1.1;
  quark["c_{sw}"]=0.3;
  U.load(argv[1]);  
  compute_em_field(U);
  for(int mu=0; mu<4; mu++)
    for(int alpha=0; alpha<4; alpha++)
      for(int i=0; i<nc; i++) {
	source=0;
	forallsites(x) 
	  for(int beta=0; beta<4; beta++)
	    source(x,alpha,i)=Gamma[mu](beta,alpha);
	mul_invQ(sink,source,U,quark);
	for(int beta=0; beta<4; beta++)
	  for(int j=0; j<nc; j++) {
	    forallsites(x) 
	      c2(x(TIME))+=real(pow(sink(x,beta,j),2));
	  }
      }
  mdp.add(c2);
  for(int t=0; t<c2.size(); t++) 
    mdp << t << "\t" << c2(t) << endl;
  mdp.close_wormholes();
  return 0;
}
