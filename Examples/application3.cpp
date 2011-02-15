// Program: application3.cpp
#include "mdp.h"
 
// /////////////////////////////////////////////////////
class scalar_field: public mdp_field<float> {
public:
  int ndim,nc;
  scalar_field(mdp_lattice &a) {
    allocate_field(a, 1);
    ndim=a.ndim;
  };
  float &operator() (mdp_site x) {
    return *(address(x));
  };
};

// /////////////////////////////////////////////////////
class ising_field: public mdp_field<int> {
public:
  int ndim,nc;
  float beta, kappa, magnetic_field;
  ising_field(mdp_lattice &a) {
    allocate_field(a, 1);
    ndim=a.ndim;
  };
  int &operator() (mdp_site x) {
    return *(address(x));
  };
  friend void set_cold(ising_field &S) {
    mdp_site x(S.lattice());
    forallsites(x) S(x)=1;
  };

  friend void montecarlo_multihit(ising_field &S, scalar_field &H, 
				  int n_iter=1, int n_hits=1) {
    int iter, parity, hit, new_spin, mu;
    float delta_action;
    mdp_site x(S.lattice());
    for(iter=0; iter<n_iter; iter++) {
      for(parity=0; parity<=1; parity++)  
        forallsitesofparity(x,parity) {
	for(hit=0; hit<n_hits; hit++) {
	  delta_action=S.kappa*(H(x)+S.magnetic_field);
	  for(mu=0; mu<S.ndim; mu++)
	    delta_action-=S(x+mu)+S(x-mu);
	  new_spin=(S.lattice().random(x).plain()>0.5)?1:-1;
	  delta_action*=S.beta*(new_spin-S(x));
	  if(delta_action<S.lattice().random(x).plain()) 
	    S(x)=new_spin;
	};
      };
      S.update(parity);
    };
  };
  friend float average_spin(ising_field &S) {
    float res=0;
    mdp_site x(S.lattice());
    forallsites(x) res+=S(x);
    mdp.add(res);
    return res/(S.lattice().nvol_gl);
  };
};

int main(int argc, char **argv) {
  mdp.open_wormholes(argc, argv);
  
  int conf,Nconfig=10;
  int mybox[]={128, 128};
  mdp_lattice mylattice(2, mybox);
  ising_field     S(mylattice); /* ths spin field +1 or -1     */
  scalar_field    H(mylattice); /* the external magnetic field */
  mdp_site        x(mylattice);
  mdp_jackboot    jb(Nconfig,1); 
  jb.plain(0);  

  S.beta=0.5;         /* inverse square temperature             */
  S.kappa=0.1;        /* coupling with the total extarnal field */
  S.magnetic_field=1; /* an extra external magnetic field       */

  for(S.beta=0.01; S.beta<0.15; S.beta+=0.001) {
    forallsites(x) {
      S(x)=1;
      H(x)=0;
    };
    S.update();
    H.update();
    
    for(conf=0; conf<Nconfig; conf++) {
      montecarlo_multihit(S,H);                  
      jb(conf,0)=average_spin(S);
    };
    mdp << "beta = " << S.beta << ", "
	<< "average plaquette = " << jb.mean() 
	<< "(" << jb.j_err() << ")\n";
  };
  mdp.close_wormholes();
  return 0;
};
