#include "fermiqcd.h"
#include "mdp_all.h"
#include "dump.h"
#include "fermiqcd_topological_charge.hpp"

class ModifiedWilsonGaugeAction : public WilsonGaugeAction {
 public:
  static gauge_stats heatbath(gauge_field &U, 
			      coefficients &coeff, 
			      int n_iter=1, int mu_set=0) {
    begin_function("WilsonGaugeAction__heatbath");
    if(U.nc==1) error("fermiqcd_gauge_algorithms/heatbath(): U(1)? (use metropolis)");
    gauge_stats stats;
    mdp_real beta, zeta;
    if(coeff.has_key("beta")) beta=coeff["beta"]; 
    else error("beta undeclared");
    if(coeff.has_key("zeta")) zeta=coeff["zeta"]; 
    else zeta=1;

    int i,j,k,iter,mu,parity;
    mdp_matrix M;
    mdp_complex a[4], tmpUik;
    site x(U.lattice());
    double time=mpi.time();

    mdp << coeff;

    forallsitesofparity(x,n_iter) {
      
      mu=mu_set;
      {
	
	for(i=0; i<U.nc-1; i++)
	  for(j=i+1; j<U.nc; j++) {
	    if(zeta==1)    M=U(x,mu)*staple_H(U,x,mu);
	    else if(mu==0) M=zeta*U(x,0)*staple_H(U,x,0);
	    else           M=((mdp_real) 1.0/zeta)*U(x,mu)*staple_H(U,x,mu);
	    a[0]=M(i,i); 
	    a[1]=M(i,j);
	    a[2]=M(j,i);
	    a[3]=M(j,j);
	    heatbath_SU2(U.lattice().random(x),beta/U.nc,a);
	    for(k=0; k<U.nc; k++) {
	      tmpUik=a[0]*U(x,mu,i,k)+a[1]*U(x,mu,j,k);
	      U(x,mu,j,k)=a[2]*U(x,mu,i,k)+a[3]*U(x,mu,j,k);
	      U(x,mu,i,k)=tmpUik;
	    }
	  }
      }
      // The next command does all the communications!
      U.update(parity, mu, U.nc*U.nc);
    }
    mdp << "\t<stats>\n\t\t<time>" << mpi.time()-time << "</time>\n\t</stats>\n";
    end_function("WilsonGaugeAction__heatbath");
    return stats;
  }
};

void save_top_charge(gauge_field &U, int code, int c1, int c2, int tmin, int tmax) {
  mdp_lattice &lattice=U.lattice();
  int L[3]={lattice.size(1),lattice.size(2),lattice.size(3)};
  mdp_lattice cube(3,L);
  mdp_field<float> Q3(cube);
  mdp_site x4(lattice);
  mdp_site x3(cube);
  mdp_field<float> Q4(lattice);
  char filename[128];

  for(int i1=0; i1<c1; i1++) {
    ape_smearing(U,0.7,c2,10);
    cout << "topological charge...\n";
    topological_charge(Q4,U);
    Q4.save("sample_topological_charge.mdp");
    cout << "3d projection...\n";
    for(int t=tmin; t<tmax; t++) {
      forallsites(x3) {
	x4.set(t,x3(0),x3(1),x3(2));
	Q3(x3)=Q4(x4);
      }
      cout << "dumping...\n";
      // sprintf(filename,"topological.%.3i.%.3i.%.3i.vtk",code,i1,t);
      sprintf(filename,"topological.%.3i.vtk",code);
      dump(Q3,0,filename);
    }  
  }
}

int main(int argc, char** argv) {
  mdp.open_wormholes(argc,argv);
  define_base_matrices("FERMILAB");
  int nc=3;
  // int L[]={10,10,10,10}; 
  int L[]={8,12,12,12};
  mdp_lattice lattice(4,L);
  mdp_site x(lattice);
  gauge_field U(lattice,nc);
  gauge_field V(lattice,nc);
  char filename[128];  
  mdp_matrix A0=mdp_random.SU(3);
  mdp_matrix A1=mdp_random.SU(3);
  mdp_matrix A2=mdp_random.SU(3);
  mdp_matrix A3=mdp_random.SU(3);

  coefficients gauge; gauge["beta"]=5.6;

  set_cold(U);
  forallsites(x) if(x(1)>=L[1]/4 && x(1)<3*L[1]/4) {
    U(x,0)=A0;
    U(x,1)=A1;
    U(x,2)=A2;
    U(x,3)=A3;
  }
  save_top_charge(U, 0, 1, 20, 0, 1);

  for(int k=1; k<1000; k++) { 
    cout << k << endl;
    ModifiedWilsonGaugeAction::heatbath(U,gauge,(k%8)/2,(k%8)%4);
    sprintf(filename,"gauge_8x12x12x12.%.5i.mdp",k);
    U.save(filename);
    V=U;
    save_top_charge(V, k, 1, 20, 0, 1);
  }
  /*
  for(int k=10000; k<10100; k++)
    save_top_charge(U, k, 1, 1, 0, 1);
  save_top_charge(U, 10100, 1, 0, 0, L[0]);
  */
  mdp.close_wormholes();
  return 0;
}
