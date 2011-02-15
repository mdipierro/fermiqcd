#include "fermiqcd.h"
#include "mdp_all.h"
#include "dump.h"
#include "fermiqcd_topological_charge.hpp"

class PunchedWilsonGaugeAction : public WilsonGaugeAction {
 public:
  static gauge_stats heatbath(gauge_field &U, 
			      coefficients &coeff, 
			      vector<mdp_site> &sites,
			      int n_iter=1) {
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

    for(iter=0; iter<n_iter; iter++) 
      for(parity=0; parity<2; parity++) 
	for(mu=0; mu<U.ndim; mu++) {
	  forallsitesofparity(x,parity) {
	    if(mu==0) for(int q=0; q<sites.size(); q++) {
	      if(x(1)==sites[q](1) &&
		 x(2)==sites[q](2) &&
		 x(3)==sites[q](3)
		 ) continue;
	    }
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


void punched_ape_smearing(gauge_field &U, 
			  vector<mdp_site> &sites,
			  mdp_real alpha=0.7, 
			  int iterations=20, 
			  int cooling_steps=10) {
  gauge_field V(U.lattice(),U.nc);
  mdp_site x(U.lattice());
  for(int iter=0; iter<iterations; iter++) {
    cout << "smearing step " << iter << "/" << iterations << endl;
    V=U;
    for(int mu=0; mu<4; mu++) {
      forallsites(x) {
	if(mu==0) for(int q=0; q<sites.size(); q++) {
	  if(x(1)==sites[q](1) &&
	     x(2)==sites[q](2) &&
	     x(3)==sites[q](3)
	     ) continue;
	}	
	U(x,mu)=(1.0-alpha)*V(x,mu);
	for(int nu=0; nu<U.ndim; nu++)
	  if(nu!=mu)
	    U(x,mu)+=(1.0-alpha)/6*
	      (V(x,nu)*V(x+nu,mu)*hermitian(V(x+mu,nu))+
	       hermitian(V(x-nu,nu))*V(x-nu,mu)*V((x-nu)+mu,nu));
	U(x,mu)=project_SU(U(x,mu),cooling_steps);
      }
    }
    U.update();
  }
}

int main(int argc, char** argv) {
  mdp.open_wormholes(argc,argv);
  define_base_matrices("FERMILAB");
  int N,nc=2;
  int L[]={10,16,10,10};
  int idx[6][3]={{0,0,1},{1,0,2},{2,0,3},{3,1,2},{4,1,3},{5,2,3}};
  mdp_lattice lattice(4,L);
  mdp_lattice cube(3,L+1);

  gauge_field U(lattice,nc);
  gauge_field V(lattice,nc);
  gauge_field W(lattice,nc);
  mdp_complex_field PA(cube,6);  
  mdp_complex_field PB(cube,6);  
  mdp_complex P, PC=0;
  mdp_field<float> Q3(cube);
  mdp_site x(lattice), y(lattice);
  mdp_site x3(cube), x3b(cube);
  mdp_complex d; 
  char filename[100];
  coefficients gauge; gauge["beta"]=2.2;
  vector<mdp_site> sites;

  int path[20][2]={{+1,0},{+1,0},{+1,0},{+1,0},{+1,0},
		   {+1,1},{+1,1},{+1,1},{+1,1},{+1,1},
		   {-1,0},{-1,0},{-1,0},{-1,0},{-1,0},
		   {-1,1},{-1,1},{-1,1},{-1,1},{-1,1}};
  
  forallsites(x3)
    for(int i=0; i<6; i++)
      PA(x3,i)=PB(x3,i)=0;
  PC=0;
  
  /*
  set_hot(V);
  for(int k=0; k<4000; k++) {
    cout << k << endl;
    WilsonGaugeAction::heatbath(V,gauge,10);
  }
  */
  

  for(int conf=0; conf<1000; conf++) {
    /* WilsonGaugeAction::heatbath(V,gauge,50);    
    sprintf(filename,"gauge_from_hot_10x16x10x10.%.3i.fixed.mdp",conf);
    V.save(filename);

    U=V;
    ape_smearing(U,0.7,20,10);
    */
    sprintf(filename,"gauge_from_hot_10x16x10x10.%.3i.fixed.mdp",conf);
    U.load(filename);    

    // GaugeFixing::fix(U,GaugeFixing::Landau,20);
 
    // for(int shift0=0; shift0<L[0]; shift0++) {
    // for(int shift1=0; shift1<L[1]; shift1++) {
    // for(int shift2=0; shift2<L[2]; shift2++) {
    for(int shift3=0; shift3<L[3]; shift3++) {
      
      x.set(3,5,L[2]/2,L[3]/2);

	d=trace(build_path(U,x,20,path));
	PC+=d;
	forallsites(x3) {
	  x.set(x3(2),x3(0),x3(1),L[3]/2);
	  for(int h=0; h<6; h++) {
	    P=trace(plaquette(U,x,idx[h][1],idx[h][2]));    
	    PA(x3,idx[h][0])+=d*P;
	    PB(x3,idx[h][0])+=P;
	  }
	}

	
	forallsites(x) for(int mu=0; mu<4; mu++) W(x,mu)=U(x-3,mu);
	W.update(); U=W;
      }
    /*
	forallsites(x) for(int mu=0; mu<4; mu++) W(x,mu)=U(x-2,mu);
	W.update(); U=W;
      }
	  forallsites(x) for(int mu=0; mu<4; mu++) W(x,mu)=U(x-1,mu);
	W.update(); U=W;
      }
	
      forallsites(x) for(int mu=0; mu<4; mu++) W(x,mu)=U(x-0,mu); 
      W.update(); U=W;
      }
    */
	
    N=(conf+1)*L[0];
    float m=0;
    float Q;
    forallsites(x3) Q3(x3)=0;
    forallsites(x3) {
      Q=0;
      for(int q=0; q<3; q++)
	Q+=(real(PA(x3,q)/PC-PB(x3,q)/N));
      for(int q=3; q<6; q++)
	Q+=(real(PA(x3,q)/PC-PB(x3,q)/N));      
      Q3(x3)+=Q;      
      //x3b.set(L[1]-x3(0),x3(1),x3(2));      Q3(x3b)+=Q;
      //x3b.set(x3(0),L[2]-x3(1),x3(2));      Q3(x3b)+=Q;
      //x3b.set(x3(0),x3(1),L[3]-x3(2));      Q3(x3b)+=Q;
      //x3b.set(L[1]-x3(0),L[2]-x3(1),x3(2));      Q3(x3b)+=Q;
      //x3b.set(L[1]-x3(0),x3(1),L[3]-x3(2));      Q3(x3b)+=Q;
      //x3b.set(x3(0),L[2]-x3(1),L[3]-x3(2));      Q3(x3b)+=Q;
      //x3b.set(L[1]-x3(0),L[2]-x3(1),L[3]-x3(2));      Q3(x3b)+=Q;      
    }
    // x3.set(0,0,0); Q3(x3)=0;
    forallsites(x3) if(Q3(x3)<m) m=Q3(x3);
    cout << "m=" << m << endl;
    forallsites(x3) Q3(x3)-=m;    
    cout << "saving vtk file\n";
    sprintf(filename,"energy_density_E2B2_XYT_Wilson5x5.%.3i.vtk", conf);
    dump(Q3,0,filename);
  }
  
  mdp.close_wormholes();
  return 0;
}
