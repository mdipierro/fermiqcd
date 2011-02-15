/*
Under development - DOES NOT WORK!!!!

class HypSmearing {
 public:
  vector<float> alpha;
  bool in(int x, vector<int> set) {
    for(int i=0; i<set.size(); i++) 
      if(x==set[i]) return true;
    return false;
  }
 public: static void smear_aux(gauge_field &U, 
			       vector<int> set,
			       int cooling_steps=10) {   
   mdp_site x(U.lattice());
   mdp_site y(U.lattice());
   mdp_matrix A(U.nc,U.nc)
   mdp_matrix staples(U.nc,U.nc)
   for(int m=0; m<U.ndim-1; m++) {
     forallsites(x)
       for(int mu=0; mu<U.ndim; mu++) {
	 for(int i=0; i<U.nc; i++)
	   for(int j=0; j<U.nc; j++) 
	     staples(i,j)=0;
	 for(int nu=0; nu<U.nc; nu++)
	   y=x+nu;
	   for(int i=0; i<U.nc; i++)
	     for(int j=0; j<U.nc; j++) {
	       A(i,j)=0;
	       for(int k=0; k<U.nc; k+) A(i,j)=U(x,nu,i,k)*U(y,mu,k,j);
	     }
	   for(int i=0; i<U.nc; i++)
	     for(int j=0; j<U.nc; j++) {	     
	       for(int k=0; k<U.nc; k+) staples(i,j)+=A(i,k)*conj(U(x+mu,nu,j,k));
	     }
	   y=x-nu;
	   for(int i=0; i<U.nc; i++)
	     for(int j=0; j<U.nc; j++) {
	       A(i,j)=0;
	       for(int k=0; k<U.nc; k+) A(i,j)=U(x-nu,nu,k,i)*U(y,mu,k,j);
	     }
	   for(int i=0; i<U.nc; i++)
	     for(int j=0; j<U.nc; j++) {	     
	       for(int k=0; k<U.nc; k+) staples(i,j)+=A(i,k)*U(y+mu,nu,k,j);
	     }	   
	 }       
	 for(int i=0; i<U.nc; i++)
	   for(int j=0; j<U.nc; j++) 
	     U(x,mu,i,j)=(1.0-alpha[m])*U(x,mu,i,j)+alpha[m]*staples(i,j);
	 U(x,mu)=ProjectSUN(U(x,mu));
       }
   }
};
*/

// from Bonnet et al. Phys Rev D 62, 094509

class ApeSmearing {
  public: static void smear(gauge_field &U, 
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
};


void compute_em_notrace_field(gauge_field &U) {
  compute_em_field(U);
  mdp_site x(U.lattice());
  forallsitesandcopies(x)
    for(int mu=0; mu<U.ndim-1; mu++)
      for(int nu=mu+1; nu<U.ndim; nu++) 
	U.em(x,mu,nu)-=8.0/3.0*I*trace(U.em(x,mu,nu));
}

void topological_charge(mdp_field<float> &Q, gauge_field &U) {
  compute_em_notrace_field(U);
  mdp_site x(U.lattice());
  forallsitesandcopies(x) {
    Q(x)=0;
    for(int i=0; i<U.nc; i++)
      for(int j=0; j<U.nc; j++)
	Q(x)+=real(U.em(x,0,1,i,j)*U.em(x,2,3,j,i)-
		   U.em(x,0,2,i,j)*U.em(x,1,3,j,i)+
		   U.em(x,0,3,i,j)*U.em(x,1,2,j,i));
  }
  Q.update();
}

float topological_charge_vtk(gauge_field &U, string filename, int t=-1) {
  mdp_field<float> Q(U.lattice()), P(U.lattice());
  mdp_site x(U.lattice());
  topological_charge(Q,U);
  Q.save_vtk(filename,t);
  cumulate_field(Q,filename).save_vtk(filename.replace(filename.rfind("."),1,".sum."),t);
  for(int t=U.lattice().size(3)-1; t>=0; t--) {    
    forallsites(x) if(x(0)==t) {
      if(t==U.lattice().size(3)-1) P(x)=Q(x);
      else P(x)=P(x+0)+Q(x);      
    }
    P.update();
  }
  P.save_vtk(filename.replace(filename.rfind(".sum."),5,".flat."),0);
  double total=0.0;
  forallsites(x) total+=Q(x);
  mdp.add(total);
  return total;
}
