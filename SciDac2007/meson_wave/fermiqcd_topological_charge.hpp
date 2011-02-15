// from Bonnet et al. Phys Rev D 62, 094509

void ape_smearing(gauge_field &U, 
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
