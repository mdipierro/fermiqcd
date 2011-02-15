
void rotate_field(fermi_field& psi, gauge_field& U, coefficients &coeff) {
  fermi_field psi_rot(psi.lattice(),psi.nspin,psi.nc);
  mdp_real d1=coeff["rotate_field:d1"];
  psi_rot=psi;
  mdp_site x(psi.lattice());
  forallsites(x)
    for(int mu=1;mu<U.ndim;mu++)
      for(int a=0;a<psi.nspin;a++)
        for(int b=0;b<psi.nspin;b++) {
          psi(x,a)+=d1*0.5*Gamma[mu](a,b)*(U(x,mu)*psi_rot(x+mu,b)-(U(x,-1,mu))*psi_rot(x-mu,b));
        }
  psi.update(); 
}

