/////////////////////////////////////////////////////////////////
/// @file fermiqcd_fermi_actions.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Basic actions for Wilson Fermions
///
/// Distributed under GPL2 License
/// 
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////

/// @brief Wilson/Clover action (SLOW: DO NOT USE IN PRODUCTION)
///
/// Example:
/// @verbatim
/// gauge_field U(lattice,nc);
/// fermi_field psi(lattice,nc);
/// fermi_field chi(lattice,nc);
/// coefficients coeff;
/// coeff["kappa_s"]=0.11;
/// coeff["kappa_t"]=0.11;
/// coeff["r_s"]=1.0;
/// coeff["r_t"]=1.0;
/// coeff["c_{sw}"]=1.0;
/// coeff["c_E"]=1.0;
/// coeff["c_B"]=1.0;
/// default_fermi_action=FermiCloverActionSlow::mul_Q;
/// if(coeff["c_{sw}"]!=0) compute_em_field(U);
/// mul_Q(chi,psi,U,coeff);
/// @endverbatim
/// Note that mul_Q(chi,psi,U,coeff) reads \f$ \chi=(/\!\!\!D[U]+m)\psi \f$ 
class FermiCloverActionSlow {
 public:
  static void mul_Q(fermi_field &psi_out, 
		    fermi_field &psi_in, 
		    gauge_field &U, 
		    coefficients &coeff,
		    int parity=EVENODD) {

    if(psi_in.nspin!=4) error("fermiqcd_fermi_algorithms/mul_Q_ONE: nspin!=4");
    if(psi_in.nc!=U.nc) error("fermiqcd_fermi_algorithms/mul_Q_ONE: gauge and spinor have different nc");
    if(parity!=EVENODD) error("parity must be EVENODD here");
    
    register int   ndim=psi_in.lattice().ndim;
    register int   nspin=psi_in.nspin;
    register int   nc=psi_in.nc;
    register mdp_real kappa_t=0;
    register mdp_real kappa_s=0;
    register mdp_real r_t;
    register mdp_real r_s;
    register mdp_real cSW;
    register mdp_real c_E;
    register mdp_real c_B;
    register mdp_real rsign;

    if(coeff.has_key("kappa")) kappa_s=kappa_t=coeff["kappa"];
    if(coeff.has_key("kappa_t")) kappa_t=coeff["kappa_t"];
    if(coeff.has_key("kappa_s")) kappa_s=coeff["kappa_s"];
    if(kappa_t==0) error("kappa_t=0 or undeclared");
    if(kappa_s==0) error("kappa_s=0 or undeclared");
    if(coeff.has_key("r_t")) r_t=coeff["r_t"];       else r_t=1;
    if(coeff.has_key("r_s")) r_s=coeff["r_s"];       else r_s=1;
    if(coeff.has_key("c_{sw}")) cSW=coeff["c_{sw}"]; else cSW=0;
    if(coeff.has_key("c_E")) c_E=coeff["c_E"];       else c_E=1;
    if(coeff.has_key("c_B")) c_B=coeff["c_B"];       else c_B=1;
    if(coeff.has_key("sign")) rsign=coeff["sign"];   else rsign=+1;

    site           x(psi_in.lattice());
    register int   a,mu,nu;
    
    mdp_matrix psi_up(nspin,nc);
    mdp_matrix psi_dw(nspin,nc);
    mdp_matrix psi_lo(nspin,nc);
    
    // =================================================================
    //  The signs are now in accord with the definitions of
    //  MacKenzie et al. hep-lat/9604004
    //    and are in agreement with  
    //  Heatlie et al., Nucl.Phys. B352 (1991) 266
    //    care here:
    //    my definition of sigma is sigma_munu=i/2[gamma_mu,gamma_nu]
    //    an extra phase when compared with Heatlie         
    // ================================================================= 
    
    psi_out=psi_in; 
    forallsites(x) {
      for(mu=0; mu<ndim; mu++) {
	for(a=0; a<nspin; a++) {
	  psi_up(a)=U(x,mu)*psi_in(x+mu,a);
	  psi_dw(a)=hermitian(U(x-mu,mu))*psi_in(x-mu,a);
	}
	if(mu==0) psi_out(x)-=(kappa_t)*((r_t-Gamma[mu]*rsign)*psi_up+
					 (r_t+Gamma[mu]*rsign)*psi_dw);
	else      psi_out(x)-=(kappa_s)*((r_s-Gamma[mu]*rsign)*psi_up+
					 (r_s+Gamma[mu]*rsign)*psi_dw);
      }
      
      if(cSW!=0)
	for(mu=0; mu<ndim-1; mu++)
	  for(nu=mu+1; nu<ndim; nu++) {
	    for(a=0; a<nspin; a++) 
	      psi_lo(a)=U.em(x,mu,nu)*psi_in(x,a);
	    if(mu==0) {
	      psi_out(x)+=(kappa_s*cSW*c_E*I)*(Sigma[mu][nu]*psi_lo);
	    } else {
	      psi_out(x)+=(kappa_s*cSW*c_B*I)*(Sigma[mu][nu]*psi_lo);
	    }
	  }
    }
  }
};

/// @brief Wilson/Clover action
///
/// Example:
/// @verbatim
/// gauge_field U(lattice,nc);
/// fermi_field psi(lattice,nc);
/// fermi_field chi(lattice,nc);
/// coefficients coeff;
/// coeff["kappa_s"]=0.11;
/// coeff["kappa_t"]=0.11;
/// coeff["r_s"]=1.0;
/// coeff["r_t"]=1.0;
/// coeff["c_{sw}"]=1.0;
/// coeff["c_E"]=1.0;
/// coeff["c_B"]=1.0;
/// default_fermi_action=FermiCloverActionFast::mul_Q;
/// if(coeff["c_{sw}"]!=0) compute_em_field(U);
/// mul_Q(chi,psi,U,coeff);
/// @endverbatim
class FermiCloverActionFast {
 public:
  static void mul_Q(fermi_field &psi_out, 
		    fermi_field &psi_in, 
		    gauge_field &U, 
		    coefficients &coeff,
		    int parity=EVENODD) {
    
    
    if(psi_in.nspin!=4) error("fermiqcd_fermi_algorithms/mul_Q_TWO: nspin!=4");
    if(psi_in.nc!=U.nc) error("fermiqcd_fermi_algorithms/mul_Q_TWO: gauge and spinor have different nc");
    if(parity!=EVENODD) error("parity must be EVENODD here");

    register int   ndim=psi_in.lattice().ndim;
    register int   nspin=psi_in.nspin;
    register int   nc=psi_in.nc;
    register mdp_real kappa_t=0;
    register mdp_real kappa_s=0;
    register mdp_real r_t;
    register mdp_real r_s;
    register mdp_real cSW;
    register mdp_real c_E;
    register mdp_real c_B;
    register mdp_real rsign;

    if(coeff.has_key("kappa")) kappa_s=kappa_t=coeff["kappa"];
    if(coeff.has_key("kappa_t")) kappa_t=coeff["kappa_t"];
    if(coeff.has_key("kappa_s")) kappa_s=coeff["kappa_s"];
    if(kappa_t==0) error("kappa_t=0 or undeclared");
    if(kappa_s==0) error("kappa_s=0 or undeclared");
    if(coeff.has_key("r_t")) r_t=coeff["r_t"];       else r_t=1;
    if(coeff.has_key("r_s")) r_s=coeff["r_s"];       else r_s=1;
    if(coeff.has_key("c_{sw}")) cSW=coeff["c_{sw}"]; else cSW=0;
    if(coeff.has_key("c_E")) c_E=coeff["c_E"];       else c_E=1;
    if(coeff.has_key("c_B")) c_B=coeff["c_B"];       else c_B=1;
    if(coeff.has_key("sign")) rsign=coeff["sign"];   else rsign=+1;
    
    site           x(psi_in.lattice());
    register int   i,j,a,mu,nu,inc,anc,bnc;
    mdp_real         coeff_kappa=0;
    mdp_real         coeff_sum=0;
    mdp_complex        coeff_dif=0;
    mdp_complex        coeff_clover=0;
    mdp_complex        psi_up, psi_dw, psi_loc;
    mdp_complex        *psi_tmp=new mdp_complex[nspin*nc]; 
    mdp_complex        *Fem,   *Fpsi_in;
    mdp_complex        *FU_up, *Fpsi_in_up;
    mdp_complex        *FU_dw, *Fpsi_in_dw;
    mdp_complex        *FHU_dw=new mdp_complex[nc*nc];
    
    // =================================================================
    //  The signs are now in accord with the definitions of
    //  MacKenzie et al. hep-lat/9604004
    //    and are in agreement with  
    //  Heatlie et al., Nucl.Phys. B352 (1991) 266
    //    care here:
    //    my definition of sigma is 
    //      Sigma[mu][nu]=iI2*[Gamma[mu],Gamma[nu]]
    //    an extra phase when compared with Heatlie         
    // ================================================================= 
    
    psi_out=psi_in;
    
    forallsites(x) {
      for(a=0; a<nspin; a++) 
	for(i=0; i<nc; i++) {
	  psi_tmp[anc=a*nc+i]=0;
	}
      for(mu=0; mu<ndim; mu++) { 
	if(mu==0) {
	  coeff_kappa=kappa_t;
	  coeff_sum=coeff_kappa*r_t;
	} else {
	  coeff_kappa=kappa_s;
	  coeff_sum=coeff_kappa*r_s;
	}
	FU_up=&U(x,mu,0,0);
	FU_dw=&U(x-mu,mu,0,0);
	for(i=0; i<nc; i++) 
	  for(j=0; j<nc; j++) 
	    FHU_dw[i*nc+j]=conj(FU_dw[j*nc+i]);
	for(a=0; a<nspin; a++) {
	  anc=a*nc;
	  Fpsi_in_up=&psi_in(x+mu,a,0);
	  Fpsi_in_dw=&psi_in(x-mu,a,0);
	  bnc=Gamma_idx[mu][a]*nc;
	  coeff_dif=coeff_kappa*Gamma_val[mu][a]*rsign;
	  for(i=0; i<nc; i++) {
	    inc=i*nc;
	    psi_up= FU_up[inc]*Fpsi_in_up[0];
	    psi_dw=FHU_dw[inc]*Fpsi_in_dw[0];
	    for(j=1; j<nc; j++) {
	      psi_up+= FU_up[inc+j]*Fpsi_in_up[j];
	      psi_dw+=FHU_dw[inc+j]*Fpsi_in_dw[j];
	    }
	    psi_tmp[anc+i]-=coeff_sum*(psi_dw+psi_up);
	    psi_tmp[bnc+i]-=coeff_dif*(psi_dw-psi_up);
	  }
	}
      }
      if(cSW!=0) {
	for(mu=0; mu<ndim-1; mu++)
	  for(nu=mu+1; nu<ndim; nu++)
	    for(a=0; a<nspin; a++) {
	      Fem=&(U.em(x,mu,nu,0,0));
	      bnc=Sigma_idx[mu][nu][a]*nc;
	      Fpsi_in=&psi_in(x,a,0);
	      if(mu==0) coeff_clover=(kappa_s*cSW*c_E*I)*Sigma_val[mu][nu][a]; 
	      else      coeff_clover=(kappa_s*cSW*c_B*I)*Sigma_val[mu][nu][a];
	      for(i=0; i<nc; i++) {
		inc=i*nc;
		psi_loc=Fem[inc]*Fpsi_in[0];
		for(j=1; j<nc; j++)
		  psi_loc+=Fem[inc+j]*Fpsi_in[j];
		psi_tmp[bnc+i]+=coeff_clover*psi_loc;
	      }
	    }
      }
      for(a=0;a<nspin; a++) 
	for(i=0; i<nc; i++) 
	  psi_out(x,a,i)+=psi_tmp[anc=a*nc+i];
    }
    delete[] psi_tmp;
    delete[] FHU_dw;
  }
};

