/////////////////////////////////////////////////////////////////
/// @file fermiqcd_dwfermi_actions.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Actions for Domain Wall fermions
///
/// Distributed under GPL2 License
/// 
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////

/// @brief domain wall action (SORRY THIS IS SLOW)
///
/// Notation from ref. hep-lat/0007038
/// Example:
/// @verbatim
/// gauge_field U(lattice,nc);
/// dwfermi_field psi(lattice,nc);
/// dwfermi_field chi(lattice,nc);
/// coefficients coeff;
/// coeff["m_f"]=0.11; // fermion mass
/// coeff["m_5"]=0.11; // mass in 5th dimension
/// default_dwfermi_action=DWFermiActionSlow::mul_Q;
/// mul_Q(chi,psi,U,coeff);
/// @endverbatim
/// Note that mul_Q(chi,psi,U,coeff) reads \f$ \chi=(/\!\!\!D[U]+m)\psi \f$ 

class DWFermiActionSlow {
 public:
  static void mul_Q(dwfermi_field &psi_out, 
		    dwfermi_field &psi_in, 
		    gauge_field &U, 
		    coefficients &coeff) {
    
    if(psi_in.nspin!=4) error("fermiqcd_dwfermi_algorithms/dwfermi_mul_Q_ONE: nspin!=4");
    if(psi_in.nc!=U.nc) error("fermiqcd_dwfermi_algorithms/dwfermi_mul_Q_ONE: gauge and spinor have different nc");
    
    register int    ndim=psi_in.lattice().ndim;
    register int    nspin=psi_in.nspin;
    register int    nc=psi_in.nc;
    register int    L5=psi_in.L5;
    register mdp_real m_5,m_f,sign;
    if(coeff.has_key("m_5")) m_5=coeff["m_5"];
    else error("coefficients m_5 undeclared");
    if(coeff.has_key("m_f")) m_f=coeff["m_f"];
    else error("coefficients m_f undeclared");
    if(coeff.has_key("sign")) sign=coeff["sign"];
    else sign=1;

    // check the sign and the 6.0 here
    register mdp_real kappa5=0.5/(m_5-6.0);
    register mdp_real kappaf=-m_f*kappa5;
    
    site           x(psi_in.lattice());
    register int   l,a,mu;
    
    mdp_matrix psi_up(nspin,nc);
    mdp_matrix psi_dw(nspin,nc);
    mdp_matrix psi_lo(nspin,nc);
    
    // =================================================================
    // ref. hep-lat/0007038
    // ================================================================= 
    
    mdp_matrix tmp;
    
    psi_out=psi_in; 
    forallsites(x) {
      for(l=0; l<L5; l++) {
	for(mu=0; mu<ndim; mu++) {
	  for(a=0; a<nspin; a++) {
	    psi_up(a)=U(x,mu)*psi_in(x+mu,l,a);
	    psi_dw(a)=hermitian(U(x-mu,mu))*psi_in(x-mu,l,a);
	  }
	  psi_out(x,l)+=kappa5*((1-sign*Gamma[mu])*psi_up+
				(1+sign*Gamma[mu])*psi_dw);
	}
	
	if(l<L5-1)
	  psi_out(x,l)+=kappa5*(1-Gamma5)*psi_in(x,l+1);
	else
	  psi_out(x,L5-1)+=kappaf*(1-Gamma5)*psi_in(x,0);
	if(l>0)
	  psi_out(x,l)+=kappa5*(1+Gamma5)*psi_in(x,l-1);
	else
	  psi_out(x,0)+=kappaf*(1+Gamma5)*psi_in(x,L5-1);
	
      }
    }
  }
};


/// @brief domain wall action fast
///
/// Notation from ref. hep-lat/0007038
/// Example:
/// @verbatim
/// gauge_field U(lattice,nc);
/// dwfermi_field psi(lattice,nc);
/// dwfermi_field chi(lattice,nc);
/// coefficients coeff;
/// coeff["m_f"]=0.11; // fermion mass
/// coeff["m_5"]=0.11; // mass in 5th dimension
/// default_dwfermi_action=DWFermiActionFast::mul_Q;
/// mul_Q(chi,psi,U,coeff);
/// @endverbatim
/// Note that mul_Q(chi,psi,U,coeff) reads \f$ \chi=(/\!\!\!D[U]+m)\psi \f$ 

class DWFermiActionFast {
 public:
  static void mul_Q(dwfermi_field &psi_out, 
		    dwfermi_field &psi_in, 
		    gauge_field &U, 
		    coefficients &coeff) {
    
    if(psi_in.nspin!=4) error("fermiqcd_dwfermi_algorithms/dwfermi_mul_Q_ONE: nspin!=4");
    if(psi_in.nc!=U.nc) error("fermiqcd_dwfermi_algorithms/dwfermi_mul_Q_ONE: gauge and spinor have different nc");
    
    register int    ndim=psi_in.lattice().ndim;
    register int    nspin=psi_in.nspin;
    register int    nc=psi_in.nc;
    register int    L5=psi_in.L5;
    register mdp_real m_5,m_f,sign;
    if(coeff.has_key("m_5")) m_5=coeff["m_5"];
    else error("coefficients m_5 undeclared");
    if(coeff.has_key("m_f")) m_f=coeff["m_f"];
    else error("coefficients m_f undeclared");
    if(coeff.has_key("sign")) sign=coeff["sign"];
    else sign=1;

    // check the sign and the 6.0 here
    register mdp_real kappa5=0.5/(m_5-6.0);
    register mdp_real kappaf=-m_f*kappa5;
    
    site           x(psi_in.lattice());
    register int   l,a,mu;
    
    mdp_matrix psi_up(nspin,nc);
    mdp_matrix psi_dw(nspin,nc);
    mdp_matrix psi_lo(nspin,nc);
    
    // =================================================================
    // ref. hep-lat/0007038
    // ================================================================= 
    
    mdp_matrix tmp;
    
    psi_out=psi_in; 
    forallsites(x) {
      for(l=0; l<L5; l++) {
	for(mu=0; mu<ndim; mu++) {
	  for(a=0; a<nspin; a++) 	    
	    for(int i=0; i<nc; i++) {
	      psi_up(a,i)=psi_dw(a,i)=0;
	      for(int j=0; j<nc; j++) {
		psi_up(a,i)+=U(x,mu,i,j)*psi_in(x+mu,l,a,j);
		psi_dw(a,i)+=conj(U(x-mu,mu,j,i))*psi_in(x-mu,l,a,j);
	      }
	      psi_out(x,l,a,i)+=kappa5*(psi_dw(a,i)+psi_up(a,i));
	      psi_out(x,l,Gamma_idx[mu][a],i)+=(kappa5*sign)*Gamma_val[mu][a]*(psi_dw(a,i)-psi_up(a,i));
	      
	    }	 
	}
	
	if(l<L5-1) 
	  for(a=0; a<nspin; a++) 	    
	    for(int i=0; i<nc; i++) {
	      psi_out(x,l,a,i)+=kappa5*psi_in(x,l+1,a,i);
	      psi_out(x,l,Gamma5_idx[a],i)-=kappa5*Gamma5_val[a]*psi_in(x,l+1,a,i);
	    } 
	else
	  for(a=0; a<nspin; a++) 	    
	    for(int i=0; i<nc; i++) {
	      psi_out(x,L5-1,a,i)+=kappaf*psi_in(x,0,a,i);
	      psi_out(x,L5-1,Gamma5_idx[a],i)-=kappaf*Gamma5_val[a]*psi_in(x,0,a,i);
	} 
	if(l>0)
	  for(a=0; a<nspin; a++) 	    
	    for(int i=0; i<nc; i++) {
	      psi_out(x,l,a,i)+=kappa5*psi_in(x,l-1,a,i);
	      psi_out(x,l,Gamma5_idx[a],i)+=kappa5*Gamma5_val[a]*psi_in(x,l-1,a,i);
	    }
	else
	  for(a=0; a<nspin; a++) 	    
	    for(int i=0; i<nc; i++) {
	      psi_out(x,0,a,i)+=kappaf*psi_in(x,L5-1,a,i);
	      psi_out(x,0,Gamma5_idx[a],i)+=kappaf*Gamma5_val[a]*psi_in(x,L5-1,a,i);
	    }
	
      }
    }
  }
};

