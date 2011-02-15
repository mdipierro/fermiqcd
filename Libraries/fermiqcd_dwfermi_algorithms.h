/////////////////////////////////////////////////////////////////
/// @file fermiqcd_dwfermi_algorithms.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// More stuff for domain wall fermions
///
/// Distributed under GPL2 License
/// 
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////

/// Projects a domain wall fermion (chi) into a wilson fermion (psi)
void project(fermi_field &psi, dwfermi_field &chi) {
  site x(chi.lattice());
  forallsites(x) {
    psi(x)=Pleft*chi(x,0)+Pright*chi(x,chi.L5-1);
  }
}

/// Projects a will fermion (psi) into a domain wall fermion (chi)
void project(dwfermi_field &chi, fermi_field &psi) {
  site x(chi.lattice());
  forallsites(x) {
    chi(x,0)=Pleft*psi(x);
    chi(x,chi.L5-1)=Pright*psi(x);
  }
}



// //////////////////////////////////////////////
// choice of the default action
// //////////////////////////////////////////////
/// Pointer to the current dwfermi action
void (*default_dwfermi_action)(dwfermi_field &,
			     dwfermi_field &,
			     gauge_field &, 
			     coefficients &) = DWFermiActionFast::mul_Q;

/// Executes the current dwfermi action
void mul_Q(dwfermi_field &psi_out,
	   dwfermi_field &psi_in,
	   gauge_field &U, 
	   coefficients &coeff) {
  (*default_dwfermi_action)(psi_out, psi_in, U, coeff);
};


/// Pointer to the current dwfermi inverter
inversion_stats (*default_dwfermi_inverter)(dwfermi_field &, 
					    dwfermi_field &, 
					    gauge_field &,
					    coefficients &,
					    mdp_real, mdp_real,int)=
	&(MinRes::inverter<dwfermi_field,gauge_field>);


/// Execute the default dwfermi inverter
inversion_stats mul_invQ(dwfermi_field &psi_out, 
			 dwfermi_field &psi_in, 
			 gauge_field &U, 
			 coefficients &coeff, 
			 mdp_real absolute_precision=dwfermi_inversion_precision,
			 mdp_real relative_precision=0,
			 int max_steps=2000) {
  return (*default_dwfermi_inverter)(psi_out, psi_in, U, coeff, absolute_precision, relative_precision,max_steps);
};
