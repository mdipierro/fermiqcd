/////////////////////////////////////////////////////////////////
/// @file fermiqcd_sdwf_algorithms.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// WORK IN PROGRESS
///
/// Distributed under GPL2 License
/// 
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////

const double MDP_SDWF_SGN=1.0;

void project(staggered_field &psi, sdwf_field &chi, gauge_field &U) {
  site x(chi.lattice());
  site y(chi.lattice());
  int L5=chi.L5;
  mdp_real phase;
  forallsites(x) {
    phase=chi.chiral_phase(x);
    y=chi.chiral_shift(x);
    psi(x)=
      0.5*chi(x,0)-(MDP_SDWF_SGN*0.5*phase)*(U.swirls(x)*chi(y,0))+
      0.5*chi(x,L5-1)+(MDP_SDWF_SGN*0.5*phase)*(U.swirls(x)*chi(y,L5-1));
  }
}

void project(staggered_field &psi, sdwf_field &chi, gauge_field &U, 
	     int sign, int L) {
  site x(chi.lattice());
  site y(chi.lattice());
  mdp_real phase;
  forallsites(x) {
    phase=chi.chiral_phase(x);
    y=chi.chiral_shift(x);
    psi(x)=
      0.5*chi(x,L)-(0.5*phase*sign)*(U.swirls(x)*chi(y,L));
  }
}

void project(sdwf_field &chi, staggered_field &psi, gauge_field &U) {
  site x(chi.lattice());
  site y(chi.lattice());
  mdp_real phase;
  forallsites(x) {
    phase=chi.chiral_phase(x);
    y=chi.chiral_shift(x);
    chi(x,0)=
      0.5*psi(x)-(MDP_SDWF_SGN*0.5*phase)*(U.swirls(x)*psi(y));
    chi(x,chi.L5-1)=
      0.5*psi(x)+(MDP_SDWF_SGN*0.5*phase)*(U.swirls(x)*psi(y));
  }
}


// ////////////////////////////////////////////////
// choice of the default sdwf action
// only one available at the moment
// ////////////////////////////////////////////////

void (*default_sdwf_action)(sdwf_field &,
			    sdwf_field &,
			    gauge_field &, 
			    coefficients &, 
			    int) = SDWFActionSlow::mul_Q;

void mul_Q(sdwf_field &psi_out,
	   sdwf_field &psi_in,
	   gauge_field &U, 
	   coefficients &coeff,
	   int parity=EVENODD) {
  (*default_sdwf_action)(psi_out, psi_in, U, coeff, parity);
}

// ////////////////////////////////////////////////
// choice of the default inversion method
// ////////////////////////////////////////////////

inversion_stats (*default_sdwf_inverter)(sdwf_field &, 
						 sdwf_field &,
						 gauge_field &, 
						 coefficients &, 
						 mdp_real, mdp_real,int) 
     =BiConjugateGradientStabilizedInverter<sdwf_field,gauge_field>;

inversion_stats mul_invQ(sdwf_field &psi_out, 
			 sdwf_field &psi_in, 
			 gauge_field &U, 
			 coefficients &coeff, 
			 mdp_real absolute_precision=sdwf_inversion_precision,
			 mdp_real relative_precision=0,
			 int max_steps=2000) {
  return (*default_sdwf_inverter)(psi_out, psi_in, U, coeff, absolute_precision, relative_precision, max_steps);
}

void compute_swirls_field(gauge_field &U) {
  int i,j,k,nc=U.nc;
  mdp_matrix A;
  U.swirls.allocate_mdp_matrix_field(U.lattice(), nc, nc);
  site x(U.lattice()), y(U.lattice());
  forallsites(x) if(x(0) % 2 == 0) {
    U.swirls(x)=0;
    // for(k=0; k<mdp_permutations(4); k++) { this gave some problems
    y=x;
    A=mdp_identity(nc);
    for(k=0; k<1; k++) {
      //  k=(int) ((float) mdp_permutations(4)*Random.plain());
      for(i=0; i<U.ndim; i++) {
	j=mdp_permutation(4,k,i);
	if(y(j)%2==0) {
	  A=A*U(y,+1,j);
	  y=y+j;
	} else {
	  A=A*U(y,-1,j);
	  y=y-j;
	}
      }
      U.swirls(x)+=A;
    }
    U.swirls(x)=U.swirls(x)/mdp_complex(1,0); // care here with 1
  }
  forallsites(x) if(x(0) % 2 == 1) {
    y=x;
    for(i=0; i<U.ndim; i++)
      if(y(i)%2==0)y=y+i;
      else         y=y-i;
    U.swirls(x)=inv(U.swirls(y));
  }
}
