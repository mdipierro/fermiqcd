/////////////////////////////////////////////////////////////////
/// @file fermiqcd_staggered_actions.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Stuff for SSE/SSE2 compile with -DSSE2
///
/// Distributed under GPL2 License
/// 
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////

/// @brief Staggered/Asqtad action (SLOW: DO NOT USE IN PRODUCTION)
///
/// Example:
/// @verbatim
/// gauge_field U(lattice,nc);
/// staggered_field psi(lattice,nc);
/// staggered_field chi(lattice,nc);
/// coefficients coeff;
/// coeff["mass"]=2.0;
/// default_staggered_action=StaggeredAsqtadActionSlow::mul_Q;
/// mul_Q(chi,psi,U,coeff);
/// @endverbatim
/// Note that mul_Q(chi,psi,U,coeff) reads \f$ \chi=(/\!\!\!D[U]+m)\psi \f$ 

class StaggeredAsqtadActionSlow {
 public:
  static void mul_Q(staggered_field &chi_out,
		    staggered_field &chi_in,
		    gauge_field &U,
		    coefficients &coeff,
		    int parity=EVENODD) {

    int   ndim=U.ndim;
    mdp_real sign, two_mass;
    if(coeff.has_key("mass")) two_mass=2.0*coeff["mass"];
    else error("coefficient mass undefined");
    if(coeff.has_key("sign")) sign=coeff["sign"];
    else sign=1;

  
    int mu;
    site x(chi_in.lattice());
    site y(chi_in.lattice());
    static mdp_matrix dslash;
        
    if(two_mass!=0) 
      forallsitesofparity(x,parity)
	chi_out(x)=two_mass*chi_in(x);
    else 
      forallsitesofparity(x,parity)
	chi_out(x)=0;
    
    if(!U.long_links.allocated()) {
      // use non-naik improved action
      forallsitesofparity(x,parity) {
	for(mu=0; mu<ndim; mu++) {
	  y=x+mu;             dslash =U(x,mu)*chi_in(y);
	  y=x-mu;             dslash-=hermitian(U(y,mu))*chi_in(y);
	  chi_out(x)+=(sign*chi_in.eta(x,mu))*dslash;
	}
      }
    } else {
      // use naik improved action
      forallsitesofparity(x,parity) {
	for(mu=0; mu<ndim; mu++) {
	  y=x+mu;             dslash =U(x,mu)*chi_in(y);
	  y=(y+mu)+mu;        dslash+=U.long_links(x,mu)*chi_in(y);
	  y=x-mu;             dslash-=hermitian(U(y,mu))*chi_in(y);
	  y=(y-mu)-mu;        dslash-=hermitian(U.long_links(y,mu))*chi_in(y);
	  chi_out(x)+=sign*chi_in.eta(x,mu)*dslash;
	}
      }
    }
  }
};

/// @brief Staggered/Asqtad action
///
/// Example:
/// @verbatim
/// gauge_field U(lattice,nc);
/// staggered_field psi(lattice,nc);
/// staggered_field chi(lattice,nc);
/// coefficients coeff;
/// coeff["mass"]=2.0;
/// default_staggered_action=StaggeredAsqtadActionFast::mul_Q;
/// mul_Q(chi,psi,U,coeff);
/// @endverbatim
class StaggeredAsqtadActionFast {
 public:
  static void mul_Q(staggered_field &chi_out,
		    staggered_field &chi_in,
		    gauge_field &U,
		    coefficients &coeff,
		    int   parity=EVENODD) {
    int   nc=U.nc;
    int   ndim=U.ndim;
    mdp_real sign, two_mass;
    if(coeff.has_key("mass")) two_mass=2.0*coeff["mass"];
    else error("coefficient mass undefined");
    if(coeff.has_key("sign")) sign=coeff["sign"];
    else sign=1;

    int i,j,mu;
    site x(chi_in.lattice());
    site x_up(chi_in.lattice());
    site x_dw(chi_in.lattice());
    mdp_complex phase;
    mdp_complex *FU_up;
    mdp_complex *FU_dw;
    mdp_complex *Fchi_up;
    mdp_complex *Fchi_dw;

    if(two_mass!=0)
      forallsitesofparity(x,parity)
	for(i=0; i<nc; i++)
	  chi_out(x,i)=two_mass*chi_in(x,i);
    else
      forallsitesofparity(x,parity)
	for(i=0; i<nc; i++)
	  chi_out(x,i)=0;
    
    if(nc==3) {
      // version optimized for su3
      forallsitesofparity(x,parity) {
	for(mu=0; mu<ndim; mu++) {
	  x_up=x+mu;
	  x_dw=x-mu;
	  FU_up=&(U(x,mu,0,0));
	  FU_dw=&(U(x_dw,mu,0,0));
	  Fchi_up=&(chi_in(x_up,0));
	  Fchi_dw=&(chi_in(x_dw,0));
	  
	  // the factors of 3 in the Naik term
	  
	  phase=(mdp_real) sign*chi_in.eta(x,mu);
	  for(i=0; i<3; i++)
	    chi_out(x,i)+=phase*
	      (FU_up[3*i+0]*Fchi_up[0]+
	       FU_up[3*i+1]*Fchi_up[1]+
	       FU_up[3*i+2]*Fchi_up[2]-
	       conj(FU_dw[3*0+i])*Fchi_dw[0]-
	       conj(FU_dw[3*1+i])*Fchi_dw[1]-
	       conj(FU_dw[3*2+i])*Fchi_dw[2]);
	  
	  // //////////////////////////////////////
	  // The follwoing is called the Naik Term
	  // ... for the Lepage improved action
	  // //////////////////////////////////////
	  // with this phase naik= -1/24*pow(u0,-2)
	  // //////////////////////////////////////
	  if(U.long_links.allocated()) {
	    x_up=(x_up+mu)+mu;
	    x_dw=(x_dw-mu)-mu;
	    FU_up=&(U.long_links(x,mu,0,0));
	    FU_dw=&(U.long_links(x_dw,mu,0,0));
	    Fchi_up=&(chi_in(x_up,0));
	    Fchi_dw=&(chi_in(x_dw,0));
	    // check the factor phase
	    for(i=0; i<3; i++)
	      chi_out(x,i)+=phase*
		(FU_up[3*i+0]*Fchi_up[0]+
		 FU_up[3*i+1]*Fchi_up[1]+
		 FU_up[3*i+2]*Fchi_up[2]-
		 conj(FU_dw[3*0+i])*Fchi_dw[0]-
		 conj(FU_dw[3*1+i])*Fchi_dw[1]-
		 conj(FU_dw[3*2+i])*Fchi_dw[2]);	  
	  }
	}
      }
    } else {
      // version non-optimized for su3
      forallsitesofparity(x,parity) {
	for(mu=0; mu<ndim; mu++) {
	  x_up=x+mu;
	  x_dw=x-mu;
	  FU_up=&(U(x,mu,0,0));
	  FU_dw=&(U(x_dw,mu,0,0));
	  Fchi_up=&(chi_in(x_up,0));
	  Fchi_dw=&(chi_in(x_dw,0));
	  phase=(mdp_real) sign*chi_in.eta(x,mu);
	  for(i=0; i<nc; i++)
	    for(j=0; j<nc; j++) 
	      chi_out(x,i)+=phase*
		(FU_up[nc*i+j]*Fchi_up[j]-conj(FU_dw[nc*j+i])*Fchi_dw[j]);
	  if(U.long_links.allocated()) {
	    x_up=(x_up+mu)+mu;
	    x_dw=(x_dw-mu)-mu;
	    FU_up=&(U.long_links(x,mu,0,0));
	    FU_dw=&(U.long_links(x_dw,mu,0,0));
	    Fchi_up=&(chi_in(x_up,0));
	    Fchi_dw=&(chi_in(x_dw,0));
	    // check the factor phase
	    for(i=0; i<nc; i++)
	      for(j=0; j<nc; j++)
		chi_out(x,i)+=phase*
		  (FU_up[3*i+j]*Fchi_up[j]-conj(FU_dw[3*j+i])*Fchi_dw[j]);
	  }
	}
      }
    }
  }
};

