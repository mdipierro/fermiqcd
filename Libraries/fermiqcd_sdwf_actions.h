/////////////////////////////////////////////////////////////////
/// @file fermiqcd_sdwf_actions.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// WORK IN PROGRESS
///
/// Distributed under GPL2 License
/// 
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////

/// @brief domain wall staggered (WORK IN PROGRESS)
class SDWFActionSlow {
 public:
  static void mul_Q(sdwf_field &chi_out,
		    sdwf_field &chi_in,
		    gauge_field &U,
		    coefficients &coeff,
		    int parity=EVENODD) {
    int   nc=U.nc;
    int   ndim=U.ndim;
    int   L5=chi_in.L5;
    int   x5;
    
    mdp_real mass; 
    mdp_real m_5=0.5;       
    mdp_real m_f=10.0;      
    mdp_real sign;

    if(coeff.has_key("m_0")) mass=coeff["m_0"];
    else error("coefficient m_0 undeclared");
    if(coeff.has_key("m_f")) m_f=coeff["m_f"];
    else error("coefficient m_f undeclared");
    if(coeff.has_key("m_5")) m_5=coeff["m_5"];
    else error("coefficient m_5 undeclared");
    if(coeff.has_key("sign")) sign=coeff["sign"];
    else sign=1;

    mdp_real two_mass=2.0*(m_f-m_5); 

    
    if(!U.swirls.allocated())
      error("fermiqcd_sdwf_algorithms/sdwf_mul_Q_TWO: no swirls?");
    
    int i,j,mu;
    site x(chi_in.lattice());
    site x_up(chi_in.lattice());
    site x_dw(chi_in.lattice());
    mdp_complex phase;
    mdp_complex *FU_up;
    mdp_complex *FU_dw;
    mdp_complex *Fchi_up;
    mdp_complex *Fchi_dw;
    
    for(x5=0; x5<L5; x5++) {
      if(two_mass!=0) 
	forallsitesofparity(x,parity)
	  for(i=0; i<nc; i++)
	    chi_out(x,x5,i)=two_mass*chi_in(x,x5,i);
      else
	forallsitesofparity(x,parity)
	  for(i=0; i<nc; i++)
	    chi_out(x,x5,i)=0;
      
      // version non-optimized for su3
      forallsitesofparity(x,parity) {
	for(mu=0; mu<ndim; mu++) {
	  x_up=x+mu;
	  x_dw=x-mu;
	  FU_up=&(U(x,mu,0,0));
	  FU_dw=&(U(x_dw,mu,0,0));
	  Fchi_up=&(chi_in(x_up,x5,0));
	  Fchi_dw=&(chi_in(x_dw,x5,0));
	  phase=(mdp_real) sign*chi_in.eta(x,mu);
	  for(i=0; i<nc; i++)
	    for(j=0; j<nc; j++) 
	      chi_out(x,x5,i)+=phase*
		(FU_up[nc*i+j]*Fchi_up[j]-conj(FU_dw[nc*j+i])*Fchi_dw[j]);
	}
	
	// #define TEST_NO_GAMMA5 
#ifndef TEST_NO_GAMMA5      
	if(x5>0) { 
	  // Pleft = (1-Gamma5) (x) 1
	  
	  phase=m_5;                         // 1 (x) 1
	  for(i=0; i<nc; i++)
	    chi_out(x,x5,i)+=phase*chi_in(x,x5-1,i);
	  
	  phase=chi_in.chiral_phase(x)*m_5; // Gamma5 (x) 1
	  x_dw=chi_in.chiral_shift(x);
	  FU_dw=&(U.swirls(x,0,0));
	  Fchi_dw=&(chi_in(x_dw,x5-1,0));
	  for(i=0; i<nc; i++)
	    for(j=0; j<nc; j++)
	      chi_out(x,x5,i)-=phase*FU_dw[nc*i+j]*Fchi_dw[j];
	}
	
	if(x5<L5-1) {
	  // Pright = (1+Gamma5) (x) 1
	  
	  phase=m_5;                         // 1 (x) 1
	  for(i=0; i<nc; i++)
	    chi_out(x,x5,i)+=phase*chi_in(x,x5+1,i);
	  
	  phase=chi_in.chiral_phase(x)*m_5; // + Gamma5 (x) 1
	  x_up=chi_in.chiral_shift(x);
	  FU_up=&(U.swirls(x,0,0));
	  Fchi_up=&(chi_in(x_up,x5+1,0));
	  for(i=0; i<nc; i++)
	    for(j=0; j<nc; j++) 
	      chi_out(x,x5,i)+=phase*FU_up[nc*i+j]*Fchi_up[j];
	}
#endif 
      }
    }
    
    forallsitesofparity(x,parity) {
      
      // Pleft = (1-Gamma5) (x) 1
      
      phase=mass;                         // 1 (x) 1
      for(i=0; i<nc; i++)
	chi_out(x,0,i)+=phase*chi_in(x,L5-1,i);
      
      phase=chi_in.chiral_phase(x)*mass; // Gamma5 (x) 1
      x_dw=chi_in.chiral_shift(x);
      FU_dw=&(U.swirls(x,0,0));
      Fchi_dw=&(chi_in(x_dw,L5-1,0));
      for(i=0; i<nc; i++)
	for(j=0; j<nc; j++)
	  chi_out(x,0,i)-=phase*FU_dw[nc*i+j]*Fchi_dw[j];
      
      // Pright = (1+Gamma5) (x) 1
      
      phase=mass;                // 1 (x) 1
      for(i=0; i<nc; i++)
	chi_out(x,L5-1,i)+=phase*chi_in(x,0,i);
      
      phase=chi_in.chiral_phase(x)*mass; // + Gamma5 (x) 1
      x_up=chi_in.chiral_shift(x);
      FU_up=&(U.swirls(x,0,0));
      Fchi_up=&(chi_in(x_up,0,0));
      for(i=0; i<nc; i++)
	for(j=0; j<nc; j++)
	  chi_out(x,L5-1,i)+phase*FU_up[nc*i+j]*Fchi_up[j];
    }
    
    /* NORMALIZE NOW
       double norm=0;
       forallsites(x)
       for(x5=0; x5<L5; x5++)
       for(i=0; i<nc; i++)
       norm+=pow(abs(chi_out(x,x5,i)),2);
       norm=sqrt(norm);
       forallsites(x)
       for(x5=0; x5<L5; x5++)
       for(i=0; i<nc; i++)
       chi_out(x,x5,i)/=norm;
    */    
  }

  

};
