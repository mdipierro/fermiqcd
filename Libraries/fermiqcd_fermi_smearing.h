/////////////////////////////////////////////////////////////////
/// @file fermiqcd_fermi_smearing.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Smearing algorithms
///
/// Distributed under GPL2 License
/// 
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////

/// @brief wupperthal smearing algotihm
///
/// Example:
/// @verbatim
/// gauge_field U(lattice,nc);
/// fermi_field psi(lattice,nc);
/// coefficient smear;
/// smear["factor"]=2;
/// smear["steps"]=10;
/// WupperthalSmearing::smear(psi,U,spear);
/// @endverbatim
class WupperthalSmearing {
 public:
  static void smear(fermi_field& psi, 
		    gauge_field& U,
		    coefficients& coeff) {

    if(coeff.has_key("factor")) 
      error("WupperthalSmearing::smear()\nCoefficient 'factor' undefined.");
    if(coeff.has_key("steps")) 
      error("WupperthalSmearing::smear()\nCoefficient 'steps' undefined.");
    mdp_real factor=coeff["factor"];
    int steps    =coeff["steps"];

    fermi_field chi(psi.lattice(),psi.nc, psi.nspin);
    site x(psi.lattice());
    int a,mu,i,j;
    for(i=0; i<steps; i++) {
      chi=psi;
      forallsites(x) 
	for(a=0; a<psi.nspin; a++)
	 for(mu=1; mu<U.ndim; mu++) 
	   psi(x,a)+=factor*(U(x,mu)*chi(x+mu,a)+U(x,-1,mu)*chi(x-mu,a));
   }
    psi.update();
  }
};


/// smears a propagator
void smearSink(fermi_propagator &S, 
	       gauge_field &U,
	       void (*smf)(fermi_field&,
			   gauge_field&,
			   coefficients&),
	       coefficients& coeff) {
  fermi_field psi(S.lattice(),S.nc, S.nspin);
  site x(psi.lattice());
  int a,b,i,j;
  for(b=0; b<S.nspin; b++)
    for(j=0; j<U.nc; j++) {
      forallsitesandcopies(x)
	for(a=0; a<S.nspin; a++)
	  for(i=0; i<U.nc; i++)
	    psi(x,a,i)=S(x,a,b,i,j);
      (*smf)(psi,U,coeff);
      forallsitesandcopies(x)
	for(a=0; a<S.nspin; a++)
	  for(i=0; i<U.nc; i++)
	    S(x,a,b,i,j)=psi(x,a,i);
    }
}
