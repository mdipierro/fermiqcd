/////////////////////////////////////////////////////////////////
/// @file fermiqcd_staggered_uml_inverter.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Various stuff for staggered fermions
///
/// Distributed under GPL2 License
/// 
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////

// ////////////////////////////////////////////////
// the MILC UML inverter
// ////////////////////////////////////////////////

/// @brief MILC staggered UML inverter (optimized bicgstab)
///
/// The algorithm is taken form hep-lat/9212007
/// This is best algorithm for staggered fermions
///
/// It inverts mul_Q(psi_out,psi_in,U,coeff) iteratively
/// @param psi_out the output field passed by reference
/// @param psi_in the input field passed by reference
/// @param U the gauge field to be passed to mul_Q
/// @param coeff the gauge parameters to be passed to mul_Q
/// @param absolute_precision the target absolute precision
/// @param relative_precision the target relative precision
/// @param max_steps the maximum number of steps
///
/// Example:
/// @verbatim
/// gauge_field U(lattice,nc);
/// staggered_field psi(lattice,nc);
/// staggered_field chi(lattice,nc);
/// coefficinets coeff;
/// coeff["kappa"]=1.12;
/// U.load("myfield");
/// psi.load("myfield_psi");
/// default_staggered_inverter=StaggeredBiCGUML::inverter;
/// default_staggered_action=StaggeredAsqtadActionFast::mul_Q;
/// mul_invQ(chi,psi,U,coeff);
/// chi.save("myfield_chi");
/// @endverbatim
class StaggeredBiCGUML {
private:
  static inversion_stats staggered_BiCG_QQh(staggered_field &psi_out, 
					    staggered_field &psi_in, 
					    gauge_field &U,
					    coefficients &coeff, 
					    mdp_real mass,
					    int parity=EVEN,
					    mdp_real absolute_precision=staggered_inversion_precision,
					    mdp_real relative_precision=0,
					    int max_steps=2000) {
    
    int i, opposite_parity=EVENODD, step=0;
    int nc=psi_in.nc;
    double beta, norm, four_mass_sq;
    double pMMp, alpha, residue, rresidue, target_residue, old_residue;
    inversion_stats stats;
    
    staggered_field r(psi_in.lattice(),nc);
    staggered_field p(psi_in.lattice(),nc);
    staggered_field t(psi_in.lattice(),nc);
    staggered_field q(psi_in.lattice(),nc);
    site x(psi_in.lattice());
    
    double time=mpi.time();
  
    if(!coeff.has_key("mass") || coeff["mass"]!=0) 
      error("coefficient mass undeclared or different from zero");
    four_mass_sq=pow(2.0*mass,2);
    
    if(parity==EVEN)    opposite_parity=ODD;
    if(parity==ODD)     opposite_parity=EVEN;
    if(parity==EVENODD) opposite_parity=EVENODD;

    // psi_out=psi_in;
    // t = M^Dagger M psi_out
    
    mdp << "BEGIN BiConjugate Gradient Inversion of Q Q^DAGGER (?) ...\n";    
    mdp << "\tstep\tresidue\t\ttime (sec)\n";
    mdp << "\t====\t=======\t\t==========\n";
    
    forallsitesofparity(x,parity)
      for(i=0; i<nc; i++) 
	p(x,i)=psi_out(x,i);
    
    p.update(parity);
    mul_Q(q, p, U, coeff, opposite_parity); 
    dagger(coeff);
    
    q.update(opposite_parity);
    mul_Q(t, q, U, coeff, parity);
    dagger(coeff);

    forallsitesofparity(x,parity)
      for(i=0; i<nc; i++) {
	t(x,i)+=four_mass_sq*p(x,i);
	r(x,i)=p(x,i)=psi_in(x,i)-t(x,i);
      }
    
    residue=sqrt(real_scalar_product(r,r,parity));
    norm   =sqrt(real_scalar_product(psi_in,psi_in,parity));
    
    target_residue=absolute_precision*norm;
    
    do {
      p.update(parity);
      mul_Q(q, p, U, coeff, opposite_parity);
      dagger(coeff);
      
      q.update(opposite_parity);
      mul_Q(t, q, U, coeff, parity);
      dagger(coeff);
      
      if(four_mass_sq!=0)
	forallsitesofparity(x,parity)
	  for(i=0; i<nc; i++)
	    t(x,i)+=four_mass_sq*p(x,i);
      
      pMMp=real_scalar_product(p,t,parity);
      
      alpha=pow(residue,2)/pMMp; 
      
      forallsitesofparity(x,parity)
	for(i=0; i<nc; i++) {
	  r(x,i)-=alpha*t(x,i);
	  psi_out(x,i)+=alpha*p(x,i);
      }
      
      old_residue=residue;
      residue=sqrt(real_scalar_product(r,r,parity));
      if(residue<target_residue) break;    
      rresidue=relative_residue(r,psi_out,parity);
      if(rresidue<relative_precision) break;    				
      beta=pow(residue/old_residue,2);
      
      forallsitesofparity(x,parity)
	for(i=0; i<nc; i++)
	  p(x,i)=r(x,i)+beta*p(x,i);
      
      mdp << step << "\t" << residue << "\t" << mpi.time()-time << '\n';

      if((ME==0) && (step>100) && (residue==old_residue))
	error("fermiqcd_staggered_uml_inverter/staggered_BiCG_QQh: not converging"); 
      step++;
  } while(step<max_steps);
    mdp << "UML BiCG converged in " << step 
	<< " iterations and " << mpi.time()-time << " seconds\n",
    
    stats.target_absolute_precision=absolute_precision;
    stats.target_relative_precision=relative_precision;
    stats.max_steps=max_steps;
    stats.absolute_precision=residue;
    stats.relative_precision=rresidue;
    stats.residue=residue;
    stats.steps=step;
    stats.mul_Q_steps=step+1;
    return stats;
  }
 public:  
  static inversion_stats inverter(staggered_field &psi_out, 
				  staggered_field &psi_in, 
				  gauge_field &U,
				  coefficients &coeff, 
				  mdp_real absolute_precision=staggered_inversion_precision,
				  mdp_real relative_precision=0,
				  int max_steps=2000) {

    site x(U.lattice());
    int i;
    staggered_field r(psi_in.lattice(),U.nc);
    mdp_real mass;
    inversion_stats stats;
    if(coeff.has_key("mass")) mass=coeff["mass"];
    else error("coefficient mass undeclared");
    coeff["mass"]=0;
    
    // set masses to zero (becasue use mul_q for d_slash)
    
    // It is important to initilize the output here
    // because staggered_BiCG_QQh uses it.
    forallsites(x)
    for(i=0; i<U.nc; i++) 
      psi_out(x,i)=0;
    psi_out.update();
    
    mul_Q(r,psi_in,U,coeff,EVEN);
    
    forallsitesofparity(x,EVEN) 
      for(i=0; i<U.nc; i++) 
	r(x,i)=-r(x,i)+2.0*mass*psi_in(x,i);
    
    stats=staggered_BiCG_QQh(psi_out,r,U,coeff,mass,EVEN,absolute_precision, relative_precision,max_steps);
    psi_out.update();
    mul_Q(r,psi_out,U,coeff,ODD);
    forallsitesofparity(x,ODD) 
      for(i=0; i<U.nc; i++) 
	psi_out(x,i)=1.0/(2.0*mass)*(psi_in(x,i)-r(x,i));
    
    // restore masses
    coeff["mass"]=mass;
    stats.mul_Q_steps++;
    return stats;
  }
};
