/////////////////////////////////////////////////////////////////
/// @file fermiqcd_fermi_algorithms.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Various algorithms for Wilson fermions
///
/// Distributed under GPL2 License
/// 
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////

/// r(x,alpha,i) = Gamma5(alpha,beta) * s(x,beta,i)
void multiply_by_gamma5(fermi_field &r, fermi_field &s) {
  site x(r.lattice());
#if defined(SSE2)  && defined(USE_DOUBLE_PRECISION) && !defined(NO_SSE2_LINALG)
  if(r.nc==3) {
    _sse_spinor* a;
    _sse_spinor* b;
    _sse_spinor c; // takes care of case r is s
    forallsites(x) {
      a=(_sse_spinor*) &s(x,0,0);      
      b=(_sse_spinor*) &r(x,0,0);
      
      _sse_double_load_up((*a).c1);
      _sse_double_vector_minus_i_mul();
      _sse_double_store_up(c.c3);

      _sse_double_load_up((*a).c3);
      _sse_double_vector_i_mul();
      _sse_double_store_up((*b).c1);
      (*b).c3=c.c3;

      _sse_double_load_up((*a).c2);
      _sse_double_vector_minus_i_mul();
      _sse_double_store_up(c.c4);

      _sse_double_load_up((*a).c4);
      _sse_double_vector_i_mul();
      _sse_double_store_up((*b).c2);
      (*b).c4=c.c4;
    }
    return;
  }
#endif
  uint i;
  mdp_complex tmp[4];
  mdp_complex c0=Gamma5_val[0];
  mdp_complex c1=Gamma5_val[1];
  mdp_complex c2=Gamma5_val[2];
  mdp_complex c3=Gamma5_val[3];
  int i0=Gamma5_idx[0];
  int i1=Gamma5_idx[1];
  int i2=Gamma5_idx[2];
  int i3=Gamma5_idx[3];
  forallsites(x) {
    for(i=0; i<3; i++) {
      tmp[i0]=c0*s(x,0,i);
      tmp[i1]=c1*s(x,1,i);
      tmp[i2]=c2*s(x,2,i);
      tmp[i3]=c3*s(x,3,i);
      r(x,0,i)=tmp[0];
      r(x,1,i)=tmp[1];
      r(x,2,i)=tmp[2];
      r(x,3,i)=tmp[3];
    }
  }
}

// //////////////////////////////////////////////
// choice of the default action
// //////////////////////////////////////////////

/// Pointer to the current Wilson/Clover action

void (*default_fermi_action)(fermi_field &,
			     fermi_field &,
			     gauge_field &, 
			     coefficients &,
			     int) = FermiCloverActionFast::mul_Q;

/// Calls the current Wilson/Clover action
void mul_Q(fermi_field &psi_out,
	   fermi_field &psi_in,
	   gauge_field &U,
	   coefficients &coeff,
	   int parity=EVENODD) {
  (*default_fermi_action)(psi_out, psi_in, U, coeff, parity);
}

/// Pointer to the current Wilson/Clover inverter
inversion_stats (*default_fermi_inverter)(fermi_field &, 
					  fermi_field &, 
					  gauge_field &, 
					  coefficients &, 
					  mdp_real,mdp_real,int)=&(MinRes::inverter<fermi_field,gauge_field>);

/// Executes the current Wilson/Clover inverter
inversion_stats mul_invQ(fermi_field &psi_out, 
			 fermi_field &psi_in, 
			 gauge_field &U, 
			 coefficients &coeff, 
			 mdp_real absolute_precision=fermi_inversion_precision, 
			 mdp_real relative_precision=0,
			 int max_steps=2000) {
  return (*default_fermi_inverter)(psi_out, 
				   psi_in, 
				   U, 
				   coeff, 
				   absolute_precision, 
				   relative_precision,
				   max_steps); 
}

/// Checks that inversion is working
mdp_real check_inversion(fermi_field &phi, 
			 gauge_field &U, 
			 coefficients &coeff) {
  begin_function("check_inversion");
  fermi_field psi(phi.lattice(), phi.nc, phi.nspin);
  fermi_field chi(phi.lattice(), phi.nc, phi.nspin);
  site x(phi.lattice());
  int a,i;
  mdp_real precision=0;
  mul_Q(psi,phi,U,coeff);
  psi.update();
  mul_invQ(chi,psi,U,coeff);
  forallsites(x)
    for(a=0; a<phi.nspin; a++)
	for(i=0; i<phi.nc; i++) 
	  precision+=real(pow(phi(x,a,i)-chi(x,a,i),2));
  mpi.add(precision);
  precision/=phi.lattice().global_volume()*phi.nc*phi.nspin;
  mdp << "Inversion precision=" << precision << '\n';
  begin_function("end_inversion");
  return precision;
}

