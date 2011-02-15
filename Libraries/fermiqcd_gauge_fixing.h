/////////////////////////////////////////////////////////////////
/// @file fermiqcd_gauge_fixing.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Gauge fixing stuff
///
/// Distributed under GPL2 License
/// 
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////

/// Structure for gaugefixing stats
class gaugefixing_stats {
 public:
  uint max_steps;
  mdp_real target_precision;
  uint steps;
  mdp_real precision;
  mdp_real action;
};

/// @brief the main gaugefixing algorithm
///
/// Example:
/// @verbatim
///    gauge_field U(lattice,nc);
///    gaugefixing_stats stats;
///    U.load("myfield");
///    stats=GaugeFixing::fix(U,GaugeFixing::Coulomb,100);
///    U.save("myfield_gaugefixed");
/// @endverbatim
class GaugeFixing {
 public:
  static const int Coulomb = 0;
  static const int Landau = 10;
    
  static void hit(gauge_field &U,
		  int mu, 
		  int parity, 
		  int i, int j,
		  mdp_real overrelaxation_boost = 1) {

  // This also works for twisted boundary even if I use generic_field<> */
  
    static mdp_real a0,a1,a2,a3,b,c,d;
    static mdp_real a0_sq,ai_sq;
    static mdp_complex x0,x1;
    int k,nu,nc=U.nc;
    int opposite_parity=EVENODD;
    generic_field<mdp_complex> W(U.lattice(),4);
    mdp_matrix U_up(nc,nc), U_dw(nc,nc);
    mdp_matrix A;
    site x(U.lattice());
    site y(U.lattice());
    switch(parity) {
    case EVEN: opposite_parity=ODD; break;
    case ODD:  opposite_parity=EVEN; break;
    }
    forallsitesofparity(x,parity) {
      a0=a1=a2=a3=0;
      for(nu=0; nu<U.ndim; nu++) if(nu!=mu) {
	U_up=U(x,nu);
	U_dw=U(x-nu,nu);
	a0+=
	  +real(U_dw(i,i))+real(U_up(i,i))
	  +real(U_dw(j,j))+real(U_up(j,j));
	a1+=
	  +imag(U_dw(j,i))-imag(U_up(j,i))
	  +imag(U_dw(i,j))-imag(U_up(i,j));
	a2+=
	  -real(U_dw(j,i))+real(U_up(j,i))
	  +real(U_dw(i,j))-real(U_up(i,j));
	a3+=
	  +imag(U_dw(i,i))-imag(U_up(i,i))
	  -imag(U_dw(j,j))+imag(U_up(j,j));
      }
      ai_sq=a1*a1+a2*a2+a3*a3;
      a0_sq=a0*a0;
      b=(overrelaxation_boost*a0_sq+ai_sq)/(a0_sq+ai_sq);
      c=sqrt(a0_sq+b*b*ai_sq);
      d=b/c;
      a0/=c;
      a1*=d;
      a2*=d;
      a3*=d;
      W(x,0)=mdp_complex(a0,a3);
      W(x,1)=mdp_complex(a2,a1);
      W(x,2)=mdp_complex(-a2,a1);
      W(x,3)=mdp_complex(a0,-a3);
    }
    W.update();
    forallsitesofparity(x,parity)
      for(nu=0; nu<U.ndim; nu++) 
	for(k=0; k<U.nc; k++) {
	  x0=U(x,nu,i,k);
	  x1=U(x,nu,j,k);
	  U(x,nu,i,k)=W(x,0)*x0+W(x,1)*x1;
	  U(x,nu,j,k)=W(x,2)*x0+W(x,3)*x1;
	}
    forallsitesofparity(x,opposite_parity)
      for(nu=0; nu<U.ndim; nu++) {
	y=x+nu;
#ifndef TWISTED_BOUNDARY
	for(k=0; k<U.nc; k++) {
	  x0=U(x,nu,k,i);
	  x1=U(x,nu,k,j);
	  U(x,nu,k,i)=conj(W(y,0))*x0+conj(W(y,1))*x1;
	  U(x,nu,k,j)=conj(W(y,2))*x0+conj(W(y,3))*x1;
	}
#else
	if(in_block(y)) 
	  for(k=0; k<U.nc; k++) {
	    x0=U(x,nu,k,i);
	    x1=U(x,nu,k,j);
	    U(x,nu,k,i)=conj(W(y,0))*x0+conj(W(y,1))*x1;
	    U(x,nu,k,j)=conj(W(y,2))*x0+conj(W(y,3))*x1;
	  } else {
	    A=identity(nc);
	    A(i,i)=W(y,0);
	    A(i,j)=W(y,1);
	    A(j,i)=W(y,2);
	    A(j,j)=W(y,3);
	    twist_boundary(A,y);
	    U(x,nu)=U(x,nu)*hermitian(A);
	  }
#endif
      }
    U.update();  
  }

  static void z3_fix(gauge_field &U, int mu) {
    int     i=0,t;
    site    x(U.lattice());
    mdp_matrix  A;
    mdp_complex phase[3]= {mdp_complex(1,0), 
		       exp(2.0*Pi*I/3.0), 
		       exp(4.0*Pi*I/3.0) };
    mdp_real alpha[3];
    for(t=0; t<U.lattice().size(mu)-1; t++) {
      A=mdp_zero(U.nc);
      forallsites(x) if(x(mu)==t) A+=U(x,mu);
      mpi.add(A);
      alpha[0]=real(trace(A*conj(phase[0])));
      alpha[1]=real(trace(A*conj(phase[1])));
      alpha[2]=real(trace(A*conj(phase[2])));
      if(alpha[0]>=alpha[1] && alpha[0]>=alpha[2]) i=0;
      if(alpha[1]> alpha[0] && alpha[1]>=alpha[2]) i=1;
      if(alpha[2]> alpha[0] && alpha[2]>=alpha[1]) i=2;
      forallsites(x) 
	if(x(mu)==t) U(x,mu)*=conj(phase[i]);
	else if(x(mu)==t+1) U(x,mu)*=phase[i];
      U.update();
    }
  }
  /// performs the gauge fixing 
  /// @param U the gauge field
  /// @param mu = GaugeFixing::Coulomb or GaugeFixing::Landau or other direction
  /// @param max_steps maximum number of gaugefixing steps
  /// @param parget_precision precision in gaugefixing
  /// @param overrelaxation_boost
  /// @param z3 if set to true fixes residual Z(n) symmatry due to lattice
  ///           torus topology
  static gaugefixing_stats fix(gauge_field &U, 
		  int mu=0, 
		  int max_steps=1, 
		  mdp_real target_precision=1e-5, 
		  mdp_real overrelaxation_boost = 1,
		  bool z3=false) {
    
    gaugefixing_stats stats;
    int step,nu,parity,i,j;
    site x(U.lattice());
    double action=0;
    double precision=0;
    mdp_matrix M(U.nc,U.nc);

    stats.max_steps=max_steps;
    stats.target_precision=target_precision;

    mdp << "step\taction\tprecision\n";
    
    for(step=0; step<max_steps; step++) {
      
      for(parity=EVEN; parity<=ODD; parity++) {
	for(i=0;     i<U.nc-1; i++)
	  for(j=i+1; j<U.nc;   j++) {
	    hit(U,mu,parity,i,j,overrelaxation_boost);
	  }
      }   
      // ******** What is the gauge action ? *********
      action=0;
      precision=0;
      forallsites(x) {
	M=0;
	for(nu=0; nu<U.ndim; nu++) if(nu!=mu) {
	  M+=U(x,nu)-U(x-nu,nu);
	  action+=real(trace(U(x,nu)+U(x-nu,nu)));
	}
	M=(M-trace(M)/U.nc);
	M=M-hermitian(M);
	for(i=0; i<U.nc; i++)
	for(j=0; j<U.nc; j++)
	  precision+=(double) pow(abs(M(i,j)),2);
      }
      mpi.add(precision);
      mpi.add(action);

      precision=sqrt(precision/(U.nc*U.nc*U.lattice().global_volume()));
      action=action/(2.0*U.nc*U.lattice().global_volume());

      mdp << step << "\t" << action << "\t" << precision << '\n';
      
      if(step !=0 && precision<target_precision) break;
      // ***********************************************
      
    }
    stats.steps=step;
    stats.precision=precision;
    stats.action=action;

    if(z3) z3_fix(U,mu);

    mdp << "steps=" << step
	<< ", action=" << action
	<< ", precison+" << precision << '\n';

    return stats;
  }
};





    





