/////////////////////////////////////////////////////////////////
/// @file fermiqcd_gauge_actions.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// All simple gauge actions
/// 
/// Distributed under GPL2 License
/// 
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////

/*
void mul_OUU(mdp_complex &A, mdp_complex &B, mdp_complex &C,
	     int nc=3, bool init=false) {
  nc=B.rows();
  for(int i=0; i<nc; i++)
    for(int j=0; j<nc; j++) {
      if(init)(&A)[i*nc+j]=0.0;
      for(int k=0; k<nc; k++)
	(&A)[i*nc+j]+=(&B)[i*nc+k]*(&C)[k*nc+j];      
    }
}
*/

/// @brief (unused)
class gauge_stats {};


/// @brief the Wilson Gauge Action
///
/// Example:
/// @verbatim
///    int ns=2, steps=10;
///    gauge_field U(lattice,nc);
///    coefficients gauge;
///    U.load("myfield.0000");
///    gauge["beta"]=6.0;
///    ImprovedGaugeAction::heatbath(U,gauge,steps);
///    U.save("myfield.0001");
/// @endverbatim
class WilsonGaugeAction {
 public:
  static void heatbath_SU2(mdp_prng &random,
			   mdp_real beta_eff, mdp_complex *a) {
    mdp_real e0,e1,e2,e3, dk, p0;
    mdp_real r1,r2,r3,r4;
    mdp_real a0,a1,a2,a3;
    mdp_complex u0,u1,u2,u3;
    mdp_complex v0,v1,v2,v3;
    mdp_real delta, phi, sin_alpha, sin_theta, cos_theta;
    e0=real(a[0])+real(a[3]);
    e1=imag(a[1])+imag(a[2]);
    e2=real(a[1])-real(a[2]);
    e3=imag(a[0])-imag(a[3]);
    dk=sqrt(e0*e0+e1*e1+e2*e2+e3*e3);
    p0=(dk*beta_eff);
    u0=mdp_complex(e0/dk,-e3/dk);
    u1=mdp_complex(-e2/dk,-e1/dk);
    u2=mdp_complex(e2/dk,-e1/dk);
    u3=mdp_complex(e0/dk,e3/dk);
    
    if(beta_eff<=0.0) 
      error("fermiqcd_gauge_algorithms/heatbath_SU2: beta is zero");
    
    do {
      do; while ((r1 = random.plain()) < 0.0001);
      r1 = -log(r1)/p0;
      do; while ((r2 = random.plain()) < 0.0001);
      r2 = -log(r2)/p0;
      r3 = cos(2.0*Pi*random.plain());
      r3 = r3*r3;
      delta = r2+r1*r3;
      r4=random.plain();
    } while (r4*r4 > (1.0-0.5*delta));
    a0=1.0-delta;
    cos_theta=2.0*random.plain()-1.0;
    sin_theta=sqrt(1.0-cos_theta*cos_theta);
    sin_alpha=sqrt(1-a0*a0);
    phi=2.0*Pi*random.plain();
    a1=sin_alpha*sin_theta*cos(phi);
    a2=sin_alpha*sin_theta*sin(phi);
    a3=sin_alpha*cos_theta;
    v0=mdp_complex(a0,a3);
    v1=mdp_complex(a2,a1);
    v2=mdp_complex(-a2,a1);
    v3=mdp_complex(a0,-a3);
    a[0]=v0*u0+v1*u2;
    a[1]=v0*u1+v1*u3;
    a[2]=v2*u0+v3*u2;
    a[3]=v2*u1+v3*u3;
  }

  /*
  static void contract_staples(mdp_matrix &M, gauge_field &U, mdp_site &U, int mu){
    int nc=U.nc;
    mdp_matrix B(nc,nc);
    mdp_site y=x+mu;
    mdp_site z;
    for(int nu=0; nu<U.ndim; nu++) if(nu!=mu) {
	z=y+nu;
	mul_OUH(M(0,0),U(y,nu,0,0),U(x+nu,mu,0,0));
	mul_OUH(M(0,0),U(y,nu,0,0),U(x+nu,mu,0,0));
      }
    
  }
  */
  static gauge_stats heatbath(gauge_field &U, 
			      coefficients &coeff, 
			      int n_iter=1) {
    begin_function("WilsonGaugeAction__heatbath");
    if(U.nc==1) error("fermiqcd_gauge_algorithms/heatbath(): U(1)? (use metropolis)");
    gauge_stats stats;
    mdp_real beta, zeta;
    if(coeff.has_key("beta")) beta=coeff["beta"]; 
    else error("beta undeclared");
    if(coeff.has_key("zeta")) zeta=coeff["zeta"]; 
    else zeta=1;

    int i,j,k,iter,mu,parity;
    mdp_matrix M(U.nc,U.nc);
    mdp_complex a[4], tmpUik;
    site x(U.lattice());
    double time=mpi.time();

    mdp << coeff;

    for(iter=0; iter<n_iter; iter++) 
      for(parity=0; parity<2; parity++) 
	for(mu=0; mu<U.ndim; mu++) {
	  forallsitesofparity(x,parity) {
	    for(i=0; i<U.nc-1; i++)
	      for(j=i+1; j<U.nc; j++) {
		if(zeta==1)    M=U(x,mu)*staple_H(U,x,mu);
		else           M=U(x,mu)*staple_H_unisotropic(U,x,mu,zeta);
		a[0]=M(i,i); 
		a[1]=M(i,j);
		a[2]=M(j,i);
		a[3]=M(j,j);
		heatbath_SU2(U.lattice().random(x),beta/U.nc,a);
		for(k=0; k<U.nc; k++) {
		  tmpUik=a[0]*U(x,mu,i,k)+a[1]*U(x,mu,j,k);
		  U(x,mu,j,k)=a[2]*U(x,mu,i,k)+a[3]*U(x,mu,j,k);
		  U(x,mu,i,k)=tmpUik;
		}
	      }
	  }
	  // The next command does all the communications!
	  U.update(parity, mu, U.nc*U.nc);
	}
    mdp << "\t<stats>\n\t\t<time>" << mpi.time()-time << "</time>\n\t</stats>\n";
    end_function("WilsonGaugeAction__heatbath");
    return stats;
  }
};

/// @brief the \f$ O(a^2)\f$ Improved Gauge Action
///
/// Example using the MILC improved action:
/// @verbatim
///    int ns=2, steps=10;
///    gauge_field U(lattice,nc);
///    coefficients gauge;
///    U.load("myfield.0000");
///    gauge["beta"]=6.0;
///    gauge["zeta"]=1.0; // MUST BE ONE
///    gauge["u_t"]=1.0;
///    gauge["u_s"]=1.0;
///    ImprovedGaugeAction::heatbath(U,gauge,steps,"MILC");
///    U.save("myfield.0001");
/// @endverbatim
/// Example using the Morningstar unisotropic improved action:
/// @verbatim
///    int ns=2, steps=10;
///    gauge_field U(lattice,nc);
///    coefficients gauge;
///    U.load("myfield.0000");
///    gauge["beta"]=6.0;
///    gauge["zeta"]=1.0; // CAN BE != ONE
///    gauge["u_t"]=1.0;
///    gauge["u_s"]=1.0;
///    ImprovedGaugeAction::heatbath(U,gauge,steps,"Morningstar");
///    U.save("myfield.0001");
/// @endverbatim
class ImprovedGaugeAction : public WilsonGaugeAction {
 private:
  static mdp_matrix rectangles_0i_H(gauge_field &U, site x, int mu) {
    mdp_matrix tmp(U.nc,U.nc);
    mdp_matrix b1(U.nc, U.nc);
    mdp_matrix b2(U.nc, U.nc);
    mdp_matrix b3(U.nc, U.nc);
    site y0(U.lattice());
    site y1(U.lattice());
    site y2(U.lattice());
    int nu;
    tmp=0;
    if(mu==0) {
      for(nu=1; nu<U.ndim; nu++) {
	y0=x+mu;
	y1=y0+nu;
	tmp+=U(y0,nu)*U(y1,nu)*hermitian(U(x,nu)*U(x+nu,nu)*U((x+nu)+nu,mu));
	y0=(x-nu)-nu;
	y1=y0+mu;
	y2=y1+nu;
	tmp+=hermitian(U(y0,mu)*U(y1,nu)*U(y2,nu))*U(y0,nu)*U(x-nu,nu);
      }
    } else {
      nu=0;
      y0=(x-mu)+nu;
      tmp+=U(x+mu,nu)*hermitian(U(x-mu,nu)*U(y0,mu)*U(x+nu,mu))*U(x-mu,mu);
      y0=x-mu;
      y1=y0-nu;
      y2=y1+mu;
      tmp+=hermitian(U(y1,mu)*U(y2,mu)*U(y2+mu,nu))*U(y1,nu)*U(y0,mu);
      y0=x+mu;
      y1=y0+nu;
      y2=y1-mu;
      tmp+=U(y0,mu)*U(y0+mu,nu)*hermitian(U(x,nu)*U(y2,mu)*U(y1,mu));
      y0=((x+mu)+mu)-nu;
      y1=y0-mu;
      y2=y1-mu;
      tmp+=U(x+mu,mu)*hermitian(U(y2,mu)*U(y1,mu)*U(y0,nu))*U(y2,nu);
      
    }
    prepare(tmp);
    return tmp;
  }
  
// if min_nu==0 then rectangles_ij computes all 6 rectanges
  
  static mdp_matrix rectangles_ij_H(gauge_field &U, site x, int mu, int min_nu=1) {
    mdp_matrix tmp(U.nc,U.nc);
    mdp_matrix b1(U.nc, U.nc);
    mdp_matrix b2(U.nc, U.nc);
    mdp_matrix b3(U.nc, U.nc);
    site y0(U.lattice());
    site y1(U.lattice());
    site y2(U.lattice());
    int nu;
    for(nu=min_nu; nu<U.ndim; nu++) if(nu!=mu) {    
      y0=x+mu;
      y1=y0+nu;
      tmp+=U(y0,nu)*U(y1,nu)*hermitian(U(x,nu)*U(x+nu,nu)*U((x+nu)+nu,mu));
      
      y0=(x-nu)-nu;
      y1=y0+mu;
      y2=y1+nu;
      tmp+=hermitian(U(y0,mu)*U(y1,nu)*U(y2,nu))*U(y0,nu)*U(x-nu,nu);
      
      y0=(x-mu)+nu;
      tmp+=U(x+mu,nu)*hermitian(U(x-mu,nu)*U(y0,mu)*U(x+nu,mu))*U(x-mu,mu);
      
      y0=x-mu;
      y1=y0-nu;
      y2=y1+mu;
      tmp+=hermitian(U(y1,mu)*U(y2,mu)*U(y2+mu,nu))*U(y1,nu)*U(y0,mu);

      y0=x+mu;
      y1=y0+nu;
      y2=y1-mu;
      tmp+=U(y0,mu)*U(y0+mu,nu)*hermitian(U(x,nu)*U(y2,mu)*U(y1,mu));

      y0=((x+mu)+mu)-nu;
      y1=y0-mu;
      y2=y1-mu;
      tmp+=U(x+mu,mu)*hermitian(U(y2,mu)*U(y1,mu)*U(y0,nu))*U(y2,nu);
    }
    prepare(tmp);
    return tmp;
  }

  // //////////////////////////////////////////////////////
  // this is slow but should make the chair correcly ...
  // see: hep-lat/0712010
  // //////////////////////////////////////////////////////

  static mdp_matrix chair_H(gauge_field &U, site x, int mu) {
    int ndim=U.ndim;
    int nu, rho;
    mdp_matrix tmp(U.nc, U.nc);
    mdp_matrix b1(U.nc, U.nc);
    mdp_matrix b2(U.nc, U.nc);
    mdp_matrix b3(U.nc, U.nc);
    site y1(U.lattice());
    site y2(U.lattice());
    site y3(U.lattice());
    site y4(U.lattice());
    site y5(U.lattice());
    tmp=0;
    for(nu=0; nu<ndim; nu++) if(nu!=mu)
      for(rho=0; rho<ndim; rho++) if((rho!=nu) && (rho!=mu)) { 
	y1=x+mu;
	y2=y1+nu;
	y3=y2+rho;
	y4=y3-mu;
	y5=y4-nu;
	tmp+=U(y1,nu)*U(y2,rho)*hermitian(U(x,rho)*U(y5,nu)*U(y4,mu));
	y1=x+mu;
	y2=y1-nu;
	y3=y2+rho;
	y4=y3-mu;
	y5=y4+nu;
	tmp+=hermitian(U(y2,nu))*U(y2,rho)
	  *hermitian(U(y4,mu))*U(y4,nu)*hermitian(U(x,rho));
	
	y1=x+mu;
	y2=y1+nu;
	y3=y2-rho;
	y4=y3-mu;
	y5=y4-nu;
	tmp+=U(y1,nu)*hermitian(U(y5,nu)*U(y4,mu)*U(y3,rho))*U(y5,rho);
	
	y1=x+mu;
	y2=y1-nu;
	y3=y2-rho;
	y4=y3-mu;
	y5=y4+nu;
	tmp+=hermitian(U(y4,mu)*U(y3,rho)*U(y2,nu))*U(y4,nu)*U(y5,rho);

      }
    return(tmp);
  }
  
  // ////////////////////////////////////////////////////////////////////
  // new_heatbath uses an improved gauge action!
  // both isotropic (param.zeta=1) and anisotropic
  // ////////////////////////////////////////////////////////////////////

  enum { iGauge_min=4 };

  static int strange_mapping(site &x) {
    int mu, type=0;
    for(mu=0; mu<x.lattice().ndim; mu++) type+=(int) pow((float) iGauge_min,mu)*(x(mu) % iGauge_min);
    return type;
  }
  
 public:
  static gauge_stats heatbath(gauge_field &U,
			      coefficients &coeff,
			      int n_iter=1,
			      string model="MILC") { 
    
    begin_function("ImprovedGaugeAction__heatbath");

    gauge_stats stats;
    mdp_real beta, zeta, u_s, u_t;

    if(coeff.has_key("beta")) beta=coeff["beta"]; else error("beta undeclared");
    if(coeff.has_key("zeta")) zeta=coeff["zeta"]; else zeta=1;
    if(coeff.has_key("u_t"))  u_t=coeff["u_t"];   else u_t=1;
    if(coeff.has_key("u_s"))  u_s=coeff["u_s"];   else u_s=1;
    
    // if(Nproc!=1)  error("improved_heatbath() does not work in parallel!");

    if(U.ndim!=4) error("fermiqcd_gauge_algorithms/improved_heatbath(): ndim!=4 (use heatbath instead)");

    /*
      if((U.lattice().size(0) % 3!=0) || (U.lattice().size(1) % 3!=0) ||
      (U.lattice().size(2) % 3!=0) || (U.lattice().size(3) % 3!=0))
      error("lattice is not divisible by 3");
    */

    if(U.nc==1)   error("fermiqcd_gauge_algorithms/improved_heatbath(): U(1)? (use metropolis instead)");
    int nc=U.nc;
    int ndim=U.ndim;
    int i,j,k,iter,mu,type;
    mdp_matrix M;
    site x(U.lattice());
    double time=mpi.time();
    mdp_complex a[4],tmpUik;
    mdp_real alpha_s;
    mdp_real c_tp=0, c_tr=0, c_sp=0, c_sr=0, c_p=0, c_r=0, c_c=0;
    
    if(model=="Morningstar") {
      
      c_tp=4.0*zeta/(3.0*pow((double)u_s*u_t,(double)2.0));
      c_tr=-1.0*zeta/(12.0*pow((double)u_s,(double)4.0)*pow((double)u_t,(double)2.0));
      c_sp=5.0/(3.0*zeta*pow((double)u_s,(double)4.0));
      c_sr=-1.0/(12.0*zeta*pow((double)u_s,(double)6.0));
      
      c_p=1.0*pow((double)u_s,(double)-4.0);
      c_r=-0.05*pow((double)u_s,(double)-6.0);
      c_c=0;
      
    } else if(model=="MILC") {
      
      if(zeta!=1) 
	error("fermiqcd_gauge_algorithms/improved_heatbath: zeta!=1");
      
      alpha_s=-4.0*log(u_s)/3.0684;
      c_p=1.0;
      c_r=-0.05*pow((double)u_s,(double)-2.0)*(1.0+0.4805*alpha_s);;
      c_c=-1.00*pow((double)u_s,(double)-2.0)*(0.03325*alpha_s);;
      
    } else {
      mdp << "Using default non-improved action" << '\n';
      stats=WilsonGaugeAction::heatbath(U,coeff,n_iter);
      end_function("ImprovedGaugeAction__heatbath");
      return stats;
    }
    
    mdp << coeff;
    
    for(iter=0; iter<n_iter; iter++)
      for(type=0; type<(int) pow((float) iGauge_min, U.ndim); type++) {
	forallsites(x) {
	  if(strange_mapping(x)==type) {
	    for(mu=0; mu<ndim; mu++) 
	      for(i=0; i<nc-1; i++)
		for(j=i+1; j<nc; j++) { 
		  if(zeta!=1) {
		    // anisotropic San Diego action
		    if(mu==0) M=U(x,mu)*(c_tp*staple_0i_H(U,x,0)+
					 c_tr*rectangles_0i_H(U,x,0));
		    else      M=U(x,mu)*(c_sp*staple_ij_H(U,x,mu)+
					 c_tp*staple_0i_H(U,x,mu)+
					 c_sr*rectangles_ij_H(U,x,mu)+
					 c_tr*rectangles_0i_H(U,x,mu));
		  } else if(c_c==0) {
		    // isotropic San Diego action
		    M=U(x,mu)*(c_p*staple_H(U,x,mu)+
			       c_r*rectangles_ij_H(U,x,mu,0));
		  } else {
		    // fully O(a^2) improved isotropic MILC action
		    M=U(x,mu)*(c_p*staple_H(U,x,mu)+
			       c_r*rectangles_ij_H(U,x,mu,0)+
			       c_c*chair_H(U,x,mu));
		  }
		  a[0]=M(i,i); 
		  a[1]=M(i,j);
		  a[2]=M(j,i);
		  a[3]=M(j,j);
		  heatbath_SU2(U.lattice().random(x),beta/U.nc,a);
		  for(k=0; k<U.nc; k++) {
		    tmpUik=a[0]*U(x,mu,i,k)+a[1]*U(x,mu,j,k);
		    U(x,mu,j,k)=a[2]*U(x,mu,i,k)+a[3]*U(x,mu,j,k);
		    U(x,mu,i,k)=tmpUik;
		  }
		}
	  }
	}
	U.update();
      }
    mdp << "\t<stats>\n\t\t<time>" << mpi.time()-time << "</time>\n\t</stats>\n";
    end_function("ImprovedGaugeAction__heatbath");
    return stats;
  }
};


