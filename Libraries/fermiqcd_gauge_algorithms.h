/////////////////////////////////////////////////////////////////
/// @file fermiqcd_gauge_algorithms.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Various stuff for gauge field
/// 
/// Distributed under GPL2 License
/// 
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////

/// make a cold gauge configuration
void set_cold(gauge_field &U) {
  begin_function("set_cold");
  mdp << "Creating a cold gauge configuration" << '\n';
  register site x(U.lattice());
  int mu;
  forallsites(x) 
    for(mu=0; mu<U.ndim; mu++) 
      U(x,mu)=mdp_identity(U.nc);
  U.update();
  end_function("set_cold");
}

/// Make a hot gauge configuration
void set_hot(gauge_field &U) { 
  begin_function("set_hot");
  mdp << "Creating a hot gauge configuration" << '\n';
  register site x(U.lattice());
  int mu;
  forallsites(x) 
    for(mu=0; mu<U.ndim; mu++)
      U(x,mu)=U.lattice().random(x).SU(U.nc);
  U.update();
  end_function("set_hot");
}

/// Check that gauge field is unitary within precision
void check_unitarity(gauge_field &U, double precision=PRECISION) {
  begin_function("check_unitarity");
  site x(U.lattice());
  int mu;
  mdp_int how_many=0;
  for(x.idx=0; x.idx<U.lattice().nvol; x.idx++) 
    for(mu=0; mu<U.ndim; mu++)
      if(max(inv(U(x,mu))-hermitian(U(x,mu)))>precision) how_many++;
  mdp.add(how_many);
  mdp << "Non unitary links found=" << how_many << '\n';
  end_function("check_unitarity");
} 

/// Compute average plaquette on plane mu-nu
mdp_real average_plaquette(gauge_field &U,int mu,int nu) {
  double tmp=0;
  site x(U.lattice());
  // U.update();
  forallsites(x) {
    tmp+=real(trace(plaquette(U,x,mu,nu)));
  }
  mdp.add(tmp);
  return tmp/(U.lattice().nvol_gl*U.nc);
} 

/// Compute average plaquette (all planes)
mdp_real average_plaquette(gauge_field &U) {
  double tmp=0;
  mdp_site x(U.lattice());
  // U.update();
  int mu,nu;
  forallsites(x) {
    for(mu=0; mu<U.ndim-1; mu++)
      for(nu=mu+1; nu<U.ndim; nu++)
	tmp+=real(trace(plaquette(U,x,mu,nu)));
  }
  mdp.add(tmp);
  return 2.0*tmp/(U.ndim*(U.ndim-1)*U.lattice().nvol_gl*U.nc);
}

/// Given a field U compute the chromo-eletro-magntic field U.em
void compute_em_field(gauge_field &U) {
  int nc=U.nc;
  site x(U.lattice());
  // It is fine to use Nmdp_matrix even if there is twist .. how lucky!
  U.em.deallocate_memory();
  U.em.allocate_em_field(U.lattice(), U.nc);
  mdp_matrix A(nc,nc);
  mdp_matrix b1(nc,nc), b2(nc,nc);
  /* 
     Fast version of code for the clover term.
     A are the four clover leafs 
  */
  int mu,nu;
  forallsites(x) 
    for(mu=0; mu<U.ndim-1; mu++)
      for(nu=mu+1; nu<U.ndim; nu++) {
	
	A=
	  U(x,mu)*U(x+mu,nu)*hermitian(U(x,nu)*U(x+nu,mu))+
	  hermitian(U(x-nu,nu))*U(x-nu,mu)*
	  U((x-nu)+mu,nu)*hermitian(U(x,mu))+
	  hermitian(U((x-mu)-nu,nu)*U(x-mu,mu))*U((x-mu)-nu,mu)*U(x-nu,nu)+
	  U(x,nu)*hermitian(U(x-mu,nu)*U((x+nu)-mu,mu))*U(x-mu,mu);
	
	U.em(x,mu,nu)=((mdp_real) 0.125)*(A-hermitian(A));

      }
  U.em.update();
}
  

// /////////////////////////////////////////////////////
// function to compute longlinks of V and attach them to
// the handle1 of U
// /////////////////////////////////////////////////////

/// For use with asqtad staggered action
/// Given field V makes a field U.long_links where (if length==2)
/// @verbatim
/// U.long_links(x,mu)=V(x,mu)*V(x+mu,mu);
/// @endverbatim
/// or (if length==3)
/// @verbatim
/// U.long_links(x,mu)=V(x,mu)*V(x+mu,mu)*V((x+mu)+mu,mu);
/// @endverbatim
void compute_long_links(gauge_field &U, gauge_field &V, int length=2) {
  if((&(U.lattice())!=&(V.lattice())) || (U.nc!=V.nc) || (U.ndim!=V.ndim))
    error("fermiqcd_gauge_auxiliary_functions/compute_long_links: U and V are not compatible lattices");
  if(V.lattice().next_next<length) 
    error("fermiqcd_gauge_auxiliary_functions/compute_long_links: next_next is not big enough");

  U.long_links.deallocate_memory();
  U.long_links.allocate_mdp_nmatrix_field(V.lattice(),U.ndim, U.nc, U.nc);
  site x(U.lattice());
  int mu;
  if(length==2)
    forallsites(x)
      for(mu=0; mu<V.ndim; mu++)
	U.long_links(x,mu)=V(x,mu)*V(x+mu,mu);
  if(length==3)
    forallsites(x)
      for(mu=0; mu<V.ndim; mu++)
	U.long_links(x,mu)=V(x,mu)*V(x+mu,mu)*V((x+mu)+mu,mu);
  U.long_links.update();
}

// ////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////
// set phases for antiperiodic boundari conditions 
// ////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////

/// To set antiperiodic boundary conditions on in direction mu
/// @verbatim
///    gauge_field U(lattice,nc);
///    // do heatbath on U
///    set_antiperiodic_phases(U,mu,true);
///    // use quarks (will have antiperiodic boundary conditions)
///    set_antiperiodic_phases(U,mu,false);
/// @endverbatim
void set_antiperiodic_phases(gauge_field &U, int mu=0, int check=true) {
  begin_function("set_antiperiodic_phases");
  site x(U.lattice());
  int  i,j;
  if(check) 
    mdp << "Setting antiperiodic boundary conditions on mu=" << mu << '\n';
  else 
    mdp << "Removing antiperiodic boundary conditions on mu=" << mu << '\n';
  forallsites(x)
    if(x(mu)==U.lattice().size(mu)-1) {
      for(i=0; i<U.nc; i++)
	for(j=0; j<U.nc; j++) 
	  U(x,mu,i,j)*=-1;
    }
  end_function("set_antiperiodic_phases");
}


/// takes a matrix M, performs a Cabibbo-Marinari cooling
/// and returns the projected matrix
mdp_matrix project_SU(mdp_matrix M, int nstep=1) {
  int i,j,k,l,step,nc=M.rowmax();
  mdp_real e0,e1,e2,e3,dk,d;
  mdp_complex dc,u0,u1,u2,u3;
  mdp_matrix B(nc,nc);
  mdp_matrix C(nc,nc);
  mdp_matrix  S(nc,nc);

  C=M;
  // /////////////////
  // preconditioning
  // /////////////////
  for(i=0; i<nc; i++) {
    for(j=0; j<i; j++) {
      dc=0;
      for(k=0; k<nc; k++)
	dc+=conj(C(k,j))*C(k,i);
      for(k=0; k<nc; k++)
	C(k,i)-=dc*C(k,j);
    }
    d=0;
    for(k=0; k<nc; k++)
      d+=pow((double)abs(C(k,i)),(double)2.0);
    d=sqrt(d);
    for(k=0; k<nc; k++)
      C(k,i)/=d;
  }
  // ////////////////////////////
  // Cabibbo Marinari Projection
  // ////////////////////////////
  for(i=0; i<nc; i++)
    for(j=0; j<nc; j++) 
      for(k=0; k<nc; k++)
	B(i,j)+=conj(M(k,i))*C(k,j);
  for(step=0; step<nstep; step++) {
    for(i=0; i<nc-1; i++)
      for(j=i+1; j<nc; j++) {
	e0=real(B(i,i))+real(B(j,j));
	e1=imag(B(i,j))+imag(B(j,i));
	e2=real(B(i,j))-real(B(j,i));
	e3=imag(B(i,i))-imag(B(j,j));
	dk=sqrt(e0*e0+e1*e1+e2*e2+e3*e3);
        u0=mdp_complex(e0,-e3)/dk;
	u1=mdp_complex(e2,-e1)/dk;
	u2=mdp_complex(-e2,-e1)/dk;
	u3=mdp_complex(e0,e3)/dk;
	// S=C;
	for(k=0; k<nc; k++) {
	  S(k,i)=C(k,i)*u0+C(k,j)*u1;
	  S(k,j)=C(k,i)*u2+C(k,j)*u3;
       	}
	if((i==nc-2) && (j==nc-1))    
	  for(k=0; k<nc; k++)
	    for(l=0; l<nc-2; l++)
	      S(k,l)=C(k,l);  
	if((i!=nc-2) || (j!=nc-1) || (step!=nstep-1)) 
	  for(k=0; k<nc; k++) {
	    C(k,i)=B(k,i)*u0+B(k,j)*u1;
	    C(k,j)=B(k,i)*u2+B(k,j)*u3;
	    B(k,i)=C(k,i);
	    B(k,j)=C(k,j);
	    C(k,i)=S(k,i);
	    C(k,j)=S(k,j);
	  }
      }  
  }
  return S;
}

/// Takes a field U and path d of length and compute the average of
/// the path on the entire lattice. Assumes computation can be done
/// locally for each site
///
/// Example:
/// @verbatim
///   int mu=0, nu=1;
///   gauge_field U(lattice,nc);
///   int d[][2]={{+1,mu},{+1,nu},{-1,mu},{-1,nu}}
///   mdp << "plaquette=" << average_path(U,4,d) << endl;
/// @endverbatim
mdp_complex average_path(gauge_field &U, int length, int d[][2]) {
  mdp_matrix_field psi1(U.lattice(),U.nc,U.nc);
  mdp_matrix_field psi2(U.lattice(),U.nc,U.nc);
  mdp_site x(U.lattice());
  mdp_complex sum=0;
  for(int i=0; i<length; i++) {
    if(i==0) 
      forallsites(x) 
	psi1(x)=U(x,d[i][0],d[i][1]);
    else
      forallsites(x) 
	psi1(x)=psi1(x)*U(x,d[i][0],d[i][1]);
    if(i<length-1) {
      psi1.update();    
      if(d[i][0]==+1)
	forallsites(x) psi2(x)=psi1(x+d[i][1]);
      else if(d[i][0]==-1)
	forallsites(x) psi2(x)=psi1(x-d[i][1]);
      psi1=psi2;
    }
  }
  forallsites(x) sum+=trace(psi1(x));
  return sum/(U.lattice().nvol_gl*U.nc);
}

/// Takes a field U, a site x, a path d of length and compute the product
/// of links amdp_int the path starting at x. Assumes computation can be done
/// locally for each site
///
/// Example:
/// @verbatim
///   int mu=0, nu=1;
///   gauge_field U(lattice,nc);
///   int d[][2]={{+1,mu},{+1,nu},{-1,mu},{-1,nu}}
///   forallsites(x)
///      cout << "plaquette(x)=" << average_path(U,x,4,d) << endl;
/// @endverbatim
mdp_matrix build_path(gauge_field &U, site x, int length, int d[][2]) {
  int nc=U.nc;
  site y(U.lattice());
  mdp_matrix tmp(nc,nc);
  tmp=U(x,d[0][0],d[0][1]);
  if(d[0][0]<0) y=x-d[0][1];
  else          y=x+d[0][1];
  for(int i=1; i<length; i++) {
    tmp=tmp*U(y,d[i][0],d[i][1]);
    if(d[i][0]<0) y=y-d[i][1];
    else          y=y+d[i][1];
  }
  return tmp;
}

void copy_path(int length, int d[][2], int c[][2]) {
  for(int i=0; i<length; i++) {
    c[i][0]=d[i][0];
    c[i][1]=d[i][1];
  }
}

void invert_path(int mu, int length, int d[][2]) {
  for(int i=0; i<length; i++)
    if(d[i][1]==mu) d[i][0]*=-1;
}

void rotate_path(int angle, int mu, int nu, int length, int d[][2]) {
  angle=(angle+360)% 360;
  if(angle==90)
    for(int i=0; i<length; i++) {
      if(d[i][1]==mu) d[i][1]=nu;
      if(d[i][1]==nu) { d[i][0]*=-1; d[i][1]=mu;}
    }
  if(angle==180)
    for(int i=0; i<length; i++) {
      if(d[i][1]==mu) d[i][0]*=-1;
      if(d[i][1]==nu) d[i][0]*=-1;
    }
  if(angle==270)
    for(int i=0; i<length; i++) {
      if(d[i][1]==mu) {d[i][0]*=-1; d[i][1]=nu; }
      if(d[i][1]==nu) d[i][1]=mu;
    }
}
 
