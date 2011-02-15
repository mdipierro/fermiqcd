/////////////////////////////////////////////////////////////////
/// @file fermiqcd_gauge_field.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains class gauge_field
///
/// Distributed under GPL2 License
/// 
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////

/// @brief the chromo-electr-magnetic field for any SU(n)
///
/// Example:
/// @verbatim
///    int nc=3; 
///    int box[]={10,8,8,8};
///    mdp_lattice lattice(4,box);
///    gauge_field U(lattice,nc);
///    mdp_site x(lattice);
///    U.load("myfield");
///    compute_em_field(U);
///    forallsites(x)
///      for(int mu=0; mu<U.ndim; mu++)
///        for(int nu=mu+1; nu<U.ndim; nu++)
///          cout << U.em(x,mu,nu) << endl;
/// @endverbatim
/// Note that U.em(x,mu,nu) is \f$ a^2 G_{\mu\nu} \f$ and 
/// it is a color matrix in SU(nc). \f$a\f$ is the lattice spacing.
class em_field: public mdp_complex_field {
 public:
  int ndim, nc, nem;
  em_field() {
    ndim=nc=nem=0;
    reset_field();
  }
  em_field(mdp_lattice &a, int nc_) {
    reset_field();
    ndim=a.ndim;
    nc=nc_;
    nem=(ndim*(ndim-1))/2;
    allocate_field(a, nem*nc*nc);
  }
  em_field(em_field &em) {
    reset_field();
    ndim=em.ndim;
    nc=em.nc;
    nem=(ndim*(ndim-1))/2;
    allocate_field(em.lattice(), nem*nc*nc);
  }
  void allocate_em_field(mdp_lattice &a, int nc_) {
    deallocate_field();
    ndim=a.ndim;
    nc=nc_;
    nem=(ndim*(ndim-1))/2;
    allocate_field(a, nem*nc*nc);
  }

  inline int ordered_index(int mu, int nu) const {
    // ////////////////////////
    // this map mu, nu -> k  //
    //           0   1    0  //
    //           0   2    1  //
    //           0   3    2  //
    //           1   2    3  //
    //           1   3    4  //
    //           2   3    5  //
    // It works for any ndim //
    // ////////////////////////
    if(mu<nu) return nu+(mu*(2*ndim-mu-3))/2-1;
    if(mu>nu) return mu+(nu*(2*ndim-nu-3))/2-1+nem;
    // error("wrong call to ordered_index() with mu>=nu");
    return -1; // error in this case!
  };
  
  inline mdp_matrix operator()(site x, int mu, int nu) {
#ifdef CHECK_ALL
    if(mu>=nu) error("em(x,mu,nu) for mu>=nu is not defined"); 
#endif
    register int k=ordered_index(mu,nu);
    return mdp_matrix(address(x,k*nc*nc),nc,nc);
  }
  inline mdp_complex& operator()(site x, int mu, int nu, int i, int j) {
#ifdef CHECK_ALL
    if(mu>=nu) error("em(x,mu,nu) for mu>=nu is not defined"); 
#endif
    register int k=ordered_index(mu,nu);
    return *(address(x,(k*nc+i)*nc+j));
  }
  inline const mdp_complex& operator()(site x, int mu, int nu, 
				       int i, int j) const {
#ifdef CHECK_ALL
    if(mu>=nu) error("em(x,mu,nu) for mu>=nu is not defined"); 
#endif
    register int k=ordered_index(mu,nu);
    return *(address(x,(k*nc+i)*nc+j));
  }

};

/// @brief the gauge field for any SU(n)
///
/// Example:
/// @verbatim
///    int nc=3; 
///    int box[]={10,8,8,8};
///    mdp_lattice lattice(4,box);
///    gauge_field U(lattice,nc);
///    mdp_site x(lattice);
///    // set_cold(U);
///    forallsites(x)
///       for(int mu=0; mu<U.ndim; mu++)
///          U(x,mu)=1;
///    U.update(); // synchronization
///    U.save("myfield");
///    U.load("myfield");
/// @endverbatim
/// Note that U(x,mu) is \f$ \exp{iaA_{\mu}} \f$ and 
/// it is a color matrix in SU(nc). \f$a\f$ is the lattice spacing.
class gauge_field: public mdp_complex_field {
 public:
  
  em_field           em;
  mdp_nmatrix_field  long_links;
  mdp_field<mdp_int>    i_jump;
  mdp_matrix_field   swirls;
  
  int ndim,nc;
  gauge_field() {
    reset_field();
  }
  gauge_field(const gauge_field &U) : mdp_complex_field(U) {
    ndim=U.ndim;
    nc=U.nc;
  }
  void operator=(const gauge_field &U) {
    ndim=U.ndim;
    nc=U.nc;
    mdp_complex_field::operator=(U);
  }
  gauge_field(mdp_lattice &a, int nc_) { 
    reset_field();
    ndim=a.ndim;
    nc=nc_;
    allocate_field(a, a.ndim*nc*nc);
  }
  void allocate_gauge_field(mdp_lattice &a, int nc_) {
    deallocate_field();
    ndim=a.ndim;
    nc=nc_;
    allocate_field(a, a.ndim*nc*nc);
  }
  inline mdp_matrix operator() (site x, int mu) {
#ifndef TWIST_BOUNDARY
    return mdp_matrix(address(x,mu*nc*nc),nc,nc);
#else 
    mdp_matrix tmp(address(x,mu*nc*nc),nc,nc);
    if(!in_block(x)) {
      mdp_matrix a;
      a=tmp;
#ifdef ADVANCED_TWIST
      twist_eat_links(a,x,mu);
#else
      twist_boundary(a,x);
#endif
      return(a);
    }
    return tmp;
#endif
  }

  inline const mdp_matrix operator() (site x, int mu) const {
#ifndef TWISTED_BOUNDARY
    return mdp_matrix(address(x,mu*nc*nc),nc,nc);
#else
    mdp_matrix tmp(address(x,mu*nc*nc),nc,nc);
    if(!in_block(x)) {
      mdp_matrix a;
      a=tmp;
#ifdef ADVANCED_TWIST
      twist_eat_links(a,x,mu);
#else
      twist_boundary(a,x);
#endif
      return(a);
    }
    return tmp;
#endif
  }

  inline mdp_complex& operator() (site x, int mu, int i, int j) {
#ifdef TWISTED_BOUNDARY
    if(!in_block(x)) 
      error("call to &operator=(site x) per x out of block");
#endif
    return *(address(x,(mu*nc+i)*nc+j));
  }

  inline const mdp_complex& operator() (site x, int mu, int i, int j) const {
#ifdef TWISTED_BOUNDARY
    if(!in_block(x)) 
      error("call to &operator=(site x) per x out of block");
#endif
    return *(address(x,(mu*nc+i)*nc+j));
  }

  inline mdp_matrix operator() (site x, int sign, int mu) {
    mdp_matrix tmp;
    if(sign==+1) tmp=(*this)(x,mu);
    if(sign==-1) tmp=hermitian((*this)(x-mu,mu));
    return tmp;
  }
  inline const mdp_matrix operator() (site x, int sign, int mu) const {
    mdp_matrix tmp;
    if(sign==+1) tmp=(*this)(x,mu);
    if(sign==-1) tmp=hermitian((*this)(x-mu,mu));
    return tmp;
  }


#ifndef TWISTED_BOUNDARY
  inline const mdp_complex operator() (site x, int sign, int mu, 
					int i, int j) const {
    if(sign==+1) return *(address(x,(mu*nc+i)*nc+j));
    if(sign==-1) return conj(*(address(x-mu,(mu*nc+j)*nc+i)));
    error("call to U(x,0,mu,i,j)");
    return mdp_complex(0,0);
  };
#endif

  // /////////////////////////////////////////////////////
  // /////////////////////////////////////////////////////
  // stuff for twisted boundary conditions
  // ignore if you do not care
  // /////////////////////////////////////////////////////
  // /////////////////////////////////////////////////////

#ifdef TWISTED_BOUNDARY
  inline friend void twist_boundary(mdp_matrix &M, site &x) {
    static int mu, block;
    static mdp_complex z=exp(mdp_complex(0,2.0*Pi/3.0));
    static mdp_complex a,b,c,d,e,f,g,h,i;
    
    for(mu=1; mu<x.lattice().ndim; mu++) {
      block=x.block[mu];
      if(block!=0) {
	a=M(0,0);    b=M(0,1);    c=M(0,2);
	d=M(1,0);    e=M(1,1);    f=M(1,2);
	g=M(2,0);    h=M(2,1);    i=M(2,2);
	if(mu==1 && block==1) {
	  M(0,0)=e; M(0,1)=f; M(0,2)=d;
	  M(1,0)=h; M(1,1)=i; M(1,2)=g;
	  M(2,0)=b; M(2,1)=c; M(2,2)=a;
	}      
	if(mu==1 && block==-1) {
	  M(0,0)=i; M(0,1)=g; M(0,2)=h;
	  M(1,0)=c; M(1,1)=a; M(1,2)=b;
	  M(2,0)=f; M(2,1)=d; M(2,2)=e;
	}      
	if(mu==2 && block==1) {
	  M(0,0)=a;     M(0,1)=b/z; M(0,2)=c/z/z;
	  M(1,0)=d*z;   M(1,1)=e;   M(1,2)=f/z;
	  M(2,0)=g*z*z; M(2,1)=h*z; M(2,2)=i;
	}      
	if(mu==2 && block==-1) {
	  M(0,0)=a;     M(0,1)=b*z; M(0,2)=c*z*z;
	  M(1,0)=d/z;   M(1,1)=e;   M(1,2)=f*z;
	  M(2,0)=g/z/z; M(2,1)=h/z; M(2,2)=i;
	}   
	if(mu==3 && block==1) {
	  M(0,0)=i;     M(0,1)=g/z; M(0,2)=h/z/z;
	  M(1,0)=c*z;   M(1,1)=a;   M(1,2)=b/z;
	  M(2,0)=f*z*z; M(2,1)=d*z; M(2,2)=e;
	}      
	if(mu==3 && block==-1) {
	  M(0,0)=e;     M(0,1)=f*z; M(0,2)=d*z*z;
	  M(1,0)=h/z;   M(1,1)=i;   M(1,2)=g*z;
	  M(2,0)=b/z/z; M(2,1)=c/z; M(2,2)=a;
	}
	if(block*block>1) error("two blocks out\n");
      }
    }
  }

  friend void define_twist_matrices() {
    begin_function("gauge_field__define_twist_matrices");
    mdp_complex z=exp(mdp_complex(0, 2.0*Pi/3.0));
    OmegaTwist[0]=identity(3);
    OmegaTwist[1]=zero(3);
    OmegaTwist[1](0,1)=1;
    OmegaTwist[1](1,2)=1;
    OmegaTwist[1](2,0)=1;
    OmegaTwist[2]=zero(3);
    OmegaTwist[2](0,0)=1.0/z;
    OmegaTwist[2](1,1)=1;
    OmegaTwist[2](2,2)=z;
    OmegaTwist[3]=zero(3);
    OmegaTwist[3](0,2)=1.0/z;
    OmegaTwist[3](1,0)=1;
    OmegaTwist[3](2,1)=z;
    end_function("gauge_field__define_twist_matrices");
  }
  
  inline friend void twist_eat_fields(mdp_matrix &M, site &x, 
				      gauge_field &omega) {
    begin_function("gauge_field__twist_eat_matrices");
    static int mu, block;
    static mdp_complex z=exp(mdp_complex(0,2.0*Pi/3.0));
    static mdp_complex a,b,c,d,e,f,g,h,i;
    
    for(mu=1; mu<x.lattice().ndim; mu++) {
      block=x.block[mu];
      if(block!=0) {
	if(mu==1 && block==1) 
	  M=omega(x,mu)*M*hermitian(omega(x,mu));
	if(mu==1 && block==-1) 
	  M=hermitian(omega(x,mu))*M*omega(x,mu);
	if(mu==2 && block==1) 
	  M=omega(x,mu)*M*hermitian(omega(x,mu));
	if(mu==2 && block==-1) 
	  M=hermitian(omega(x,mu))*M*omega(x,mu);
	if(mu==3 && block==1) 
	  M=omega(x,mu)*M*hermitian(omega(x,mu));
	if(mu==3 && block==-1) 
	  M=hermitian(omega(x,mu))*M*omega(x,mu);
	if(block*block>1) error("two blocks out\n");
      }
    }
    end_function("gauge_field__twist_eat_matrices");
  }
#endif
};

