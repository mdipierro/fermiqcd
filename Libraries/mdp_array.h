/////////////////////////////////////////////////////////////////
/// @file mdp_array.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains declaration of class mdp_array
///
/// Licensed under GPL2 license
/// Read attached license in file mdp_license.pdf
/// This file cannot be distributed without file mdp_license.pdf
//////////////////////////////////////////////////////////////////

/// @brief generic container for multidimensional arrays
///
/// Example:
/// @verbatim
///    mdp_array<float,3> a(5,5,5);
///    a(0,0,0)=3.15;
/// @endverbatim
template <class T, uint nc_>
class mdp_array {
 private:
  enum {FREE, HARD} flag;
  enum {IMAX=5};
  T    *m;        // data 
  uint  nc, c[IMAX];
  uint imax;      // total n. of elements
 public:

  inline const int& ndim() const {
    return nc;
  }

  inline T *address() {
    return m;
  }  

  inline uint* size_address() {
    return c;
  }  

  inline T& operator[] (const uint i) {
#if defined(CHECK_ALL)
    if(i>=size()) error("mdp_array::operator[]\nIndex out of bounds");
#endif
    return m[i];
  }

  inline const T& operator[] (const uint i) const {
#if defined(CHECK_ALL)
    if(i>=size()) error("mdp_array::operator[]\nIndex out of bounds");
#endif
    return m[i];
  }

  inline uint length(const uint i) const {
    return size(i);
  }

  inline uint length() const {
    return imax;
  }

  inline uint size(uint i) const {
#if defined(CHECK_ALL)
    if(i>=size()) error("mdp_array::size(...)\nIndex out of bounds");
#endif
    return c[i];
  }

  inline uint size() const {
    return imax;
  }

  inline void dimension(const uint* p) {
    if(flag!=FREE) error("mdp_array::dimension(...)\nCannot redimension a HARD object");
    register uint i;
    for(imax=1, i=0; i<nc; i++) imax*=(c[i]=p[i]); 
    for(i=nc; i<IMAX; i++) c[i]=1; 
    if(m!=0) delete[] m;
    m=new T[imax];
  }

  inline void dimension(const uint c0_=1, const uint c1_=1, const uint c2_=1, const uint c3_=1, const uint c4_=1) {
    if(m!=0) delete[] m;
    c[0]=c0_; imax=c0_;
    c[1]=c1_; imax*=c1_; 
    c[2]=c2_; imax*=c2_;
    c[3]=c3_; imax*=c3_;
    c[4]=c4_; imax*=c4_;
    m=new T[imax];
  }

  inline mdp_array(const uint c0_=1, const uint c1_=1, const uint c2_=1, const uint c3_=1, const uint c4_=1) { 
    flag=FREE;
    nc=nc_;
    m=0;
    dimension(c0_,c1_,c2_,c3_,c4_);
  }

  inline mdp_array(const uint* p) {
    flag=FREE;
    nc=nc_;
    m=0;
    dimension(p);
  }

  inline mdp_array(const T *m0, 
		   const uint c0_=1, 
		   const uint c1_=1, 
		   const uint c2_=1, 
		   const uint c3_=1, 
		   const uint c4_=1) { 
    flag=HARD;
    nc=nc_;
    m=0;
    dimension(c0_, c1_, c2_, c3_, c4_);
  }

  inline mdp_array(const T *m0,
		   const uint* p) {
    flag=HARD;
    nc=nc_;
    m=m0;
    imax=c[0]=p[0];
    uint i;
    for(i=1; i<nc; i++) imax*=(c[i]=p[i]);
    for(i=nc; i<IMAX; i++) c[i]=1;
  }

  inline mdp_array(const mdp_array &a) {
    flag=FREE;
    uint i;
    m=0;
    nc=nc_;
    if(nc!=a.nc) error("mdp_array::mdp_array(...)\nIncompatible size()");
    if(imax!=a.imax) dimension(a.imax);
    for(i=0; i<IMAX; i++) c[i]=a.c[i];
    for(i=0; i<imax; i++) m[i]=a.m[i];
  }   

  virtual ~mdp_array() {
    nc=imax=0;
    if(flag!=HARD && m!=0) delete[] m;
  }

  inline void operator= (const mdp_array& a) {
    uint i;
    nc=a.nc;
    if(nc!=a.nc) error("mdp_array::operator=(...)\nIncompatible size()");
    if(imax!=a.imax) dimension(a.imax);
    for(i=0; i<IMAX; i++) c[i]=a.c[i];
    for(i=0; i<imax; i++) m[i]=a.m[i];
  }

  inline friend void prepare (const mdp_array &a) {
    // do nothing...
  }

  inline friend mdp_array operator+ (const mdp_array& a, const mdp_array& b) {
    mdp_array tmp(a.c);
    for(register uint i=0; i<a.imax; i++) tmp[i]=a[i]+b[i];
    return tmp;
  }

  inline friend mdp_array operator- (const mdp_array& a, const mdp_array& b) {
    mdp_array tmp(a.c);
    for(register uint i=0; i<a.imax; i++) tmp[i]=a[i]-b[i];
    return tmp;
  }

  template<class T2>
  inline friend mdp_array operator* (T2 x, const mdp_array& a) {
    mdp_array tmp(a.c);
    for(register uint i=0; i<a.imax; i++) tmp[i]=a[i]*x;
    return tmp;
  }

  // Implementation of generic unary operator
  inline friend mdp_array applytoall (const mdp_array& a, 
				      T (*fptr)(T,void*), void *x=0) {
    mdp_array tmp(a.c);
    for(register uint i=0; i<a.imax; i++) 
      tmp[i]=(*fptr)(a[i],x);
    return tmp;
  }

  // implementation of generic binary operator
  inline friend mdp_array applytoall (const mdp_array& a, const mdp_array& b, 
				       T (*fptr)(T,T,void*), void* x=0) {
    mdp_array tmp(a.c);
    for(register uint i=0; i<a.imax; i++) 
      tmp[i]=(*fptr)(a[i], b[i],x);
    return tmp;
  }

  inline T &operator() (const uint i0, 
			const uint i1=0, 
			const uint i2=0, 
			const uint i3=0, 
			const uint i4=0) {
#if defined(CHECK_ALL)
    if ((i1!=0 && nc<2) ||
	(i2!=0 && nc<3) ||
	(i3!=0 && nc<4) ||
	(i4!=0 && nc<IMAX))
      error("mdp_array::operator()(...)\nIncompatible size()");
    if(i0>=c[0]) error("mdp_array::operator()\nIndex out of bounds");
    if(i1>=c[1]) error("mdp_array::operator()\nIndex out of bounds");
    if(i2>=c[2]) error("mdp_array::operator()\nIndex out of bounds");
    if(i3>=c[3]) error("mdp_array::operator()\nIndex out of bounds");
    if(i4>=c[4]) error("mdp_array::operator()\nIndex out of bounds");
#endif 
    return m[(((i0*c[1]+i1)*c[2]+i2)*c[3]+i3)*c[4]+i4];
  }

  inline const T &operator() (const uint i0, 
			      const uint i1=0, 
			      const uint i2=0, 
			      const uint i3=0, 
			      const uint i4=0) const {
#if defined(CHECK_ALL)
    if ((i1!=0 && nc<2) ||
	(i2!=0 && nc<3) ||
	(i3!=0 && nc<4) ||
	(i4!=0 && nc<IMAX))
      error("mdp_array::operator()(...)\nIncompatible size()");
    if(i0>=c[0]) error("mdp_array::operator()\nIndex out of bounds");
    if(i1>=c[1]) error("mdp_array::operator()\nIndex out of bounds");
    if(i2>=c[2]) error("mdp_array::operator()\nIndex out of bounds");
    if(i3>=c[3]) error("mdp_array::operator()\nIndex out of bounds");
    if(i4>=c[4]) error("mdp_array::operator()\nIndex out of bounds");
#endif
    return m[(((i0*c[1]+i1)*c[2]+i2)*c[3]+i3)*c[4]+i4];
  }

  friend ostream& operator<< (ostream& os, const mdp_array& a) {
    uint i;
    os << "{";

    switch(a.ndim()) {      
    case 0:
      break;
    case 1:
      for(i=0; i<a.length(); i++) 
	if(i==0) os << " " << a[i];
	else  os << ", " << a[i];
      break;
    default:
      os << " Sorry. Option not implemented";
      break;
    }
    
    os << " }";
    return os;
  }
};
