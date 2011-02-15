/////////////////////////////////////////////////////////////////
/// @file mdp_matrix.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains class mdp_matrix
///
/// Licensed under GPL2 license
/// Read attached license in file mdp_license.pdf
/// This file cannot be distributed without file mdp_license.pdf
//////////////////////////////////////////////////////////////////

/// @brief matrices of complex numbers
///
/// Example:
/// @verbatim
///    mdp_matrix A,B;
///    A.dimension(3,3);
///    A(0,0)=A(1,1)=A(2,2)=A(1,2)=1.0+I/2;
///    B=A+inv(A)+exp(A+5);
/// @endverbatim
class mdp_matrix {
private:
  /// MATRIXOPTIMIZE is still under development and does not work yet
#ifdef MATRIXOPTIMIZE
  map<int,deque<mdp_complex*> > stored;
#endif
  enum {FREE, HARD} flag;
  uint irows,icols,imax;
  mdp_complex *m;

  uint& rows(); 
  uint& cols();

public:
  void allocate();
  void reallocate();
  void deallocate();
  void dimension(const uint, const uint);
  mdp_matrix();
  mdp_matrix(const mdp_matrix& a);
  mdp_matrix(const uint r, const uint c);
  mdp_matrix(mdp_complex* z, const uint r, const uint c);
  virtual ~mdp_matrix();

  const mdp_matrix& operator=(const mdp_matrix& a);
  mdp_complex& operator[] (const uint i);
  const mdp_complex& operator[] (const uint i) const;

  mdp_complex& operator() (const uint i, const uint j);
  const mdp_complex& operator() (const uint i, const uint j) const;
  mdp_matrix operator() (const uint i);

  friend void prepare(mdp_matrix&); // does nothing, here for compatibility

  mdp_complex* address();
  const uint rows() const; 
  const uint cols() const;
  const uint size() const;
  uint rowmax() const; // compatibility function 
  uint colmax() const; // compatibility function

  friend mdp_matrix operator+ (const mdp_matrix& a);
  friend mdp_matrix operator- (const mdp_matrix& a);
  
  friend mdp_matrix operator+ (const mdp_matrix& a, const mdp_matrix& b);
  friend mdp_matrix operator- (const mdp_matrix& a, const mdp_matrix& b); // sse2
  friend mdp_matrix operator* (const mdp_matrix& a, const mdp_matrix& b); // sse2
  friend mdp_matrix operator/ (const mdp_matrix& a, const mdp_matrix& b);

  friend mdp_matrix operator+ (const mdp_matrix& a, mdp_complex b);
  friend mdp_matrix operator- (const mdp_matrix& a, mdp_complex b);
  friend mdp_matrix operator* (const mdp_matrix& a, mdp_complex b);
  friend mdp_matrix operator/ (const mdp_matrix& a, mdp_complex b);

  friend mdp_matrix operator+ (mdp_complex a, const mdp_matrix& b);
  friend mdp_matrix operator- (mdp_complex a, const mdp_matrix& b);
  friend mdp_matrix operator* (mdp_complex a, const mdp_matrix& b);
  friend mdp_matrix operator/ (mdp_complex a, const mdp_matrix& b);

  friend mdp_matrix operator+ (const mdp_matrix& a, mdp_real b);
  friend mdp_matrix operator- (const mdp_matrix& a, mdp_real b);
  friend mdp_matrix operator* (const mdp_matrix& a, mdp_real b);
  friend mdp_matrix operator/ (const mdp_matrix& a, mdp_real b);

  friend mdp_matrix operator+ (mdp_real a, const mdp_matrix& b);
  friend mdp_matrix operator- (mdp_real a, const mdp_matrix& b);
  friend mdp_matrix operator* (mdp_real a, const mdp_matrix& b);
  friend mdp_matrix operator/ (mdp_real a, const mdp_matrix& b);

  mdp_matrix operator+=(const mdp_matrix& a);  
  mdp_matrix operator-=(const mdp_matrix& a);
  mdp_matrix operator*=(const mdp_matrix& a);
  mdp_matrix operator/=(const mdp_matrix& a);

  mdp_matrix operator+=(mdp_complex a);
  mdp_matrix operator-=(mdp_complex a);
  mdp_matrix operator*=(mdp_complex a);
  mdp_matrix operator/=(mdp_complex a);

  mdp_matrix operator+=(mdp_real a);
  mdp_matrix operator-=(mdp_real a);
  mdp_matrix operator*=(mdp_real a);
  mdp_matrix operator/=(mdp_real a);

  void operator=(mdp_complex a);
  void operator=(mdp_real a);

  friend mdp_matrix inv(const mdp_matrix& a);
  friend mdp_matrix pow(const mdp_matrix& a, uint b);
  friend mdp_matrix exp(const mdp_matrix& a);
  friend mdp_matrix log(const mdp_matrix& a);
  friend mdp_matrix sin(const mdp_matrix& a);
  friend mdp_matrix cos(const mdp_matrix& a);

  friend mdp_matrix mdp_identity();
  friend mdp_matrix mdp_zero();

  friend mdp_real max(const mdp_matrix& a);
  friend mdp_matrix submatrix(const mdp_matrix& a, uint i, uint j);
  friend mdp_complex det(const mdp_matrix& a);
  friend mdp_complex trace(const mdp_matrix& a);
  friend mdp_matrix hermitian(const mdp_matrix& a);
  friend mdp_matrix transpose(const mdp_matrix& a);
  friend mdp_matrix conj(const mdp_matrix& a);

};

ostream& operator<< (ostream& os, const mdp_matrix& a) {
  uint i,j; // AW - unsigned added
  for(i=0; i<a.rows(); i++) {
    if(i==0) os << "[[";
    else     os << " [";
    for(j=0; j<a.cols(); j++) 
      if(j==0) os << " " << a(i,j);
      else os << ", " << a(i,j) << " ";
    if(i==(a.rows()-1)) os << "]]\n";
    else              os << "], \n";
  }
  return os;
}

// this is a compatibility function
void print(const mdp_matrix& a) {
  cout << a;
}

void mdp_matrix::dimension(const uint a, const uint b) {
  flag=FREE;
  rows()=a; 
  cols()=b; 
  imax=a*b;
  reallocate();
}

mdp_matrix::mdp_matrix() {
  flag=FREE;
  rows()=0;  
  cols()=0;  
  imax=0;   
  allocate();
}

mdp_matrix::mdp_matrix(const mdp_matrix& a) {
  flag=FREE;
  rows()=a.rows(); 
  cols()=a.cols(); 
  imax=a.imax;
  allocate();
  for(uint i=0; i<imax; i++) m[i]=a[i];
}

mdp_matrix::mdp_matrix(const uint a, const uint b) {
  flag=FREE;
  rows()=a; 
  cols()=b; 
  imax=a*b; 
  allocate();
}

mdp_matrix::mdp_matrix(mdp_complex* z, const uint a, const uint b) {
  flag=HARD;
  rows()=a; 
  cols()=b; 
  imax=a*b;
  m=z;
#if defined(USE_DOUBLE_PRECISION) && defined(MATRIX_SSE2)
  _sse_check_alignment((void*) m, 0xf);
#endif
}

mdp_matrix::~mdp_matrix() {
  if(flag!=HARD) deallocate();
}

inline void mdp_matrix::allocate() {
  if(flag==HARD) error("mdp_matrix::allocate()\nCannot allocate a HARD object");
#ifdef MATRIXOPTIMIZE
  deque<mdp_complex*> &q=stored[imax];
  if(q.size()) { 
    m=q[0]; 
    q.pop_front(); 
  } else
#endif 
    m=new mdp_complex[imax];
  if(imax!=0 && m==0) error("mdp_matrix::allocate()\nOut of memory");
  memset(m,0,imax*sizeof(mdp_complex));
}

inline void mdp_matrix::deallocate() {
  if(flag==HARD) error("mdp_matrix::deallocate()\nCannot allocate a HARD object");
#ifdef MATRIXOPTIMIZE
  stored[imax].push_front(m); 
#else
  delete[] m;
#endif
  m=0; 
}

inline void mdp_matrix::reallocate() {
  if(flag==HARD) error("mdp_matrix::reallocate()\nCannot allocate a HARD object");
  deallocate();
  allocate();
}
 
inline const mdp_matrix& mdp_matrix::operator=(const mdp_matrix& x) {
  uint i=0;
  if(rows()!=x.rows() || cols()!=x.cols()) {
    rows()=x.rows(); cols()=x.cols(); imax=x.imax;
    reallocate();
  }
  for(; i<imax; i++) m[i]=x[i];
  return *this;
}

inline mdp_complex& mdp_matrix::operator[] (uint i) {
  return m[i];
}

inline const mdp_complex& mdp_matrix::operator[] (uint i) const {
  return m[i];
}


inline mdp_complex& mdp_matrix::operator() (uint i, uint j) {
  return m[i*cols()+j];
}

inline const mdp_complex& mdp_matrix::operator() (uint i, uint j) const {
  return m[i*cols()+j];
}

inline mdp_matrix mdp_matrix::operator() (const uint i) {
  return mdp_matrix(m+i*cols(), cols(),1);
}

inline void prepare(mdp_matrix& a) {
}

inline mdp_complex* mdp_matrix::address() {
  return m;
}

inline uint& mdp_matrix::rows() {
  return irows;
}

inline uint& mdp_matrix::cols() {
  return icols;
}

inline const uint mdp_matrix::rows() const {
  return irows;
}

inline const uint mdp_matrix::cols() const {
  return icols;
}

inline const uint mdp_matrix::size() const {
  return imax;
}

inline uint mdp_matrix::rowmax() const {
  return irows;
}

inline uint mdp_matrix::colmax() const {
  return icols;
}

inline mdp_matrix operator+ (const mdp_matrix& a) {
  return a;
}


inline mdp_matrix operator- (const mdp_matrix& a) {
  mdp_matrix tmp(a.rows(), a.cols());
  uint i,j;
  for(i=0; i<a.rows(); i++)
    for(j=0; j<a.cols(); j++)
      tmp(i,j)=-a(i,j);
  return tmp;
}



inline mdp_matrix operator+ (const mdp_matrix& x, 
			     const mdp_matrix& y) {
  mdp_matrix z(x.rows(),x.cols());
#if defined(CHECK_ALL)
  if(x.rows()!=y.rows() || x.cols()!=y.cols()) 
    error("mdp_matrix::operator+()\nWrong argument size");
#endif
  for(register uint i=0; i<z.imax; i++) z[i]=x[i]+y[i];
  return z;
}

inline mdp_matrix operator- (const mdp_matrix& x, const mdp_matrix& y) {
  if(x.rows()!=y.rows() || x.cols()!=y.cols()) 
    error("mdp_matrix::operator-()\nwrong argument size");
  mdp_matrix z(x.rows(),x.cols());
  for(register uint i=0; i<z.imax; i++) z[i]=x[i]-y[i];
  return z;
}

inline mdp_matrix operator* (const mdp_matrix& x, const mdp_matrix& y) {
  register uint i, j, k;
  if(x.cols()!=y.rows()) 
    error("mdp_matrix::operator*()\nwrong argument size");
  mdp_matrix z(x.rows(),y.cols());
#if defined(MATRIX_SSE2) && defined(USE_DOUBLE_PRECISION)
  if(x.rows()==x.cols() && y.rows()==3) {
    register int i, cols()=y.cols(), cols()2=2*cols();
    _sse_su3    *u  =(_sse_su3*) x.m;   
    _sse_double *in =(_sse_double*) y.m; 
    _sse_double *out=(_sse_double*) z.m;
    for(i=0; i<cols(); i++, in++, out++) {
      _sse_double_load_123(*in, *(in+cols()), *(in+cols()2));
      _sse_double_su3_multiply(*u);
      _sse_double_store_up_123(*out, *(out+cols()), *(out+cols()2));
    }
    return z;
  }
#endif
#if defined(MATRIX_SSE2) && !defined(USE_DOUBLE_PRECISION)
  if(x.rows()==x.cols() && y.imax==3) {
    _sse_su3        *u  =(_sse_su3*) x.m;   
    _sse_su3_vector *in =(_sse_su3_vector*) y.m; 
    _sse_su3_vector *out=(_sse_su3_vector*) z.m;    
    _sse_float_pair_load(*in, *in);
    _sse_float_su3_multiply(*u);
    _sse_float_pair_store_up(*out, *out);
    return z;
  }
#endif
  for(i=0; i<x.rows(); i++)
    for(j=0; j<y.cols(); j++) {
      z(i,j)=x(i,0)*y(0,j);
      for(k=1; k<x.cols(); k++) z(i,j)+=x(i,k)*y(k,j);	
    }
  return z;
}

inline mdp_matrix operator/ (const mdp_matrix& a, 
			     const mdp_matrix& b) {
  return a*inv(b);
}
    
inline mdp_matrix operator+ (const mdp_matrix& a,  mdp_complex b) {
  if(a.cols()!=a.rows()) error("mdp_matrix::operator+(...)\nmdp_matrix is not squared");
  mdp_matrix tmp;
  register uint i;
  tmp=a;
  for(i=0; i<a.cols(); i++) tmp(i,i)+=b;
  return tmp;
}

inline mdp_matrix operator- (const mdp_matrix& a, mdp_complex b) {
  if(a.cols()!=a.rows()) error("mdp_matrix::operator-(...)\nmdp_matrix is not squared");
  mdp_matrix tmp;
  register uint i;
  tmp=a;
  for(i=0; i<a.cols(); i++) tmp(i,i)-=b;
  return tmp;
}

inline mdp_matrix operator* (const mdp_matrix& y, mdp_complex x) {
  register uint i;
  mdp_matrix z(y.rows(),y.cols());
#if defined(MATRIX_SSE2)
  if(y.rows()==3) {
    static _sse_float  factor1 ALIGN16;   
    static _sse_float  factor2 ALIGN16;   
    static _sse_double factor3 ALIGN16;   
    static _sse_double factor4 ALIGN16;   
    _sse_su3_vector *in =(_sse_su3_vector*) y.m; 
    _sse_su3_vector *out=(_sse_su3_vector*) z.m;    
#if defined(USE_DOUBLE_PRECISION)
    factor3.c1=factor3.c2=x.imag();
    factor4.c1=factor4.c2=x.real()/x.imag();
    for(i=0; i<y.cols(); i++, in++, out++) {
      _sse_double_load(*in);
      _sse_double_vector_mulc(factor3,factor4);
      _sse_double_store(*out);
    }
#else 
    factor1.c1=factor1.c2=factor1.c3=factor1.c4=x.imag();
    factor2.c1=factor2.c2=factor2.c3=factor2.c4=x.real()/x.imag();
    for(i=0; i<y.cols()-1; i+=2, in+=2, out+=2) {
      _sse_float_pair_load(*in, *(in+1));
      _sse_float_vector_mulc(factor1,factor2);
      _sse_float_pair_store(*out, *(out+1));
      
    }
    if(i==y.cols()-1) {
      _sse_float_pair_load(*in, *in);
      _sse_float_vector_mulc(factor1,factor2);
      _sse_float_pair_store(*out, *out);
    }
#endif
    return z;		    
  }
#endif
  for(i=0; i<y.imax; i++) z[i]=x*y[i];
  return z;
}

inline mdp_matrix operator/ (const mdp_matrix& a, mdp_complex b) {
  return a*(1.0/b);
}
;
inline mdp_matrix operator+ (mdp_complex b, const mdp_matrix& a) {
  return a+b;
}

inline mdp_matrix operator- (mdp_complex b, const mdp_matrix& a) {
  if(a.cols()!=a.rows()) error("mdp_matrix::operator-(...)\nmdp_matrix is not squared");
  mdp_matrix tmp(a.rows(), a.cols());
  uint i,j; 
  for(i=0; i<a.rows(); i++) {
    for(j=0; j<a.cols(); j++) tmp(i,j)=-a(i,j);
    tmp(i,i)+=b;
  }
  return tmp;
}

inline mdp_matrix operator* (mdp_complex x, const mdp_matrix& y) {
  return y*x;
}

inline mdp_matrix operator/ (mdp_complex b, const mdp_matrix& a) {
  return b*inv(a);
}

inline mdp_matrix operator+ (const mdp_matrix& a, mdp_real b) {
  if(a.cols()!=a.rows()) error("mdp_matrix::operator+(...)\nmdp_matrix is not squared");
  mdp_matrix tmp;
  uint i;
  tmp=a;
  for(i=0; i<a.cols(); i++) tmp(i,i).real()+=b;
  return tmp;  
}

inline mdp_matrix operator- (const mdp_matrix& a, mdp_real b) {
  if(a.cols()!=a.rows()) error("mdp_matrix::operator-(...)\nmdp_matrix is not squared");
  mdp_matrix tmp;
  uint i;
  tmp=a;
  for(i=0; i<a.cols(); i++) tmp(i,i).real()-=b;
  return tmp;  
}

inline mdp_matrix operator* (const mdp_matrix& y, mdp_real x) {
  register uint i;
  mdp_matrix z(y.rows(),y.cols());
#if defined(MATRIX_SSE2)
  if(y.rows()==3) {
    static _sse_float  factor1 ALIGN16;   
    static _sse_double factor2 ALIGN16;   

    _sse_su3_vector *in =(_sse_su3_vector*) y.m; 
    _sse_su3_vector *out=(_sse_su3_vector*) z.m;    
#if defined(USE_DOUBLE_PRECISION)
    factor2.c1=factor2.c2=x;
    for(i=0; i<y.cols(); i++, in++, out++) {
      _sse_double_load(*in);
      _sse_double_vector_mul(factor2);
      _sse_double_store(*out);
    }
#else     
    factor1.c1=factor1.c2=factor1.c3=factor1.c4=x;
    for(i=0; i<y.cols()-1; i+=2, in+=2, out+=2) {
      _sse_float_pair_load(*in, *(in+1));
      _sse_float_vector_mul(factor1);
      _sse_float_pair_store(*out, *(out+1));
    }
    if(i==y.cols()-1) {
      _sse_float_pair_load(*in, *in);
      _sse_float_vector_mul(factor1);
      _sse_float_pair_store(*out, *out);
    }
#endif
    return z;		    
  }
#endif
  for(i=0; i<y.imax; i++) z[i]=x*y[i];
  return z;
}

inline mdp_matrix operator/ (const mdp_matrix& a, mdp_real b) {
  return a*(1.0/b);
}

inline mdp_matrix operator+ (mdp_real b, const mdp_matrix& a) {
  return a+b;
}

inline mdp_matrix operator- (mdp_real b, const mdp_matrix& a) {
  if(a.cols()!=a.rows()) error("mdp_matrix::operator-(...)\nmdp_matrix is not squared");
  mdp_matrix tmp(a.rows(), a.cols());
  uint i,j; 
  for(i=0; i<a.rows(); i++) {
    for(j=0; j<a.cols(); j++) tmp(i,j)=-a(i,j);
    tmp(i,i)+=b;
  }
  return tmp;
}

inline mdp_matrix operator* (mdp_real a, const mdp_matrix& b) {
  return b*a;
}

inline mdp_matrix mdp_identity(uint i) {
  mdp_matrix tmp(i,i);
  register uint r,c;
  for(r=0; r<i; r++)
    for(c=0; c<i; c++)
      tmp(r,c)=(r==c)?mdp_complex(1,0):mdp_complex(0,0);
  return tmp;
}

inline mdp_matrix mdp_zero(uint i) {
  mdp_matrix tmp(i,i);
  register uint r,c;
  for(r=0; r<i; r++)
    for(c=0; c<i; c++)
      tmp(r,c)=(mdp_complex) 0;
  return tmp;
}

inline mdp_real max(const mdp_matrix& a) {
  register uint i;
  register double x=0,y=0;
  for(i=0; i<a.imax; i++) if ((y=abs(a[i]))>x) x=y;
  return x;
}

inline mdp_matrix submatrix (const mdp_matrix& a, uint i, uint j) {
#ifdef CHECK_ALL
  if (((a.rows()<2) || (a.cols()<2)) || 
      ((a.rows()-1<i) || (a.cols()-1<j))) error("submatrix(...)\nWrong dimensions in submatrix");
#endif
  mdp_matrix tmp(a.rows()-1,a.cols()-1);
  register uint r,c;
  for(r=0; r<i; r++)        
    for(c=0; c<j; c++)        tmp(r,c)  =a(r,c);
  for(r=i+1; r<a.rows(); r++) 
    for(c=0; c<j; c++)        tmp(r-1,c)=a(r,c);
  for(r=0; r<i; r++)        
    for(c=j+1; c<a.cols(); c++) tmp(r,c-1)=a(r,c);
  for(r=i+1; r<a.rows(); r++) 
    for(c=j+1; c<a.cols(); c++) tmp(r-1,c-1)=a(r,c);  
  return tmp;
}

inline mdp_complex det(const mdp_matrix& a) {
#ifdef CHECK_ALL
  if (a.rows()!=a.cols()) error("det(...)\nmdp_matrix is not squared");
#endif
  if (a.rows()==0) return 0;
  if (a.rows()==1) return a(0,0);
  register uint i,j,k,l;
  mdp_matrix A;
  A=a;
  mdp_complex tmp, pivot, x=mdp_complex(1,0);
  for(i=0; i<A.cols(); i++) {
    for(j=i; (A(i,j)==mdp_complex(0,0)) && (j<A.cols()); j++);
    if (j==A.cols()) return 0;
    if (i!=j) {
      for(k=0; k<A.rows(); k++) {
	tmp=A(k,j);
	A(k,j)=A(k,i);
	A(k,i)=tmp;
      }
      x*=-A(i,i);
    } else x*=A(i,i);
    for(k=i+1; k<A.rows(); k++) {
      pivot=A(k,i)/A(i,i);
      for(l=i; l<A.cols(); l++) A(k,l)-=pivot*A(i,l);
    }
  }
  return x;
}

inline mdp_matrix inv(const mdp_matrix& a) {
#ifdef CHECK_ALL
  if ((a.rows()!=a.cols()) || (a.rows()==0)) error("inv(...)\nmdp_matrix is not squared");
#endif
  mdp_matrix tma, tmp;
  mdp_complex x,pivot;
  register uint r,c,i,rmax;
  tma=a;
  tmp=mdp_identity(a.rows());
  for(c=0; c<a.cols(); c++) {
    rmax=c;
    pivot=tma(c,c);
    for(r=c+1; r<a.rows(); r++) 
      if(abs(tma(r,c))>abs(pivot)) {
	rmax=r;
	pivot=tma(r,c);
      }
    for(i=0; i<a.cols(); i++) {
      x=tma(rmax,i);
      tma(rmax,i)=tma(c,i);
      tma(c,i)=x/pivot;
      x=tmp(rmax,i);
      tmp(rmax,i)=tmp(c,i);
      tmp(c,i)=x/pivot;
    }
    for(r=0; r<a.rows(); r++) 
      if(r!=c) {
	pivot=tma(r,c);
	for(i=0; i<a.cols(); i++) {
	  tma(r,i)-=pivot*tma(c,i);
	  tmp(r,i)-=pivot*tmp(c,i);
	}
      }
  }
  return tmp;
}

inline mdp_matrix pow(const mdp_matrix& a, int i) {
#ifdef CHECK_ALL
  if (a.rows()!=a.cols()) error("pow(...)\nmdp_matrix is not squared");
#endif
  mdp_matrix tmp;
  tmp=mdp_identity(a.cols());
  uint j=(i<0)?-i:i;
  for(;j>0; j--) tmp=tmp*a;
  if (i<0) tmp=inv(tmp);
  return tmp;
}

inline mdp_matrix exp(const mdp_matrix& a) {
#ifdef CHECK_ALL
  if (a.rows()!=a.cols()) error("exp(...)\nmdp_matrix is not squared");
#endif
  mdp_matrix tmp;
  mdp_matrix term;
  register uint i=1;
  term=a;
  tmp=mdp_identity(a.rows());
  tmp+=a;
  do {
    term=(1./++i)*term*a;
    tmp+=term;
  } while (max(term)>mdp_precision);
  return tmp;
}


inline mdp_matrix mdp_matrix::operator+=(const mdp_matrix& a) {
  (*this)=(*this)+a;
  return *this;
}

inline mdp_matrix mdp_matrix::operator-=(const mdp_matrix& a) {
  (*this)=(*this)-a;
  return *this;
}

inline mdp_matrix mdp_matrix::operator*=(const mdp_matrix& a) {
  (*this)=(*this)*a;
  return *this;
}

inline mdp_matrix mdp_matrix::operator/=(const mdp_matrix& a) {
  (*this)=(*this)/a;
  return *this;
}


inline mdp_matrix mdp_matrix::operator+=(mdp_complex a) {
  (*this)=(*this)+a;
  return *this;
}

inline mdp_matrix mdp_matrix::operator-=(mdp_complex a) {
  (*this)=(*this)-a;
  return *this;
}

inline mdp_matrix mdp_matrix::operator*=(mdp_complex a) {
  (*this)=(*this)*a;
  return *this;
}

inline mdp_matrix mdp_matrix::operator/=(mdp_complex a) {
  (*this)=(*this)/a;
  return *this;
}

inline mdp_matrix mdp_matrix::operator+=(mdp_real a) {
  (*this)=(*this)+a;
  return *this;
}

inline mdp_matrix mdp_matrix::operator-=(mdp_real a) {
  (*this)=(*this)-a;
  return *this;
}

inline mdp_matrix mdp_matrix::operator*=(mdp_real a) {
  (*this)=(*this)*a;
  return *this;
}

inline mdp_matrix mdp_matrix::operator/=(mdp_real a) {
  (*this)=(*this)/a;
  return *this;
}

inline void mdp_matrix::operator=(mdp_complex a) {
  uint i,j;
  for(i=0; i<rows(); i++)
    for(j=0; j<cols(); j++)
      (*this)(i,j)=(i==j)?a:0;
}

inline void mdp_matrix::operator=(mdp_real a) {
  uint i,j;
  for(i=0; i<rows(); i++)
    for(j=0; j<cols(); j++)
      (*this)(i,j)=(i==j)?mdp_complex(a,0):0;
}


mdp_matrix log(const mdp_matrix& a) {
#ifdef CHECK_ALL
  if (a.rows()!=a.cols()) error("log(...)\nmdp_matrix is not squared");
#endif
  mdp_matrix tmp,b,c,t1;
  register uint i=1, j=1;
  b=mdp_identity(a.cols());   
  b=a-b;
  c=b;
  tmp=b;
  do { 
    c=c*b;
    t1=((mdp_real)(i=-i)/(j+=1))*c;
    tmp+=t1;
  } while (max(t1)>mdp_precision);
  return tmp;
}

mdp_matrix sin(const mdp_matrix& a) {
#ifdef CHECK_ALL
  if (a.rows()!=a.cols()) error("sin(...)\nmdp_matrix is not squared");
#endif
  mdp_matrix tmp, t1;
  register uint i=1;
  t1=a;
  tmp=t1;
  // pruintf("\n");
  do { 
    t1=((mdp_real) -1.0/(++i)/(++i))*t1*a*a;
    tmp+=t1;
  } while (max(t1)>mdp_precision);
  return tmp;
}

mdp_matrix cos(const mdp_matrix& a) {
#ifdef CHECK_ALL
  if (a.rows()!=a.cols()) error("cos(...)\nmdp_matrix is not squared");
#endif
  mdp_matrix tmp,t1;
  register uint i=0;
  t1=mdp_identity(a.rows());
  tmp=t1;
  do { 
    t1=((mdp_real) -1.0/(++i)/(++i))*t1*a*a;
      tmp+=t1;
  } while (max(t1)>mdp_precision);
  return tmp;
}

inline mdp_complex trace(const mdp_matrix& a) {
#ifdef CHECK_ALL
  if (a.rows()!=a.cols()) error("trace(...)\nmdp_matrix is not squared");
#endif
  register uint c;
  mdp_complex x;
  for(c=0; c<a.cols(); c++) x+=a(c,c);
  return x;
}

inline mdp_matrix transpose(const mdp_matrix& a) {
  mdp_matrix tmp(a.cols(), a.rows());
  register uint r,c;
  for(r=0; r<a.rows(); r++) 
    for(c=0; c<a.cols(); c++) 
      tmp(c,r)=a(r,c);
  return tmp;
}

inline mdp_matrix hermitian(const mdp_matrix& a) {
  mdp_matrix tmp(a.cols(), a.rows());
  
#if defined(MATRIX_SSE2) && defined(USE_DOUBLE_PRECISION)
  if(a.cols()==3 && a.rows()==3) {
    _sse_double_hermitian_su3((_sse_double*) tmp.address(),
			      (_sse_double*) a.address());
    return tmp;
  }
#endif  
  register uint r,c;
  for(r=0; r<a.rows(); r++) 
    for(c=0; c<a.cols(); c++) 
      tmp(c,r)=conj(a(r,c));
  return tmp;
}


inline mdp_matrix conj(const mdp_matrix& a) {
  mdp_matrix tmp(a.rows(), a.cols());
  register uint r,c;
  for(r=0; r<a.rows(); r++) 
    for(c=0; c<a.cols(); c++) 
      tmp(r,c)=conj(a(r,c));
  return tmp;
}


