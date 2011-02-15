/////////////////////////////////////////////////////////////////
/// @file mdp_complex.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains delcaration of class mdp_complex for complex numbers
///
/// Licensed under GPL2 license
/// Read attached license in file mdp_license.pdf
/// This file cannot be distributed without file mdp_license.pdf
//////////////////////////////////////////////////////////////////

#if defined(DO_NOT_USE_MDP_COMPLEX)
#include "complex.h"
typedef complex<mdp_real> mdp_complex;
#else

/// @brief portable complex numbers
///
/// Example:
/// @verbatim
///    mdp_complex x=3+5*I;
///    cout << x.read() << "," << x.imag() << endl;
///    cout << sin(x) << endl;
/// @endverbatim
class mdp_complex {
 public:
  mdp_real re;
  mdp_real im;
  
  mdp_real& real() { return re; }
  mdp_real& imag() { return im; }
  
  const mdp_real& real() const { return re; }
  const mdp_real& imag() const { return im; }

  inline mdp_complex(const mdp_real a=0.0, const mdp_real b=0.0) {
    re=a;
    im=b;
  }
  
  inline mdp_complex(const mdp_complex& c) {
    re=c.real();
    im=c.imag();
  }

  inline bool operator== (const mdp_complex& c) {
    return ((re==c.real()) && (im==c.imag()));
  }

  inline bool operator!= (const mdp_complex& c) {
    return ((re!=c.real()) || (im!=c.imag()));
  }  

  inline friend mdp_real real(const mdp_complex& c) {
    return c.real();
  }

  inline friend mdp_real imag(const mdp_complex& c) {
    return c.imag();
  }

  inline friend mdp_real abs(const mdp_complex& c) {
    return sqrt(c.real()*c.real()+c.imag()*c.imag());
  }

  inline friend mdp_real arg(const mdp_complex& c) {
    return phase(c);
  }
  
  inline friend mdp_complex pow(const mdp_complex& c, mdp_real z) {
    mdp_real r=pow((double) abs(c),(double) z),p=(mdp_real) z*phase(c);
    return mdp_complex(r*cos(p),r*sin(p));
  }

  inline friend mdp_complex sqrt(const mdp_complex& c) {
    mdp_real r=sqrt(abs(c)),p=0.5*arg(c);
    return mdp_complex(r*cos(p),r*sin(p));
  }

  inline friend mdp_complex exp(const mdp_complex& c) {
    mdp_real exp_c_re=exp(c.real());
    return mdp_complex(exp_c_re*cos(c.imag()), exp_c_re*sin(c.imag()));
  }

  inline friend mdp_complex sin(const mdp_complex& c) {
    return mdp_complex(cosh(c.imag())*sin(c.real()), sinh(c.imag())*cos(c.real()));
  }

  inline friend mdp_complex cos(const mdp_complex& c) {
    return mdp_complex(cosh(c.imag())*cos(c.real()), -sinh(c.imag())*sin(c.real()));
  }

  inline friend mdp_complex times_i(const mdp_complex& c) {
    return mdp_complex(-c.imag(), c.real());
  }

  inline friend mdp_complex times_minus_i(const mdp_complex& c) {
    return mdp_complex(c.imag(), -c.real());
  }

  inline friend mdp_complex operator- (const mdp_complex& c) {
    return mdp_complex(-c.real(), -c.imag());
  }

  inline friend mdp_complex operator+ (const mdp_complex& c) {
    return c;
  }

  inline void operator+= (const mdp_complex& c) { 
    re+=c.real();
    im+=c.imag();
  }

  inline void operator-= (const mdp_complex& c) {
    re-=c.real();
    im-=c.imag();
  }

  inline void operator*= (const mdp_complex& c);

  inline void operator/= (const mdp_complex& c);

  inline friend mdp_real phase(const mdp_complex& c) {
    return atan2(c.imag(), c.real());
  }

  inline friend mdp_complex conj(const mdp_complex& a) {
    return mdp_complex(a.real(), -a.imag());
  }
  
  inline void operator+= (const mdp_real c) {
    re+=c;
  }
  
  inline void operator-= (const mdp_real c) {
    re-=c;
  }
  
  inline void operator*= (const mdp_real c) {
    re*=c;
    im*=c;
  }
  
  inline void operator/= (const mdp_real c) {
    re/=c;
    im/=c;
  }

};


inline mdp_complex operator+(const mdp_complex& a, const mdp_complex& b) {
  return mdp_complex(a.real()+b.real(), a.imag()+b.imag());
}

inline mdp_complex operator-(const mdp_complex& a, const mdp_complex& b) {
  return mdp_complex(a.real()-b.real(), a.imag()-b.imag());
}

inline mdp_complex operator*(const mdp_complex& a, const mdp_complex& b) {
  return mdp_complex(a.real()*b.real()-a.imag()*b.imag(),a.real()*b.imag()+a.imag()*b.real());
}

inline mdp_complex operator/(const mdp_complex& a, const mdp_complex& b) {
  mdp_real den=b.real()*b.real()+b.imag()*b.imag();
  return mdp_complex((a.real()*b.real()+a.imag()*b.imag())/den,(a.imag()*b.real()-a.real()*b.imag())/den);
}

inline void mdp_complex::operator*= (const mdp_complex& c) {
  (*this)=(*this)*c;
}

inline void mdp_complex::operator/= (const mdp_complex& c) {
  (*this)=(*this)/c;
}



inline mdp_complex operator+ (const mdp_complex& a, const int c) {
  return mdp_complex(a.real()+c,a.imag());
}

inline mdp_complex operator- (const mdp_complex& a, const int c) {
  return mdp_complex(a.real()-c,a.imag());
}

inline mdp_complex operator* (const mdp_complex& a, const int c) {
  return mdp_complex(a.real()*c,a.imag()*c);
}

inline mdp_complex operator/ (const mdp_complex& a, const int c) {
  return mdp_complex(a.real()/c,a.imag()/c);
  }

inline mdp_complex operator+ (const int a, const mdp_complex& c) {
  return mdp_complex(c.real()+a,c.imag());
}

inline mdp_complex operator- (const int a, const mdp_complex& c) {
  return mdp_complex(c.real()-a,c.imag());
}

inline mdp_complex operator* (const int a, const mdp_complex& c) {
  return mdp_complex(c.real()*a,c.imag()*a);
}

inline mdp_complex operator/ (const int a, const mdp_complex& c) {
  mdp_real den=(c.real()*c.real()+c.imag()*c.imag())/a;
  return mdp_complex((c.real())/den,(-c.imag())/den);
}


inline mdp_complex operator+ (const mdp_complex& a, const float c) {
  return mdp_complex(a.real()+c,a.imag());
}

inline mdp_complex operator- (const mdp_complex& a, const float c) {
  return mdp_complex(a.real()-c,a.imag());
}

inline mdp_complex operator* (const mdp_complex& a, const float c) {
  return mdp_complex(a.real()*c,a.imag()*c);
}

inline mdp_complex operator/ (const mdp_complex& a, const float c) {
  return mdp_complex(a.real()/c,a.imag()/c);
  }

inline mdp_complex operator+ (const float a, const mdp_complex& c) {
  return mdp_complex(c.real()+a,c.imag());
}

inline mdp_complex operator- (const float a, const mdp_complex& c) {
  return mdp_complex(c.real()-a,c.imag());
}

inline mdp_complex operator* (const float a, const mdp_complex& c) {
  return mdp_complex(c.real()*a,c.imag()*a);
}

inline mdp_complex operator/ (const float a, const mdp_complex& c) {
  mdp_real den=(c.real()*c.real()+c.imag()*c.imag())/a;
  return mdp_complex((c.real())/den,(-c.imag())/den);
}


inline mdp_complex operator+ (const mdp_complex& a, const double c) {
  return mdp_complex(a.real()+c,a.imag());
}

inline mdp_complex operator- (const mdp_complex& a, const double c) {
  return mdp_complex(a.real()-c,a.imag());
}

inline mdp_complex operator* (const mdp_complex& a, const double c) {
  return mdp_complex(a.real()*c,a.imag()*c);
}

inline mdp_complex operator/ (const mdp_complex& a, const double c) {
  return mdp_complex(a.real()/c,a.imag()/c);
  }

inline mdp_complex operator+ (const double a, const mdp_complex& c) {
  return mdp_complex(c.real()+a,c.imag());
}

inline mdp_complex operator- (const double a, const mdp_complex& c) {
  return mdp_complex(c.real()-a,c.imag());
}

inline mdp_complex operator* (const double a, const mdp_complex& c) {
  return mdp_complex(c.real()*a,c.imag()*a);
}

inline mdp_complex operator/ (const double a, const mdp_complex& c) {
  mdp_real den=(c.real()*c.real()+c.imag()*c.imag())/a;
  return mdp_complex((c.real())/den,(-c.imag())/den);
}

#endif


const mdp_complex I=mdp_complex(0,1);

inline mdp_real abs2(const mdp_complex &a) {
  return real(a)*real(a)+imag(a)*imag(a);
}

ostream& operator<< (ostream& os, const mdp_complex& a) {
  os << a.real() << "+" << a.imag() << "I";
  return os;
}


