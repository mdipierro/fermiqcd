/////////////////////////////////////////////////////////////////
/// @file mdp_prng.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Class mdp_prng (the random number generator of MDP)
///
/// Licensed under GPL2 license
/// Read attached license in file mdp_license.pdf
/// This file cannot be distributed without file mdp_license.pdf
//////////////////////////////////////////////////////////////////

/// @brief Marsaglia's random number generator (same as UKQCD)
///
/// You should not instantiate this class because:
/// - there is a global object mdp_random
/// - each field "lattice" has a parallel generator "lattice.random(x)"
/// Example:
/// @verbatim
///    // print a uniform number in (0,1)
///    cout << mdp_random.plain() << endl;
///    // print a gaussian number
///    cout << mdp_random.gaussian() << endl;
///    // print a random SU(10) matrix
///    cout << mdp_random.SU(10) << endl;
/// @endverbatim
class mdp_prng {
 private:
  float u[98];
  float c;
  float cd;
  float cm;
  int ui;
  int uj;
public:
  /// return a uniform random number in (0,1)
  inline float plain() {
    float luni;			/* local variable for Float */      
    luni = u[ui] - u[uj];
    if (luni < 0.0)
      luni += 1.0;
    u[ui] = luni;
    if (--ui == 0)
      ui = 97;
    if (--uj == 0)
      uj = 97;
    if ((c -= cd) < 0.0)
      c += cm;
    if ((luni -= c) < 0.0)
      luni += 1.0;
    return ((float) luni);
  }

  ///////////////////////////////////////////////////////////////////////
  //      initializei: this takes a single integer in the range
  //		0 <= ijkl <= 900 000 000
  //	and produces the four smaller integers needed for start. It is
  //	based on the ideas contained in the RMARIN subroutine in
  //		F. James, "A Review of Pseudorandom Number Generators",
  //			Comp. Phys. Commun. Oct 1990, p.340
  //	To reduce the modifications to the existing code, seed_uni now
  //	takes the role of a preprocessor for rstart.
  //
  //	This is useful for the parallel version of the code as James
  //	states that any integer ijkl will produce a statistically
  //	independent sequence of random numbers.
  //
  //     Very funny. 
  //     If that statement was worth anything he would have provided
  //     a proof to go with it. spb 12/12/90 
  ///////////////////////////////////////////////////////////////////////

  void initialize(mdp_int ijkl) {
    int i, j, k, l, ij, kl;
    if( (ijkl < 0) || (ijkl > 900000000) )
      error("Wrong initialization for random number generator");
    ij = ijkl/30082;
    kl = ijkl - (30082 * ij);
    i = ((ij/177) % 177) + 2;
    j = (ij % 177) + 2;
    k = ((kl/169) % 178) + 1;
    l = kl % 169;
    if( (i <= 0) || (i > 178) )
      error("Wrong initialization for random number generator");	
    if( (j <= 0) || (j > 178) )
      error("Wrong initialization for random number generator");	
    if( (k <= 0) || (k > 178) )
      error("Wrong initialization for random number generator");	
    if( (l < 0) || (l > 168) )
      error("Wrong initialization for random number generator");	
    if (i == 1 && j == 1 && k == 1)
      error("Wrong initialization for random number generator");	
    int ii, jj, m;
    float s, t;
    for (ii = 1; ii <= 97; ii++) {
      s = 0.0;
      t = 0.5;
      for (jj = 1; jj <= 24; jj++) {
	m = ((i*j % 179) * k) % 179;
	i = j;
	j = k;
	k = m;
	l = (53*l+1) % 169;
	if (l*m % 64 >= 32) s += t;
	t *= 0.5;
      }
      u[ii] = s;
    }
    c  = 362436.0   / 16777216.0;
    cd = 7654321.0  / 16777216.0;
    cm = 16777213.0 / 16777216.0;
    ui = 97; /*  There is a bug in the original Fortran version */
    uj = 33; /*  of UNI -- i and j should be SAVEd in UNI()     */
  }
  mdp_prng(mdp_int k=0) {
    if(k==0) initialize(ME);
  }
  /// returns a gaussian random number
  float gaussian(float sigma=1) {
#ifdef AIX
    static int i; // assumes i is set to zero by default
#else
    static int i=0;
#endif
    static float r;
    float x,y;
    if(i==0) {
      x=(float) sqrt(-2.0*log(plain()));
      y=(float) 2.0*Pi*plain();
      i=1;
      r=sigma*x*((float) cos(y));
      return sigma*x*((float) sin(y));
    } else {
      i=0;
      return r;
    }
  }
  /// draws a random float in (0,1) from a distribution using accept-reject
  double distribution(float (*fp)(float, void *), void *a=0) {
    float x,y;
    do {
      x=plain();
      y=plain();
    } while ((*fp)(x, a)>=y);
    return x;
  }
  /// returns a random SU(n) matrix using Cabibbo-Marinari
  mdp_matrix SU(int n) {
    mdp_matrix tmp, small;
    register int i,j;
    register float alpha, sin_alpha;
    register float a0, a1, a2, a3;
    register float phi, cos_theta, sin_theta;
    if(n==1) { 
      tmp.dimension(1,1);
      alpha=2.0*Pi*this->plain();
	tmp(0,0)=mdp_complex(cos(alpha),sin(alpha));
      return tmp;
    }
    tmp=mdp_identity(n);
    for(i=0; i<(n-1); i++)
      for(j=i+1; j<n; j++) {
        alpha=Pi*this->plain(); 
	phi=2.0*Pi*this->plain();
	cos_theta=2.0*this->plain()-1.0;
	sin_theta=sqrt(1.0-cos_theta*cos_theta);
	sin_alpha=sin(alpha);
	a0=cos(alpha);
	a1=sin_alpha*sin_theta*cos(phi);
	a2=sin_alpha*sin_theta*sin(phi);
	a3=sin_alpha*cos_theta;
	small=mdp_identity(n);
	small(i,i)=mdp_complex(a0,a3);
	small(i,j)=mdp_complex(a2,a1);
	small(j,i)=mdp_complex(-a2,a1);
	small(j,j)=mdp_complex(a0,-a3);
	tmp=small*tmp;
      }
    return tmp;
  }
  /// skip n numbers from the sequence
  void skip(int n) {
    for(int i=0; i<n; i++) this->plain();
  }
} mdp_random; /// the global random number generator
