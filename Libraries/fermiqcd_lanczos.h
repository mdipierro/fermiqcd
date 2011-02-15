/////////////////////////////////////////////////////////////////
/// @file fermiqcd_lanczos.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Lanczos routine
///
/// Distributed under GPL2 License
/// 
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////

/// @brief Lanczos algorithms
/// 
/// Example:
/// @verbatim
/// mdp_gauge U(lattice,nc);
/// fermi_field psi(lattice,nc);
/// coefficients coeff;
/// coeff["kappa"]=1.12;
/// for(int k=0; k<100; k++)
///    mdp << Lanczos::step(psi,U,coeff) << endl;
/// @endverbatim
/// return mdp_complex(alpha,eta)
template<class fieldT>
class Lanczos {
 public:
  static mdp_complex step(fieldT &psi, 
		      gauge_field &U, 
		      coefficients &coeff,
		      bool force=false,
		      bool output_check=false) { 
    begin_function("Lanczos__step");
    static bool init=true;
    static fieldT p(psi);
    static fieldT q(psi);
    static fieldT r(psi);
    static double alpha, beta;
    static site x(psi.lattice());

    if(init || force) {
      mdp << "Initializing Lanczos vectors\n";

      // ///////////////////////////
      // initialize Lanczos vectors
      // ///////////////////////////

      psi.update();
      q=psi;
      double norm = sqrt(norm_square(q));
      q /= norm;
      p = 0.0;
      r = 0.0;    
      beta = 1.0;  
      init=false;
    }
    
    // ///////////////////////////
    // here the Lanczos algorithms
    // ///////////////////////////

    r=p;
    q/=beta;
    p=q;
    r*=-beta;
    q=r;

    p.update();
    mul_Q(r,p,U,coeff);  
    multiply_by_gamma5(r,r);
    q += r; 
    alpha = real_scalar_product(p,q);
    mdp_add_scaled_field(q,-alpha,p);
    beta = sqrt(norm_square(q));

    if(output_check) {
      // this prints some data for check
      static fieldT s(psi);
      double pp, qq;
      mdp_complex pq, qr, ps;
      site x(psi.lattice());
      
      pp=norm_square(p);
      qq=norm_square(q);
      pq=p*q;
      
      mul_Q(s,q,U,coeff);
      multiply_by_gamma5(s,s);
      s.update();
      qr=q*r;
      
      mul_Q(r,p,U,coeff);
      multiply_by_gamma5(r,r);
      r.update();
      ps=p*s;
      
      mdp << "|p|=" << pp << '\n';
      mdp << "|q|=" << qq << '\n';
      mdp << "&ltp|q&gt=" << pq << '\n';
      mdp << "&ltq|Q|p&gt=" << qr << '\n';
      mdp << "&ltp|Q|q&gt^*=" << ps << '\n';
      mdp << "alpha=" << alpha << '\n';
      mdp << "beta=" << beta << '\n';

    }
    end_function("Lanczos__step");   
    return mdp_complex(alpha, beta);
  }
};















