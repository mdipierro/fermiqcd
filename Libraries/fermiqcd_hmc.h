// Experimental, based on code developed by Simon Catterall
// 

template<class GaugeClass, class FermiClass>
class HMC {
  
private:
  // these are all pointers here for speed tricks
  GaugeClass * to_U;
  GaugeClass * to_V; // for alternate representation!
  FermiClass * to_F;
  GaugeClass * to_f_U;
  GaugeClass * to_p_U;
  GaugeClass * to_old_U;
  GaugeClass * to_old_f_U;
  FermiClass * to_f_F;
  FermiClass * to_p_F;
  FermiClass * to_old_F;
  FermiClass * to_old_f_F;
  
public:
  static const int FUNDAMENTAL = 0;
  static const int SYMMETRIC = 1;
  static const int ANTISYMMETRIC = 2;
  static const int SO4 = 14;
  coefficients coeff;
  double bs, bs_old;
  double fs, fs_old;
  mdp_real s_old;
  int numgen;
  int accepted;
  int steps;
  vector<mdp_matrix> S;
  vector<mdp_matrix> lambda;
  
  HMC(GaugeClass &U, FermiClass &F, coefficients &coeff) {
    
    this->coeff = coeff;
    this->accepted = 0;
    this->steps = 0;

    this->to_F = &F;
    this->to_U = &U;
    this->to_V = &U; // initialize should change this

    this->to_f_U = new GaugeClass(U.lattice(),U.nc);
    this->to_p_U = new GaugeClass(U.lattice(),U.nc);
    this->to_old_U = new GaugeClass(U.lattice(),U.nc);
    this->to_old_f_U = new GaugeClass(U.lattice(),U.nc);

    this->to_f_F = new FermiClass(F.lattice(),F.nc);
    this->to_p_F = new FermiClass(F.lattice(),F.nc);
    this->to_old_F = new FermiClass(F.lattice(),F.nc);
    this->to_old_f_F = new FermiClass(F.lattice(),F.nc);

    cout << to_p_U->lattice().ndim << endl;

    initialize();
  }
  
  ~HMC() {
    delete this->to_f_U;
    delete this->to_p_U;
    delete this->to_old_U;
    delete this->to_old_f_U;
    delete this->to_f_F;
    delete this->to_p_F;
    delete this->to_old_F;
    delete this->to_old_f_F;
    if(this->to_V != this->to_U) delete this->to_V;
  }
  
  void step() {
    GaugeClass &U = *(this->to_U);
    GaugeClass &V = *(this->to_V);
    FermiClass &F = *(this->to_F);
    GaugeClass &f_U = *(this->to_f_U);
    FermiClass &f_F = *(this->to_f_F);

    GaugeClass &p_U = *(this->to_p_U);
    GaugeClass &old_U = *(this->to_old_U);
    GaugeClass &old_f_U = *(this->to_old_f_U);
    FermiClass &p_F = *(this->to_p_F);
    FermiClass &old_F = *(this->to_old_F);
    FermiClass &old_f_F = *(this->to_old_f_F);

    double h_old, k_old, k_new, s_new, h_new;
    
    // IS THIS NEEDD?
    // U.update();
    // F.update();
    // f_U.update();
    // f_F.update();

    if(steps == 0) {
      bs_old = bs; // how initializes these?
      fs_old = fs;    
      s_old = compute_action(U, V, F);
      compute_force(U, f_U, F, f_F);
    }

    compute_gaussian_momenta(p_U);
    set_gaussian(p_F);
    k_old = compute_kinetic_energy(p_U, p_F);

    /// CHECKED UP TO HERE!

    h_old = s_old + k_old;
    old_U = U;
    old_F = F;
    old_f_U = f_U;
    old_f_F = f_F;

    for(int i = 0; i<coeff["trajectory_length"]; i++) 
      compute_fields_evolution(U, p_U, f_U, F, p_F, f_F);

    s_new = compute_action(U, V, F);
    k_new = compute_kinetic_energy(p_U, p_F);
    h_new = s_new + k_new;

    cout << "h_new " << h_new << endl;

    //metropolis test
    float random_number = mdp_random.plain();
    mdp.broadcast(random_number, 0);
    if(random_number<exp(h_old - h_new)) {
      // ACCEPT
      mdp << "DH: " << abs(h_new - h_old)/h_new << endl;
      s_old = s_new;
      bs_old = bs;
      fs_old = fs;
      accepted++;
    } else {
      // REJECT
      bs = bs_old;
      fs = fs_old;
      U = old_U;
      F = old_F;
      f_F = old_f_F;
      f_U = old_f_U;
    }
    steps++;
  }

  mdp_real acceptance_rate() {
    return (mdp_real) accepted/steps;
  }

  // methods that are specific for the particular action 
  void initialize() {    
    GaugeClass &U = * to_U;    
    if(coeff["representation"] == FUNDAMENTAL) {
      SU_Generators g(U.nc);
      numgen = g.ngenerators;
      to_V = to_U; // default      
      lambda.resize(g.ngenerators);
      for(int a=0; a<lambda.size(); a++) lambda[a] = g.lambda[a];
    } else if(coeff["representation"] == SYMMETRIC) {
      SU_Generators g(U.nc);
      numgen = g.ngenerators;
      int dimrep = U.nc*(U.nc+1)/2;
      lambda.resize(g.ngenerators);
      for(int a=0; a<lambda.size(); a++) lambda[a] = g.lambda[a];
      to_V = new GaugeClass(U.lattice(), dimrep);
      mdp_matrix tmp(U.nc, U.nc);
      S.resize(g.ngenerators);

      tmp(0,0)=1.0/sqrt(2.0)*mdp_complex(1.0,0.0);
      tmp(1,1)=1.0/sqrt(2.0)*mdp_complex(1.0,0.0);
      S[0]=tmp;
      tmp(0,0)=mdp_complex(0.0,0.0);
      tmp(1,1)=mdp_complex(0.0,0.0);

      tmp(0,1)=1.0/sqrt(2.0)*mdp_complex(1.0,0.0);
      tmp(1,0)=1.0/sqrt(2.0)*mdp_complex(1.0,0.0);
      S[1]=tmp;
      tmp(0,1)=mdp_complex(0.0,0.0);
      tmp(1,0)=mdp_complex(0.0,0.0);

      tmp(0,0)=1.0/sqrt(2.0)*mdp_complex(1.0,0.0);
      tmp(1,1)=1.0/sqrt(2.0)*mdp_complex(-1.0,0.0);
      S[2]=tmp;
    } else if(coeff["representation"] == ANTISYMMETRIC) {
      throw string("ANTISYMMETRIC representation is not defined yet");   
    } else if(coeff["representation"] == SO4) {
      SO_Generators g(4);
      numgen = g.ngenerators;
      lambda.resize(g.ngenerators);
      for(int a=0; a<lambda.size(); a++) lambda[a] = g.lambda[a];
      to_V = to_U;
    }    
    /*
    ampdeg = 9.9997112279957390e - 02;
    amp[0] = 3.6832229992796258e - 07; shift[0] = 3.7549480881878877e - 09;
    amp[1] = 1.3567666284582589e - 06; shift[1] = 4.3373206800920752e - 08;
    amp[2] = 5.1757466437096689e - 06; shift[2] = 2.8970668616661478e - 07;
    amp[3] = 2.0060578172377753e - 05; shift[3] = 1.7970117665113235e - 06;
    amp[4] = 7.7976055655092961e - 05; shift[4] = 1.1016220840281374e - 05;
    amp[5] = 3.0323983201324125e - 04; shift[5] = 6.7403510935204510e - 05;
    amp[6] = 1.1793570136758038e - 03; shift[6] = 4.1228407663619111e - 04;
    amp[7] = 4.5868079395172696e - 03; shift[7] = 2.5216729791207432e - 03;
    amp[8] = 1.7839421514438226e - 02; shift[8] = 1.5423383071080004e - 02;
    amp[9] = 6.9386638849859295e - 02; shift[9] = 9.4337434071853923e - 02;
    amp[10] = 2.6997414529708952e - 01; shift[10] = 5.7713151658913675e - 01;
    amp[11] = 1.0526731536884490e + 00; shift[11] = 3.5350396633271388e + 00;
    amp[12] = 4.1584233028628317e + 00; shift[12] = 2.1815101171813343e + 01;
    amp[13] = 1.7800823020581991e + 01; shift[13] = 1.4102992696626504e + 02;
    amp[14] = 1.2795681699057995e + 02; shift[14] = 1.2544425313051306e + 03;
    */
  }

  mdp_real compute_gaussian_momenta(GaugeClass &U) {
    mdp_site x(U.lattice());
    mdp_real re, im;
    forallsites(x) {
      for(int mu = 0; mu<U.ndim; mu++) {
	U(x, mu) = 0;
	for(int a = 0; a<numgen; a++) {
	  re = U.lattice().random(x).gaussian();
	  // im = U.lattice().random(x).gaussian();
	  // U(x, mu) += ((re + I * im)/sqrt(2)) * lambda[a];
	  U(x, mu) += re * lambda[a];
	}
      }
    }
    U.update();
  }

  void set_gaussian(FermiClass &F) {
    mdp_site x(F.lattice());
    mdp_real re, im;
    forallsites(x) {
      for(int alpha = 0; alpha<F.nspin; alpha++) {
	for(int i = 0; i<F.nc; i++) {
	  re = F.lattice().random(x).gaussian();
	  im = F.lattice().random(x).gaussian();
	  F(x, alpha, i) = (re + I * im)/sqrt(2);
	}
      }
    }
    F.update();
  }

  mdp_real compute_kinetic_energy(GaugeClass &p_U, FermiClass &p_F) {
    mdp_complex tmp = 0;
    mdp_site x(p_U.lattice());
    forallsites(x)
      for(int mu = 0; mu<p_U.ndim; mu++)
	tmp -= 0.5 * trace(p_U(x,mu) * p_U(x,mu));
    mdp.add(tmp);
    tmp += p_F * p_F; // CHECK THE U.nc
    return tmp.real();
  }

  void compute_effective_links(GaugeClass &U, GaugeClass &V) {
    if(coeff["representation"] == FUNDAMENTAL && &U != &V)
      V = U;
    else if(coeff["representation"] == SYMMETRIC) {
      mdp_site x(U.lattice());
      forallsitesandcopies(x) {
	for(int mu = 0; mu<U.ndim; mu++) {
	  for(int a = 0; a<V.nc; a++) {
	    for(int b = 0; b<V.nc; b++) {
	      V(x, mu, a, b) = trace(S[a] * U(x, mu) * S[b] * transpose(U(x, mu)));
	    }
	  }
	}
      }
    } else if(coeff["representation"] == ANTISYMMETRIC) {
      throw string("ANTISYMMETRIC representation is not defined yet");
    } else if(coeff["representation"] == SO4) {
      // nothing to do
    }    
  }

  mdp_real compute_action(GaugeClass &U, GaugeClass &V, FermiClass &F) {
    double action_gauge = 0.0, action_fermi = 0.0;
    mdp_site x(U.lattice());
    FermiClass inverse_F(F);
    double c = coeff["beta"] *(U.ndim *(U.ndim - 1) * U.lattice().nvol_gl)/2;
    action_gauge = - c * average_plaquette(U);
    if (coeff["dynamical_fermions"]>0) {
      compute_effective_links(U, V);
      CG2::inverter(inverse_F, F, V, coeff,
		    coeff["cg_absolute_precision"],
		    coeff["cg_relative_precision"],
		    coeff["cg_max_steps"], true);
      action_fermi += real_scalar_product(F, inverse_F);
    }
    return action_gauge + action_fermi;
  }
  void compute_fields_evolution(GaugeClass &U,
				GaugeClass &p_U,
				GaugeClass &f_U,
				FermiClass &F,
				FermiClass &p_F,
				FermiClass &f_F) {
			
    int dynamical_fermions = coeff["dynamical_fermions"];
    GaugeClass nf_U(f_U);
    FermiClass nf_F(f_F);
    mdp_real dt = coeff["timestep"];
    mdp_site x(U.lattice());
    GaugeClass U_temp(U);
    FermiClass F_temp(p_F);
    
    // leap frog algorithm
    // F_temp = p_F;
    // F_temp.update();
    mdp_add_scaled_field(F_temp, 0.5 * dt, f_F);
    // F_temp.update();
    mdp_add_scaled_field(F, dt, F_temp);
    F.update();

    forallsites(x)
      for(int mu = 0; mu<U.ndim; mu++)
	U(x, mu) = exp(dt * p_U(x, mu) + 0.5 * dt * dt * f_U(x, mu)) * U(x, mu);
    U.update();

    compute_force(U, nf_U, F, nf_F);
    
    F_temp = 0;
    // F_temp.update();
    mdp_add_scaled_field(F_temp, 0.5 * dt, nf_F);
    // F_temp.update();
    // f_F.update();
    mdp_add_scaled_field(F_temp, 0.5 * dt, f_F);
    // F_temp.update();
    mdp_add_scaled_field(p_F, 1, F_temp);
    p_F.update();
    
    forallsites(x)
      for(int mu=0; mu<U_temp.ndim; mu++)
	U_temp(x,mu) = 0;
    // U_temp.update();
    mdp_add_scaled_field(U_temp, 0.5 * dt, nf_U);
    // U_temp.update();
    
    mdp_add_scaled_field(U_temp, 0.5 * dt, f_U);
    // U_temp.update();
    mdp_add_scaled_field(p_U, 1, U_temp);
    p_U.update();
    
    f_U = nf_U;
    f_U.update();
    f_F = nf_F;
    f_F.update();
  }
  
  void compute_force(GaugeClass &U,
		     GaugeClass &f_U,
		     FermiClass &F,
		     FermiClass &f_F) {
    
    int dynamical_fermions = coeff["dynamical_fermions"];
    mdp_site x(U.lattice());
    mdp_matrix staple(U.nc, U.nc);
    GaugeClass Udag(U.lattice(), U.nc), utmp(U.lattice(), U.nc);
    FermiClass sol(F.lattice(), F.nc, F.nspin);
    FermiClass psol(F.lattice(), F.nc, F.nspin);
    GaugeClass &V = *to_V;

    forallsitesandcopies(x)
      for(int mu=0; mu<U.ndim; mu++)
	Udag(x, mu) = hermitian(U(x, mu));

    forallsites(x) {
      for(int mu = 0; mu<U.ndim; mu++) {
	staple = 0;
	for(int nu = 0; nu<U.ndim; nu++)
	  if(nu != mu)
	    staple = staple + U(x + mu, nu) * Udag(x + nu, mu) * Udag(x, nu) +
	      Udag(x + mu - nu, nu) * Udag(x - nu, mu) * U(x - nu, nu);
	staple = U(x, mu) * staple - hermitian(U(x, mu) * staple);
	staple -= trace(staple) * mdp_identity(U.nc)/U.nc; //is this right?
	f_U(x, mu) = - coeff["beta"]/(2.0 * U.nc) * staple;
      }
    }
    f_U.update();

    if (coeff["dynamical_fermions"] > 0) {
      compute_effective_links(U, V);
      CG2::inverter(sol, F, V, coeff,
		    coeff["cg_absolute_precision"],
		    coeff["cg_relative_precision"],
		    coeff["cg_max_steps"], true);
      mul_Q(psol, sol, V, coeff);
      psol/= 2.0 * coeff["kappa"];
      psol.update();
      compute_fermion_forces(U, utmp, sol, psol);

      f_F = 0;
      f_F.update();
      mdp_add_scaled_field(f_F, - 1, sol);
      f_F.update();
      forallsites(x)
	for(int mu = 0; mu<U.ndim; mu++)
	  f_U(x, mu) -= utmp(x, mu);
      f_U.update();
    }
  }  
  void compute_fermion_forces(GaugeClass &U,
			      GaugeClass &f_U,
			      FermiClass &sol,
			      FermiClass &psol) {

    mdp_site x(U.lattice());
    int fnc = U.nc *(U.nc + 1)/2;
    mdp_matrix dum(U.nc, U.nc);
    mdp_matrix tmp1(U.nc, U.nc);
    mdp_matrix tmp2(U.nc, U.nc);
    mdp_matrix stemp;
    GaugeClass Utr(U);
    GaugeClass Udag(U);
    GaugeClass Udagtr(U);

    forallsites(x) {
      for(int mu = 0; mu<U.ndim; mu++) {
	Udag(x, mu) = hermitian(U(x, mu));
	Utr(x, mu) = transpose(U(x, mu));
	Udagtr(x, mu) = transpose(Udag(x, mu));
      }
    }
    Udag.update();
    Utr.update();
    Udagtr.update();

    if(coeff["representation"] == FUNDAMENTAL) {
      forallsites(x) {
	for(int mu = 0; mu<U.ndim; mu++) {
	  dum = 0;	  
	  for(int a = 0; a<U.nc; a++) {
	    for(int b = 0; b<U.nc; b++) {
	      stemp = 0.5 * hermitian(spinor(psol, x, a)) 
		* ((1 - Gamma[mu]) * spinor(sol, x + mu, b));
	      tmp1(b, a) = stemp(0, 0);
	      stemp = 0.5 * hermitian(spinor(psol, x + mu, a)) 
		* ((1 + Gamma[mu]) * spinor(sol, x, b));
	      tmp2(b, a) = stemp(0, 0);
	    }
	  }
	  dum = tmp2 * Udag(x, mu) - U(x, mu) * tmp1;
	  dum -= hermitian(dum);
	  dum -= trace(dum) * mdp_identity(U.nc)/U.nc;
	  f_U(x, mu) = dum;
	}
      }
    } else if(coeff["representation"] == SYMMETRIC) {
      forallsites(x) {
	for(int mu = 0; mu<U.ndim; mu++) {
	  dum = 0;
	  for(int a = 0; a<fnc; a++) {
	    for(int b = 0; b<fnc; b++) {
	      stemp = hermitian(spinor(psol, x, a)) *((1 - Gamma[mu]) * spinor(sol, x + mu, b));
	      dum -= stemp(0, 0) * U(x, mu) * S[b] * Utr(x, mu) * S[a];
	      stemp = hermitian(spinor(psol, x + mu, a)) *((1 + Gamma[mu]) * spinor(sol, x, b));
	      dum += stemp(0, 0) * S[b] * Udagtr(x, mu) * S[a] * Udag(x, mu);
	    }
	  }
	  dum -= hermitian(dum);
	  dum -= trace(dum) * mdp_identity(U.nc)/U.nc;
	  f_U(x, mu) = dum;
	}
      }
    } else if(coeff["representation"] == SO4) {
      mdp_real f;
      forallsites(x) {
	for(int mu = 0; mu<U.ndim; mu++) {
	  dum=0;         
	  for(int g=0;g<numgen;g++) {
	    f=2.0*trace(lambda[g]*tmp2*Udag(x,mu)).real()
	      -2.0*trace(lambda[g]*U(x,mu)*tmp1).real();
	    dum=dum-f*lambda[g];
	  }
	  f_U(x, mu) = dum;
	}
      }
    } else throw string("representation not supported");
    f_U.update();
  }

  static mdp_matrix spinor(FermiClass &psi, mdp_site x, int b) {
    mdp_matrix temp(psi.nspin,1);
    for(int i=0; i<psi.nspin; i++) temp(i,0)=psi(x,i,b);
    return temp;
  }
};

