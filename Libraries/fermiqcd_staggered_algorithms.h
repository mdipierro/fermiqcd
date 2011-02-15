/////////////////////////////////////////////////////////////////
/// @file fermiqcd_staggered_algorithms.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Various stuff for staggered fermions
///
/// Distributed under GPL2 License
/// 
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////

mdp_matrix Omega4x4(mdp_site x) {
  mdp_matrix M(4,4);
  M=1;
  if(x(3) % 2==1) M=Gamma[3]*M;
  if(x(2) % 2==1) M=Gamma[2]*M;
  if(x(1) % 2==1) M=Gamma[1]*M;
  if(x(0) % 2==1) M=Gamma[0]*M;
  return M;
}

/// Pointer to current Staggered/Asqtad action
void (*default_staggered_action)(staggered_field &,
				 staggered_field &,
				 gauge_field &, 
				 coefficients &, int) = StaggeredAsqtadActionFast::mul_Q;


/// Executes current Staggered/Asqtad action
void mul_Q(staggered_field &psi_out,
	   staggered_field &psi_in,
	   gauge_field &U, 
	   coefficients &coeff,
	   int parity=EVENODD) {
  (*default_staggered_action)(psi_out, psi_in, U, coeff, parity);
}

// ////////////////////////////////////////////////
// choice of the default inversion method
// ////////////////////////////////////////////////

/// Pointer to current Staggered/Asqtad inverter
inversion_stats (*default_staggered_inverter)(staggered_field &, 
					      staggered_field &,
					      gauge_field &, 
					      coefficients &,
					      mdp_real, mdp_real,int) 
  =&(BiCGStab::inverter<staggered_field,gauge_field>);

/// Executes current Staggered/Asqtad inverter
inversion_stats mul_invQ(staggered_field &psi_out, 
			 staggered_field &psi_in, 
			 gauge_field &U, 
			 coefficients &coeff, 
			 mdp_real absolute_precision=staggered_inversion_precision,
			 mdp_real relative_precision=0, int max_steps=2000) {
  return (*default_staggered_inverter)(psi_out, psi_in, U, coeff, absolute_precision, relative_precision, max_steps);
}

// ////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////
// Lepage improved fat links for improved staggered action
// ////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////

/// Takes a plaquette and a type of action and returns a 1D array
/// with weights of paths required to build fat links for the action
/// @see lepage_improved_links()
mdp_array<mdp_real,1> lepage_coefficients(mdp_real plaquette, char type[]) {
  begin_function("lepage_coefficients");

  mdp_real u0=pow((double)plaquette,(double)0.25);
  mdp_array<mdp_real,1> c(6); 

  c[0]=c[1]=c[2]=c[3]=c[4]=c[5]=0.0;
  if(strcmp(type,"None")==0) {
    c[0]=1;
  }
  if(strcmp(type,"Staple+Naik")==0) {
    c[0]=9.0/(8*4);
    c[1]=9.0/(8*8)*pow(u0,-2);
    c[2]=c[3]=c[4]=0;
    c[5]=-9.0/(8*27);
  }
  if(strcmp(type,"Fat3")==0) {
    error("Not implemented");
  }
  if(strcmp(type,"Fat5")==0) {
    c[0]=1.0/7;
    c[1]=1.0/(7*2)*pow(u0,-2);
    c[2]=1.0/(7*8)*pow(u0,-4);
    c[3]=c[4]=0;
    c[5]=0;
  }
  if(strcmp(type,"Fat7")==0) {
    c[0]=1.0/8;
    c[1]=1.0/(8*2)*pow(u0,-2);
    c[2]=1.0/(8*8)*pow(u0,-4);
    c[3]=1.0/(8*48)*pow(u0,-6);
    c[4]=0;
    c[5]=0;
  }
  if(strcmp(type,"Full")==0) {
    c[0]=5.0/8.0;
    c[1]=1.0/16.0*pow(u0,-2); 
    c[2]=1.0/64.0*pow(u0,-4);
    c[3]=1.0/(8.0*48.0)*pow(u0,-6);
    c[4]=-1.0/16.0*pow(u0,-4);
    c[5]=-1.0/24.0*pow(u0,-2);
  }
  
  mdp << "\tImprovement type=" << type << '\n';
  mdp << "Coefficients:\n";
  mdp << "u_0=" << u0 << '\n';
  mdp << "link=" << c[0] << '\n';
  mdp << "3staple=" << c[1] << '\n';
  mdp << "5staple=" << c[2] << '\n';
  mdp << "7staple=" << c[3] << '\n';
  mdp << "lepage=" << c[4] << '\n';
  mdp << "naik=" << c[5] << '\n';

  end_function("lepage_coefficients");

  return c;
}



// ///////////////////////////////////////////////
// Note we cannot use Nmdp_matrixField here
// because it does not implement the twist
// we use vectors of gauge_fields instead!
// ///////////////////////////////////////////////

/// Takes a gauge field U and a set of coefficients as computed by
/// lepage_coefficients() and fills the gauge field V with fat links and
/// Long links
///
/// Example:
/// @verbatim
/// gauge_field U(lattice,nc);
/// gauge_field V(lattice,nc);
/// U.load("myfield");
/// float p=1.0; // the average plaquette
/// lepage_improved_links(V,U,lepage_coefficients(p,"Full"),false);
/// /// now use V instead of U for staggered actions and inverters
/// @endverbatim
/// Note that the type of action can be
/// - "Full" for full as asqtad
/// - "Staple+Naik"
/// - "Fat3"
/// - "Fat5"
/// - "Fat7"
/// Also note that if project==true the fat links are projected back to SU(nc)
void lepage_improved_links(gauge_field &V, 
			   gauge_field &U, 
			   mdp_array<mdp_real,1> c, 
			   int project=false) {

  begin_function("lepage_improved_links");

  if(U.ndim!=4) error("fermiqcd_staggered_auxiliary_functions/lepage_improved_links: ndim!=4 (contact <mdp@fnal.gov>)");

  // this is a list of permitations of 0,1,2,3

  
  const int epsilon[24][5] = {{0, 1, 2, 3, 0}, {0, 1, 3, 2, 3},
			      {0, 2, 1, 3, 1}, {0, 2, 3, 1, 4},
			      {0, 3, 1, 2, 2}, {0, 3, 2, 1, 5},
			      {1, 0, 2, 3, 0}, {1, 0, 3, 2, 3},
			      {1, 2, 0, 3, 1}, {1, 2, 3, 0, 4},
			      {1, 3, 0, 2, 2}, {1, 3, 2, 0, 5},
			      {2, 0, 1, 3, 0}, {2, 0, 3, 1, 3},
			      {2, 1, 0, 3, 1}, {2, 1, 3, 0, 4},
			      {2, 3, 0, 1, 2}, {2, 3, 1, 0, 5},
			      {3, 0, 1, 2, 0}, {3, 0, 2, 1, 3},
			      {3, 1, 0, 2, 1}, {3, 1, 2, 0, 4},
			      {3, 2, 0, 1, 2}, {3, 2, 1, 0, 5}}; 
  // mu, nu, rho, sig, im2 

  
  int nc=U.nc;
  int ndim=U.ndim;
  int mu,nu,rho, idx, imn, im2, i,j;
  site x(U.lattice());
  site y(U.lattice());
  mdp_matrix b1(nc,nc), b2(nc,nc);

  mdp << "Allocating temporary vectors...";

  gauge_field *Delta1=new gauge_field[(ndim-1)];
  gauge_field *Delta2=new gauge_field[(ndim-1)*(ndim-2)];

  mdp << "done.\n";

  for(imn=0; imn<(ndim-1); imn++)
    Delta1[imn].allocate_gauge_field(U.lattice(),nc);
  for(im2=0; im2<(ndim-1)*(ndim-2); im2++)
    Delta2[im2].allocate_gauge_field(U.lattice(),nc);

  // ////////////////////////////////////////////////////////////////////
  // These coefficints are in agreement with Irginos, Toussaint and Sugar
  // hep-lat/9903032
  // equivalent but easier than the Lepage version
  // ////////////////////////////////////////////////////////////////////

  
  mdp << "Computing Fat Links for Lepage improved action:\n";

  // LINK
  mdp << "link...\n";

  forallsites(x)
    for(mu=0; mu<ndim; mu++)
      V(x,mu)=c[0]*U(x,mu);

  // 3 STAPLE
  mdp << "3staple...\n";

  forallsites(x)
    for(mu=0; mu<ndim; mu++)
      for(nu=0; nu<ndim; nu++) if(nu!=mu) {
	imn=(nu<mu)?nu:(nu-1);
	Delta1[imn](x,mu) =U(x,nu)*U(x+nu,mu)*hermitian(U(x+mu,nu));
	y=x-nu;
	Delta1[imn](x,mu)+=hermitian(U(y,nu))*U(y,mu)*U(y+mu,nu);
	V(x,mu)+=c[1]*Delta1[imn](x,mu);
      }

  
  for(imn=0; imn<(ndim-1); imn++)
    Delta1[imn].update();  

  // 5 STAPLE
  mdp << "5staple...\n";

  forallsites(x) 
    for(idx=0; idx<24; idx++) {
      mu =epsilon[idx][0];
      nu =epsilon[idx][1];
      rho=epsilon[idx][2];
      im2=epsilon[idx][4];
      imn=(nu<mu)?nu:(nu-1);
      Delta2[im2](x,mu) =U(x,rho)*Delta1[imn](x+rho,mu)*hermitian(U(x+mu,rho));
      y=x-rho;
      Delta2[im2](x,mu)+=hermitian(U(y,rho))*Delta1[imn](y,mu)*U(y+mu,rho);
      V(x,mu)+=c[2]*Delta2[im2](x,mu);
    }
  for(im2=0; im2<(ndim-1)*(ndim-2); im2++)
    Delta2[im2].update();
  
  // 7 STAPLE
  mdp << "7staple...\n";

  if(c[3]!=0) 
    forallsites(x) 
      for(idx=0; idx<24; idx++) {
	mu =epsilon[idx][0];
	rho=epsilon[idx][3];
	im2=epsilon[idx][4];
	b2=U(x,rho)*Delta2[im2](x+rho,mu)*hermitian(U(x+mu,rho));
	y=x-rho;
	b2+=hermitian(U(y,rho))*Delta2[im2](y,mu)*U(y+mu,rho);
	V(x,mu)+=c[3]*b2;
      }
  
  // LEPAGE TERM UP
  mdp << "lepage...\n";

  if(c[4]!=0) {
    forallsites(x)
      for(mu=0; mu<ndim; mu++)
	for(nu=0; nu<ndim; nu++) if(nu!=mu) {
	  imn=(nu<mu)?nu:(nu-1);
	  Delta1[imn](x,mu)=U(x,nu)*U(x+nu,mu)*hermitian(U(x+mu,nu));
	}
    for(imn=0; imn<(ndim-1); imn++)
      Delta1[imn].update();

    forallsites(x)
      for(mu=0; mu<ndim; mu++)
	for(nu=0; nu<ndim; nu++) if(nu!=mu) {
	  imn=(nu<mu)?nu:(nu-1);
	  b2=U(x,nu)*Delta1[imn](x+nu,mu)*hermitian(U(x+mu,nu));
	  V(x,mu)+=c[4]*b2;
	}
    
    // LEPAGE TERM DOWN
    forallsites(x)
      for(mu=0; mu<ndim; mu++)
	for(nu=0; nu<ndim; nu++) if(nu!=mu) {
	  imn=(nu<mu)?nu:(nu-1);
	  y=x-nu;
	  Delta1[imn](x,mu)=hermitian(U(y,nu))*U(y,mu)*U(y+mu,nu);
	}
    for(imn=0; imn<(ndim-1); imn++)
      Delta1[imn].update();
    
    forallsites(x)
      for(mu=0; mu<ndim; mu++)
	for(nu=0; nu<ndim; nu++) if(nu!=mu) {
	  imn=(nu<mu)?nu:(nu-1);
	  y=x-nu;
	  b2=hermitian(U(y,nu))*Delta1[imn](y,mu)*U(y+mu,nu);
	  V(x,mu)+=c[4]*b2;
	}
  }
  if(project==true) {
    mdp << "projecting to SU(n)...\n";

    forallsites(x)
      for(mu=0; mu<4; mu++)
	V(x,mu)=project_SU(V(x,mu));
  }
  V.update();

  if(c[5]!=0) {
    mdp << "naik...\n";
    compute_long_links(V,U,3);
    mdp << "normalizing naik...\n";

    forallsitesandcopies(x)
      for(mu=0; mu<U.ndim; mu++)
	for(i=0; i<U.nc; i++)
	  for(j=0; j<U.nc; j++)
	    V.long_links(x,mu,i,j)*=c[5];
  }

  cout << "Freeing temporary vectors...";

  for(imn=0; imn<(ndim-1); imn++)
    Delta1[imn].deallocate_memory();
  for(im2=0; im2<(ndim-1)*(ndim-2); im2++)
    Delta2[im2].deallocate_memory();

  cout << "done.\n";

  end_function("lepage_improved_links");
}

void staggered_rephase(gauge_field &U, staggered_field &chi) {
  
  begin_function("staggered_rephase");
  
  site x(U.lattice());
  int  mu,i,j;
  forallsites(x) {
    for(mu=0; mu<U.ndim; mu++) {
      for(i=0; i<U.nc; i++)
	for(j=0; j<U.nc; j++) 
	  U(x,mu,i,j)=U(x,mu,i,j)*chi.eta(x,mu);
    }
    // this takes car of antiperiodic bc
    if(x(0)==U.lattice().size(0)-1) {
      for(i=0; i<U.nc; i++)
	for(j=0; j<U.nc; j++) 
	  U(x,0,i,j)*=-1;
    }
  }

  end_function("staggered_rephase");

}

