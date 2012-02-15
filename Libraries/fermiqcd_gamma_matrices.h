/////////////////////////////////////////////////////////////////
/// @file fermiqcd_gamma_matrices.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Declares:
/// - GammaI = \f$ 1 \f$
/// - Gamma[mu] = \f$ \gamma^\mu \f$
/// - Gamma5 = \f$ \gamma^5 \f$
/// - Sigma[mu][nu] = \f$ \sigma^{\mu\nu} \f$
/// - sigma[i] = \f$ \sigma^{i} \f$ (generators SU(2))
/// - Lambda[i] = \f$ \lambda^{i} \f$ (generators SU(3))
///
/// Distributed under GPL2 License
/// 
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////


// modified 20 Nov 2009 to include MILC convention

#define GAMMA_MATRICES

// ////////////////////////////////////////////////////////////////////
// Gamma matrices and relatives!
// ////////////////////////////////////////////////////////////////////

mdp_complex Gamma_val[4][4], Sigma_val[4][4][4];
int     Gamma_idx[4][4], Sigma_idx[4][4][4];
mdp_complex Gamma5_val[4], GammaxGamma5_val[4][4];
int     Gamma5_idx[4], GammaxGamma5_idx[4][4];
mdp_complex Gamma_valr[4][4], Sigma_valr[4][4][4];
int     Gamma_idxr[4][4], Sigma_idxr[4][4][4];
mdp_complex Gamma5_valr[4], GammaxGamma5_valr[4][4];
int     Gamma5_idxr[4], GammaxGamma5_idxr[4][4];
int     G16_idx[16][4];
mdp_complex G16_val[16][4];

mdp_matrix Gamma[4], Gamma1, Gamma5, Pleft, Pright, Lambda[9];
mdp_matrix Sigma[4][4], sigma[4];

/// @brief define convetion for Gamma matrices and bases of SU(2) and SU(3)
///
/// At the beginning of any FermiQCD program you MUST call
/// @verbatim
/// define_base_matrices("FERMILAB");
/// @endverbatim
/// Possible convections are:
/// - "FERMILAB" (ok for lattice qcd)
/// - "MILC" (ok for lattice qcd)
/// - "UKQCD" (ok for lattice qcd)
/// - "Minkowsky-Dirac" (not ok for lattice qcd)
/// - "Minkowsky-Chiral" (not ok for lattice qcd)
/// Convention can be changed within the program.
void define_base_matrices(string convention="FERMILAB") {
  begin_function("define_base_matrices");
  int i,j,a,b,mu,nu; // was int i,j,k;
  for(i=0; i<4; i++) Gamma[i].dimension(4,4);
  Gamma5.dimension(4,4);

  if(convention=="UKQCD") {
    Gamma[0](0,0)=+1;
    Gamma[0](1,1)=+1;
    Gamma[0](2,2)=-1;
    Gamma[0](3,3)=-1;
    
    Gamma[1](0,3)=+I;
    Gamma[1](1,2)=+I;
    Gamma[1](2,1)=-I;
    Gamma[1](3,0)=-I;
    
    Gamma[2](0,3)=+1;
    Gamma[2](1,2)=-1;
    Gamma[2](2,1)=-1;
    Gamma[2](3,0)=+1;
    
    Gamma[3](0,2)=+I;
    Gamma[3](1,3)=-I;
    Gamma[3](2,0)=-I;
    Gamma[3](3,1)=+I;
    Gamma5=Gamma[1]*Gamma[2]*Gamma[3]*Gamma[0];
  }

  if(convention=="FERMILAB") {
    Gamma[0](0,0)=+1;
    Gamma[0](1,1)=+1;
    Gamma[0](2,2)=-1;
    Gamma[0](3,3)=-1;

    Gamma[1](0,3)=+1;
    Gamma[1](1,2)=+1;
    Gamma[1](2,1)=+1;
    Gamma[1](3,0)=+1;
    
    Gamma[2](0,3)=-I;
    Gamma[2](1,2)=+I;
    Gamma[2](2,1)=-I;
    Gamma[2](3,0)=+I;
    
    Gamma[3](0,2)=+1;
    Gamma[3](1,3)=-1;
    Gamma[3](2,0)=+1;
    Gamma[3](3,1)=-1;

    Gamma5=Gamma[0]*Gamma[1]*Gamma[2]*Gamma[3];
  }

    if(convention=="MILC") {
    Gamma[0](0,2)=+1;
    Gamma[0](1,3)=+1;
    Gamma[0](2,0)=+1;
    Gamma[0](3,1)=+1;

    Gamma[1](0,3)=+I;
    Gamma[1](1,2)=+I;
    Gamma[1](2,1)=-I;
    Gamma[1](3,0)=-I;
    
    Gamma[2](0,3)=-1;
    Gamma[2](1,2)=+1;
    Gamma[2](2,1)=+1;
    Gamma[2](3,0)=-1;
    
    Gamma[3](0,2)=+I;
    Gamma[3](1,3)=-I;
    Gamma[3](2,0)=-I;
    Gamma[3](3,1)=+I;

    Gamma5=Gamma[1]*Gamma[2]*Gamma[3]*Gamma[0];

  }

  if(convention=="Minkowsky-Dirac") {
    Gamma[0](0,0)=+1;
    Gamma[0](1,1)=+1;
    Gamma[0](2,2)=-1;
    Gamma[0](3,3)=-1;
    
    Gamma[1](0,3)=+1;
    Gamma[1](1,2)=+1;
    Gamma[1](2,1)=-1;
    Gamma[1](3,0)=-1;
    
    Gamma[2](0,3)=-I;
    Gamma[2](1,2)=+I;
    Gamma[2](2,1)=+I;
    Gamma[2](3,0)=-I;
    
    Gamma[3](0,2)=+1;
    Gamma[3](1,3)=-1;
    Gamma[3](2,0)=-1;
    Gamma[3](3,1)=+1;
    Gamma5=I*Gamma[0]*Gamma[1]*Gamma[2]*Gamma[3];
  }

  if(convention=="Minkowsky-Chiral") {
    Gamma[0](0,3)=-1;
    Gamma[0](1,2)=-1;
    Gamma[0](2,1)=-1;
    Gamma[0](3,0)=-1;
    
    Gamma[1](0,3)=+1;
    Gamma[1](1,2)=+1;
    Gamma[1](2,1)=-1;
    Gamma[1](3,0)=-1;
    
    Gamma[2](0,3)=-I;
    Gamma[2](1,2)=+I;
    Gamma[2](2,1)=+I;
    Gamma[2](3,0)=-I;
    
    Gamma[3](0,2)=+1;
    Gamma[3](1,3)=-1;
    Gamma[3](2,0)=-1;
    Gamma[3](3,1)=+1;
    Gamma5=I*Gamma[0]*Gamma[1]*Gamma[2]*Gamma[3];
  }

  if(convention=="CHIRAL") {
    Gamma[0](0,3)=-1;
    Gamma[0](1,2)=-1;
    Gamma[0](2,1)=-1;
    Gamma[0](3,0)=-1;
    
    Gamma[1](0,3)=+I;
    Gamma[1](1,2)=+I;
    Gamma[1](2,1)=-I;
    Gamma[1](3,0)=-I;
    
    Gamma[2](0,3)=+1;
    Gamma[2](1,2)=-1;
    Gamma[2](2,1)=-1;
    Gamma[2](3,0)=+1;
    
    Gamma[3](0,2)=+I;
    Gamma[3](1,3)=-I;
    Gamma[3](2,0)=-I;
    Gamma[3](3,1)=+I;
    Gamma5=Gamma[0]*Gamma[1]*Gamma[2]*Gamma[3];
  }
  Gamma1=mdp_identity(4);
  Pleft  = mdp_complex(.5,0)*(Gamma[0]*Gamma[0]-Gamma5);
  Pright = mdp_complex(.5,0)*(Gamma[0]*Gamma[0]+Gamma5);

  /* sigma matrices SU(2) generators */

  for(i=0; i<4; i++)
    sigma[i].dimension(2,2);
  sigma[0](0,0)=+1;
  sigma[0](1,1)=+1;
  sigma[1](0,1)=+1;
  sigma[1](1,0)=+1;
  sigma[2](0,1)=-I;
  sigma[2](1,0)=+I;
  sigma[3](0,0)=+1;
  sigma[3](1,1)=-1;

  /* Lambda matrices SU(3) generators */

  for(i=0; i<9; i++) Lambda[i].dimension(3,3);

  Lambda[0](0,0)=+1;
  Lambda[0](1,1)=+1;
  Lambda[0](2,2)=+1;

  Lambda[1](0,1)=+1;
  Lambda[1](1,0)=+1;

  Lambda[2](0,1)=-I;
  Lambda[2](1,0)=+I;

  Lambda[3](0,0)=+1;
  Lambda[3](1,1)=-1;
  
  Lambda[4](0,2)=+1;
  Lambda[4](2,0)=+1;

  Lambda[5](0,2)=-I;
  Lambda[5](2,0)=+I;

  Lambda[6](1,2)=+1;
  Lambda[6](2,1)=+1;

  Lambda[7](1,2)=-I;
  Lambda[7](2,1)=+I;

  Lambda[8](0,0)=pow(3.0,-0.5);
  Lambda[8](1,1)=pow(3.0,-0.5);
  Lambda[8](2,2)=-2.0*pow(3.0,-0.5);

  /* Sigma[mu][nu] */

  mdp_matrix tmp;
  
  for(i=0; i<4; i++)
    for(j=0; j<4; j++) {
      Sigma[i][j].dimension(4,4);
      if(i!=j)
	Sigma[i][j]=I*Gamma[i]*Gamma[j];
      else
	Sigma[i][j]=0;
    }

  /* Filling arrays for fast Gamma multiplication */

  for(a=0; a<4; a++) 
    for(b=0; b<4; b++) {
      if(Gamma5(a,b)!=mdp_complex(0,0)) {
	Gamma5_val[b]=Gamma5(a,b);
	Gamma5_idx[b]=a;
	Gamma5_valr[a]=Gamma5(a,b);
	Gamma5_idxr[a]=b;	
      }
    }

  for(mu=0; mu<4; mu++) {
    tmp=Gamma[mu]*Gamma5;
    for(a=0; a<4; a++) 
      for(b=0; b<4; b++) {
	if(Gamma[mu](a,b)!=mdp_complex(0,0)) {
	  Gamma_val[mu][b]=Gamma[mu](a,b);
	  Gamma_idx[mu][b]=a;
	  Gamma_valr[mu][a]=Gamma[mu](a,b);
	  Gamma_idxr[mu][a]=b;
	}
	if(tmp(a,b)!=mdp_complex(0,0)) {
	  GammaxGamma5_val[mu][b]=tmp(a,b);
	  GammaxGamma5_idx[mu][b]=a;
	  GammaxGamma5_valr[mu][a]=tmp(a,b);
	  GammaxGamma5_idxr[mu][a]=b;
	}
      }
  }

  /* Filling arrays for fast Sigma multiplication */
  
  for(mu=0; mu<4; mu++) 
    for(nu=0; nu<4; nu++) 
      for(a=0; a<4; a++) 
	for(b=0; b<4; b++) {
	  if(Sigma[mu][nu](a,b)!=mdp_complex(0,0)) {
	    Sigma_val[mu][nu][b]=Sigma[mu][nu](a,b);
	    Sigma_idx[mu][nu][b]=a;
	    Sigma_valr[mu][nu][a]=Sigma[mu][nu](a,b);
	    Sigma_idxr[mu][nu][a]=b;
	  }
	}

  // Filling G16;
  
  for(a=0; a<4; a++) {
    G16_idx[0][a]=a;
    G16_val[0][a]=1;
    G16_idx[1][a]=Gamma5_idx[a];
    G16_val[1][a]=Gamma5_val[a];
    for(mu=0; mu<4; mu++) {
      G16_idx[2+mu][a]=Gamma_idx[mu][a];
      G16_val[2+mu][a]=Gamma_val[mu][a];
      G16_idx[6+mu][a]=Gamma_idx[mu][Gamma5_idx[a]];
      G16_val[6+mu][a]=Gamma_val[mu][Gamma5_idx[a]]*Gamma5_val[a];
    }
    for(mu=1; mu<4; mu++) {
      G16_idx[9+mu][a]=Gamma_idx[0][Gamma_idx[mu][a]];
      G16_val[9+mu][a]=Gamma_val[0][Gamma_idx[mu][a]]*Gamma_val[mu][a];
      G16_idx[12+mu][a]=Gamma_idx[0][Gamma_idx[mu][Gamma5_idx[a]]];
      G16_val[12+mu][a]=Gamma_val[0][Gamma_idx[mu][Gamma5_idx[a]]]*
	Gamma_val[mu][Gamma5_idx[a]]*Gamma5_val[a];
    }
  }
  end_function("define_base_matrices");
}

mdp_matrix parse_gamma(string g) {
  if(g=="I") return Gamma1;
  if(g=="0") return Gamma[0];
  if(g=="1") return Gamma[1];
  if(g=="2") return Gamma[2];
  if(g=="3") return Gamma[3];
  if(g=="5") return Gamma5;
  if(g=="05") return Gamma[0]*Gamma5;
  if(g=="15") return Gamma[1]*Gamma5;
  if(g=="25") return Gamma[2]*Gamma5;
  if(g=="35") return Gamma[3]*Gamma5;
  if(g=="01") return Gamma[0]*Gamma[1];
  if(g=="02") return Gamma[0]*Gamma[2];
  if(g=="03") return Gamma[0]*Gamma[3];
  if(g=="12") return Gamma[1]*Gamma[2];
  if(g=="13") return Gamma[1]*Gamma[3];
  if(g=="23") return Gamma[2]*Gamma[3];
  throw string("undefined gamma structure");
}
