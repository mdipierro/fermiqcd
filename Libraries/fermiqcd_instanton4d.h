/////////////////////////////////////////////////////////////////
/// @file fermiqcd_fermi_algorithms.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// class for building a single instanton gauge configuration
///
/// Distributed under GPL2 License
/// 
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////      

class Instanton4D {
private:   
  static int epsilon123(int i, int j, int k) {
    if(i==j || j==k || i==k) return 0;
    if(i==1 && j==2 && k==3) return 1;
    if(i==1 && j==3 && k==2) return -1;
    if(i==2 && j==1 && k==3) return -1;
    if(i==2 && j==3 && k==1) return 1;
    if(i==3 && j==1 && k==2) return 1;
    if(i==3 && j==2 && k==1) return -1;
  }
public:
  vector<mdp_real> p; // location of the instanton
  int nc;
  int sub_i, sub_j;
  mdp_real charge; // this is 1/g
  mdp_real lambda;
  mdp_matrix eta[4][4];
  Instanton4D(int nc, int sub_i, int sub_j, mdp_real charge, mdp_real lambda, vector<mdp_real> &p) {
    this->nc=nc;
    this->sub_i=sub_i;
    this->sub_j=sub_j;
    this->lambda=lambda;
    this->charge=charge;
    this->p=p;
    mdp_matrix T(nc,nc);
    mdp_matrix sigma_rot1[4];    
    mdp_matrix sigma_rot2[4];
    mdp_matrix sigma_rot3[4];

    sigma_rot1[3]=sigma[3];
    sigma_rot1[1]=cos(alpha)*sigma[1]+sin(alpha)*sigma[2];
    sigma_rot1[2]=-sin(alpha)*sigma[1]+cos(alpha)*sigma[2];

    sigma_rot2[1]=sigma_rot1[1];
    sigma_rot2[2]=cos(beta)*sigma_rot1[2]+sin(beta)*sigma_rot1[3];
    sigma_rot2[3]=-sin(beta)*sigma_rot1[2]+cos(beta)*sigma_rot1[3];

    sigma_rot3[3]=sigma_rot2[3];
    sigma_rot3[1]=cos(gamma)*sigma_rot2[1]+sin(gamma)*sigma_rot2[2];
    sigma_rot3[2]=-sin(gamma)*sigma_rot2[1]+cos(gamma)*sigma_rot2[2];
      
    for(int mu=0; mu<4; mu++)
      for(int nu=0; nu<4; nu++) {
	this->eta[mu][nu].dimension(nc,nc);      
	this->eta[mu][nu]=0;
	for(int a=1; a<4; a++) {
	  T=0;
	  T(sub_i,sub_i)=sigma_rot3[a](0,0);
	  T(sub_i,sub_j)=sigma_rot3[a](0,1);
	  T(sub_j,sub_i)=sigma_rot3[a](1,0);
	  T(sub_j,sub_j)=sigma_rot3[a](1,1);
	  if(a>0 && mu>0 && nu>0) {
	    this->eta[mu][nu]+=epsilon123(a,mu,nu)*T;
	  } else if (mu==0 && a>0 && a==nu) {
	    this->eta[mu][nu]-=T; // T[a]*(-1)*delta(a,nu)
	  } else if (nu==0 && a>0 && a==mu) {
	    this->eta[mu][nu]+=T; // T[a]*delta(a,mu)
	  }
	}
      }
    /*
    for(int mu=0; mu<4; mu++)
      for(int nu=0; nu<4; nu++) {
	this->eta[mu][nu].dimension(nc,nc);      
	for(int i=0; i<nc; i++)
	  for(int j=0; j<nc; j++)
	    if(i==j && (i!=sub_i && i!=sub_j)) {
	      this->eta[mu][nu](i,j) = 1;
	    } else if((i!=sub_i && i!=sub_j) || (j!=sub_i && j!=sub_j)) {
	      this->eta[mu][nu](i,j) = 0;
	    } else if(mu==0 && nu==0) {
	      this->eta[mu][nu](i,j) = 0;
	    } else {
	      int i0=(i==sub_i?0:1);
	      int j0=(j==sub_i?0:1);
	      for(int a=1; a<4; a++) {
		if(mu==0) {
		  //this->eta[mu][nu](i,j) += ((a==nu)?-1:0);
		  this->eta[mu][nu](i,j) += sigma[a](i0,j0)*((a==nu)?-1:0);
		} else if(nu==0) {
		  //this->eta[mu][nu](i,j) += ((a==mu)?+1:0);
		  this->eta[mu][nu](i,j) += sigma[a](i0,j0)*((a==mu)?+1:0);
		} else {
		  // this->eta[mu][nu](i,j) += epsilon123(a,mu,nu);
		  this->eta[mu][nu](i,j) += sigma[a](i0,j0)*epsilon123(a,mu,nu);
		}	
	      }
	    }
      }
    */
  }
  mdp_matrix operator()(mdp_site &x, int mu) {
    int v[4];
    mdp_lattice &lattice=x.lattice();
    for(int nu=0; nu<4; nu++) v[nu] = x(nu)-this->p[nu];
    float d2 = v[0]*v[0]+v[1]*v[1]+v[2]*v[2]+v[3]*v[3];
    mdp_matrix A(this->nc,this->nc);
    A=0;
    for(int nu=0; nu<4; nu++)
      A += this->eta[mu][nu]*v[nu];
    return (2.0*this->charge/(d2+this->lambda*this->lambda))*A;
  }
};

/* usage example

int main(int argc, char** argv) {
  define_base_matrices("FERMILAB");
  mdp.open_wormholes(argc,argv);    
  int nc=3;
  int L[]={10,40,10,10};
  mdp_lattice lattice(4,L,default_partitioning0,torus_topology,0,2,false);
  gauge_field U(lattice,nc);
  char filename[128];
  mdp_site x(lattice);
  set_cold(U);
  vector<mdp_real> p(4);
  for(int mu=0; mu<4; mu++) p[mu]=L[mu]/2;

  for(int k=0; k<40; k++) {
    p[1]=10.0+0.25*k;
    Instanton4D A1(3,0,1, 1.0,3,p);
    p[1]=30.0-0.25*k;
    Instanton4D A2(3,0,1, 1.0,3,p);
    
    forallsites(x)
      for(int mu=0; mu<U.ndim; mu++)
	U(x,mu)=exp(-I*(A1(x,mu)+A2(x,mu)));
    
    mdp << "k=" << k << endl;
    //if(k>0) ApeSmearing::smear(U,0.7,1,10);
    sprintf(filename,"%s.%i.vtk",argv[0],k);
    mdp << "top=" << topological_charge_vtk(U,filename) << endl;    
  }
  mdp.close_wormholes();
  return 0;
}

*/
