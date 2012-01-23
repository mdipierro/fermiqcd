// based on code from James Hetrick
/* EXAMPLE OF USAGE
#include "fermiqcd.h"

int main(int argc, char** argv) {
  mdp.open_wormholes(argc,argv);
  define_base_matrices("FERMILAB");
  int L[] = {8,32,32,32};
  mdp_lattice lattice(4,L);
  gauge_field U(lattice,3);
  InstantonGenerator4D generator;
  for(int k=0; k<100; k++) {
    cout << k << endl;
    vector<SingleInstanton4D> instantons;
    instantons.push_back(SingleInstanton4D(7.5,15.5,15.5,10+0.1*k,1,+1)); // instanton
    instantons.push_back(SingleInstanton4D(7.5,15.5,16.5,20-0.1*k,1,-1)); // anti-instanton
    generator.generate(U,instantons);
    check_unitarity(U);
      mdp_site x(lattice);
      x.set(0,0,0,0);
      cout << U(x,0) << endl;
    topological_charge_vtk(U,string("engineered.topological.")+tostring(k)+".vtk",0);
  }
  mdp.close_wormholes();
  return 0;
}
*/

class SingleInstanton4D {
 public:
  float x[4];
  float rho;
  int charge; // +1 or -1;
  SingleInstanton4D(float x0, float x1, float x2, float x3, float rho, int charge) {
    x[0]=x0;
    x[1]=x1;
    x[2]=x2;
    x[3]=x3;
    this->rho=rho;
    this->charge=charge;
  }
};

class InstantonGenerator4D {
 private:
  mdp_lattice *lattice;
  mdp_matrix tau[4], bar[4];
  
  void init_tau_and_bar() {
    tau[0] = sigma[0]; // indentity
    tau[1] = I*sigma[1];
    tau[2] = I*sigma[2];
    tau[3] = I*sigma[3];
    bar[0] = sigma[0]; // indentity
    bar[1] = -I*sigma[1];
    bar[2] = -I*sigma[2];
    bar[3] = -I*sigma[3];
  }

  Matrix make_singular_instanton(float xl[4], int mu, SingleInstanton4D &instanton) {
    Matrix A;    
    float x2, x[4], rho2;
    for(int i=0; i<4; i++) x[i] = min((float)xl[i]-instanton.x[i],
				      (float)(*lattice).size(i)-xl[i]+instanton.x[i]);
    x2 = x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+x[3]*x[3];
    
    rho2 = pow(instanton.rho,2);
    if(instanton.charge>0) 
      A=(-0.5*rho2/(x2+rho2)/x2)*
	(bar[mu]*(x[0]*tau[0]+x[1]*tau[1]+x[2]*tau[2]+x[3]*tau[3])
	 -(x[0]*bar[0]+x[1]*bar[1]+x[2]*bar[2]+x[3]*bar[3])*tau[mu]);
    else
      A=(-0.5*rho2/(x2+rho2)/x2)*
	(tau[mu]*(x[0]*bar[0]+x[1]*bar[1]+x[2]*bar[2]+x[3]*bar[3])
	 -(x[0]*tau[0]+x[1]*tau[1]+x[2]*tau[2]+x[3]*tau[3])*bar[mu]);
    return A;
  }

  Matrix make_su2_link(mdp_site xn, int mu, vector<SingleInstanton4D> &instantons) {
    Matrix A(2,2);
    Matrix P = sigma[0];
    float x[4], a[3], norm_a;
    int steps = 5;
    float dx = 1.0/steps;
    float precision = 0.0001;
    for(int i=0; i<4; i++) x[i]=xn(i);
    for(int m=0; m<steps; m++) {
      x[mu] = (0.5+m)/steps + xn(mu);
      for(int k=0; k<instantons.size(); k++)
	A += make_singular_instanton(x,mu, instantons[k]);
      a[0] = imag(A(0,1))*dx;
      a[1] = real(A(0,1))*dx;
      a[2] = imag(A(0,0))*dx;
      norm_a = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
      if(norm_a > precision)
	P *= sin(norm_a)/norm_a*(a[0]*tau[1]+a[1]*tau[2]+a[2]*tau[3])+cos(norm_a)*tau[0];      
      /*
	else P*=1;
      */
    }
    return P;
  }

 public:
  void generate(gauge_field &U, vector<SingleInstanton4D> &instantons) {
    init_tau_and_bar();
    lattice = &U.lattice();
    mdp_site x(U.lattice());
    Matrix A;
    if(U.ndim!=4) throw string("instantons only in 4D");
    forallsites(x) {
      for(int mu=0; mu<U.ndim; mu++) {
	A = make_su2_link(x, mu, instantons);
	U(x,mu)=1;
	for(int i=0; i<2; i++)
	  for(int j=0; j<2; j++)
	    U(x,mu,i,j)=A(i,j);
      }      
    }
  }
};
