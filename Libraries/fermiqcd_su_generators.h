/////////////////////////////////////////////////////////////////
/// @file fermiqcd_su_generators.h
/// @version 11-3-2009
/// @author Simon Catterall and Massimo Di Pierro
///                                                                          
//////////////////////////////////////////////////////////////////
// Example:
// #include "fermiqcd.h"
// int main() {
//  SU_Generators g(3);
//  for(int a=0; a<g.ngenerators; a++)
//    cout << "g=" << g.lambda[a] << endl;
//  return 0;
// }

class SU_Generators {
public:
  vector<mdp_matrix> lambda;
  int n;
  int ngenerators;
  SU_Generators(int n) {
    this->n = n;
    this->ngenerators = (n*n-1);
    lambda.resize(ngenerators);
    
    for(int a=0; a<ngenerators; a++)
      lambda[a] = build_matrix(a) ;      
  }

  mdp_matrix build_matrix(int a) {
    
    int pos1 = (n*(n-1)/2);
    int pos2 = (n*(n-1));    
    mdp_matrix mat(n,n); 
    mdp_complex mult=0;    
    vector<mdp_complex> vec(ngenerators);
    for(int i=0; i<ngenerators; i++)
      vec[i] = (i==a)?1:0;
    
    for(int i=0, a=0; i<n; i++) {
      mat(i,i) = 0;
      for(int j=i+1; j<n; j++) {
	mat(i,j) = 0.5*(vec[a] - I*vec[pos1+a]);
	mat(j,i) = 0.5*(vec[a] + I*vec[pos1+a]);
	a++;
      }
    }
    for(int i=0; i<n-1; i++) {      
      mult = vec[pos2+i] * (1.0/sqrt(2.+2./(1.+i))/(1.+i));
      for(int j=0; j<i+1; j++) {
	mat(j,j) += mult;
      }      
      mat(i+1,i+1) -= (1+i)*mult; 
    }    
    return I * sqrt(2) * mat;
  }
};

