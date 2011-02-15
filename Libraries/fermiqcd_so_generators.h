/////////////////////////////////////////////////////////////////
/// @file fermiqcd_su_generators.h
/// @version 11-3-2009
/// @author Simon Catterall and Massimo Di Pierro
///                                                                          
//////////////////////////////////////////////////////////////////
// Example:
// #include "fermiqcd.h"
// int main() {
//  SO_Generators g(4);
//  for(int a=0; a<g.ngenerators; a++)
//    cout << "g=" << g.lambda[a] << endl;
//  return 0;
// }

class SO_Generators {
public:
  vector<mdp_matrix> lambda;
  int n;
  int ngenerators;
  SO_Generators(int n) {
    this->n = n;
    this->ngenerators = n*(n-1)/2;
    lambda.resize(ngenerators);
    mdp_matrix temp(n,n);
    mdp_complex z = 1.0/sqrt(2.0);
    int k=0;
    for(int j=0; j<n-1; j++)
      for(int i=j+1; i<n; i++) {
	temp=0;
	temp(i,j)=z;
	temp(i,j)=-z;
	lambda[k++] = temp;
      }
  }
};

