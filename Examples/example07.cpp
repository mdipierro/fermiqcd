// Program: example07.cpp
#include "mdp.h"
 
const int n=100; 

float f1(float *x, void *a) { return x[0]/x[1]; }

float f2(float *x, void *a) { return x[0]*x[1]; }

int main() {
   mdp_matrix A;
   mdp_jackboot jb(n,2);
   int i;
   for(i=0; i<n; i++) { 
      A=Random.SU(6)+Random.gaussian();
      jb(i,0)=real(det(inv(A))); 
      jb(i,1)=real(det(A)); 
   }
   jb.f=f1;
   cout << "Result x[0]/x[1] = " << jb.mean() << '\n';
   cout << "Jackknife error  = " << jb.j_err() << '\n';
   cout << "Bootstrap error  = " << jb.b_err(100) << '\n';
   jb.f=f2;
   cout << "Result x[0]*x[1] = " << jb.mean() << '\n';
   cout << "Jackknife error  = " << jb.j_err() << '\n';
   cout << "Bootstrap error  = " << jb.b_err(100) << '\n';
   return 0;
}
