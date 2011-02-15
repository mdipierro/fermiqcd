// Program: example04.cpp
#include "mdp.h"

mdp_array<mdp_matrix,3> initialize() {
  mdp_array<mdp_matrix,3> d(20,20,20);
  int i,j,k;	
  for(i=0; i<20; i++)
    for(j=0; j<20; j++)
      for(k=0; k<20; k++) {
        d(i,j,k).dimension(2,2);	
        d(i,j,k)(0,0)=k;
        d(i,j,k)(0,1)=i;
        d(i,j,k)(1,0)=j;
        d(i,j,k)(1,1)=k;
      }
  return(d);
}

mdp_array<mdp_matrix,3> f(mdp_array<mdp_matrix,3>& c) {
  mdp_array<mdp_matrix,3> d(c.size(0),c.size(1),c.size(2));
  int i,j,k;
  for(i=0; i<c.size(0); i++)
    for(j=0; j<c.size(1); j++)
      for(k=0; k<c.size(2); k++)
        d(i,j,k)=sin(c(i,j,k));
  return(d);
}

int main() {
  mdp_array<mdp_matrix,3> a, b;
  a=initialize();
  b=f(a);

  int i=1, j=2, k=3;
  cout << a(i,j,k) << '\n';
  cout << b(i,j,k) << '\n';
  return 0;
}
