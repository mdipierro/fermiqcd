// Program: example03.cpp
#include "mdp.h"  

mdp_matrix cube(Matrix X) {  
   mdp_matrix Y;
   Y=X*X*X;
   return Y;
}

int main() {
   mdp_matrix A,B;
   A=Random.SU(3);
   B=cube(A)*exp(A)+inv(A);
   cout << A << '\n';
   cout << B << '\n';
   return 0;
}
