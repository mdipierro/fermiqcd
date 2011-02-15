// Program: example06.cpp
#include "mdp.h"

float Q(float x, void *a) {
   return sin(Pi*x);
}

int main() {
   mdp_random_generator random;
   int i,N=100;
   float a,b,average=0, sigma=0.3, a_bar=1;
   for(i=0; i<N; i++) {
      a=(sigma*random.gaussian()+a_bar);
      b=random.distribution(Q);
      average+=a+b;
      cout << "average=" << average/(i+1) << '\n';
   }
   return 0;
}
