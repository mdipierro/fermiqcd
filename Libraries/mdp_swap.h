/////////////////////////////////////////////////////////////////
/// @file mdp_swap.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains swap function
///
/// Licensed under GPL2 license
/// Read attached license in file mdp_license.pdf
/// This file cannot be distributed without file mdp_license.pdf
//////////////////////////////////////////////////////////////////

template<class T>
void swap(T &a, T &b) {
  T c;
  c=a;
  a=b;
  b=c;
}

template<class T>
void swap(T* a, T* b, int n) {
  int i;
  T c;
  for(i=0; i<n; i++) {
    c=a[i];
    a[i]=b[i];
    b[i]=c;
  }
}
