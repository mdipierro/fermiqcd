/////////////////////////////////////////////////////////////////
/// @file mdp_permutations.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Functions to compute permutations
///
/// Licensed under GPL2 license
/// Read attached license in file mdp_license.pdf
/// This file cannot be distributed without file mdp_license.pdf
//////////////////////////////////////////////////////////////////

// this is my favourite piece of code...


// this is just n!
mdp_int mdp_permutations(int n) {
  mdp_int a=1;
  for(;n;n--) a*=n;
  return a;
}

// this sorts the first k elements of map[] assuming 
// the fisrt k-1 are already sorted
void mdp_permutation_sort(int map[], int k) {
  int i, tmp;
  for(i=k-1; i>=0; i--)
    if(map[i]>map[i+1]) {
      tmp=map[i];
      map[i]=map[i+1];
      map[i+1]=tmp;
    }
}

/// Returns j-th element of the k-th permutations of n numbers
/// For example if n=4
/// [0123] k=0 
/// [0132] k=1
/// ...
/// [3210] k=23
/// Returns -1 on error when (i>n || k>n_permutations(n))
int mdp_permutation(int n, int k, int i) {
  int* map=new int[i+1];
  int j,l,m;

  if(i>n || k>mdp_permutations(n)) return -1;

  for(j=0; j<=i; j++) {
    map[j] = (k % mdp_permutations(n-j)) / mdp_permutations(n-1-j);
    mdp_permutation_sort(map, j-1);
    for(l=0; l<j; l++)
      if(map[l]<=map[j]) map[j]++;    
    if(i==j) {
      m=map[j];
      delete[] map;
      return m; 
    }
  }
  delete[] map;
  return -1;
}
