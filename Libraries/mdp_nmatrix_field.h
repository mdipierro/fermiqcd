/////////////////////////////////////////////////////////////////
/// @file mdp_nmatrix_field.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains class mdp_nmatrix_field
///
/// Licensed under GPL2 license
/// Read attached license in file mdp_license.pdf
/// This file cannot be distributed without file mdp_license.pdf
//////////////////////////////////////////////////////////////////

/// @brief field of vectors of matrices
///
/// Example:
/// @verbatim
///    int box[]={10,10,10};
///    mdp_lattice lattice(3,box);
///    mdp_nmatrix_field h(lattice,10,3,3);
///    mdp_site x(lattice);
///    forallsites(x)
///      for(int i=0; i<10; i++)
///        h(x,i)=lattice.random(x).SU(3);
/// @endverbatim
class mdp_nmatrix_field: public mdp_field<mdp_complex> {
public:
  uint rows, columns, matrices, imax, imax2;
  mdp_nmatrix_field() {
    rows=columns=matrices=imax=imax2=0;
    mdp_field<mdp_complex>::reset_field();
  }
  mdp_nmatrix_field(mdp_nmatrix_field &field) {
    mdp_field<mdp_complex>::reset_field();
    rows=field.rows;
    columns=field.columns;
    matrices=field.matrices;
    imax=field.imax;
    imax2=field.imax2;
    allocate_field(field.lattice(),field.imax);
  }
  /// declares a field object that at each site as vector of n ixj matrices
  mdp_nmatrix_field(mdp_lattice &a, int n, int i, int j) {
    mdp_field<mdp_complex>::reset_field();
    rows=i;
    columns=j;
    matrices=n;
    imax=i*j*n;
    imax2=i*j;
    allocate_field(a,imax);
  }
  /// dynamically allocates a field object that at each site as vector of n ixj matrices
  void allocate_mdp_nmatrix_field(mdp_lattice &a, int n, int i, int j) {
    deallocate_field();
    rows=i;
    columns=j;
    matrices=n;
    imax=i*j*n;
    imax2=i*j;
    allocate_field(a,imax);
  }  
  /// returns the n-th matrix stored at site x
  mdp_matrix operator() (mdp_site x, int n) {
    return mdp_matrix(address(x,n*imax2),rows,columns);
  }
  /// returns the (i,j) component of the n-th matrix stored at site x
  mdp_complex& operator() (mdp_site x, int n, int i, int j) {
    return address(x,n*imax2)[i*columns+j];
  }
  const mdp_complex& operator() (mdp_site x, int n, int i, int j) const {
    return address(x,n*imax2)[i*columns+j];
  }
};
