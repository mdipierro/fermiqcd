/////////////////////////////////////////////////////////////////
/// @file mdp_matrix_field.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains class mdp_matrix_field
///
/// Licensed under GPL2 license
/// Read attached license in file mdp_license.pdf
/// This file cannot be distributed without file mdp_license.pdf
//////////////////////////////////////////////////////////////////

/// @brief a field of matrices
/// 
/// Example:
/// @verbatim
///    int box[]={10,10,10};
///    mdp_lattice lattice(3,box);
///    mdp_matrix_field h(lattice,5,5);
///    mdp_site x(lattice);
///    forallsites(x)
///       h(x)=lattice.random(x).SU(5);
/// @endverbatim
class mdp_matrix_field: public mdp_field<mdp_complex> {
public:
  int rows, columns, imax;
  mdp_matrix_field() {
    rows=columns=imax=0;
    mdp_field<mdp_complex>::reset_field();
  }
  mdp_matrix_field(mdp_matrix_field &field) {
    mdp_field<mdp_complex>::reset_field();
    rows=field.rows;
    columns=field.columns;
    imax=field.imax;
    allocate_field(field.lattice(),field.imax);
  }
  mdp_matrix_field(mdp_lattice &a, int i, int j) {
    mdp_field<mdp_complex>::reset_field();
    rows=i;
    columns=j;
    imax=i*j;
    allocate_field(a,imax);
  }
  void allocate_mdp_matrix_field(mdp_lattice &a, int i, int j) {
    deallocate_field();
    rows=i;
    columns=j;
    imax=i*j;
    allocate_field(a,imax);
  }
  mdp_matrix operator() (mdp_site x) {
    return mdp_matrix(address(x),rows,columns);
  }
  mdp_complex& operator() (mdp_site x, int i, int j) {
    return address(x)[i*columns+j];
  }
  const mdp_complex& operator() (mdp_site x, int i, int j) const {
    return address(x)[i*columns+j];
  }
};
