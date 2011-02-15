/////////////////////////////////////////////////////////////////
/// @file mdp_vector_field.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains class mdp_vector_field
///
/// Licensed under GPL2 license
/// Read attached license in file mdp_license.pdf
/// This file cannot be distributed without file mdp_license.pdf
//////////////////////////////////////////////////////////////////

/// @brief a field of vectors of complex numbers
/// 
/// Example:
/// @verbatim
///    int box[]={10,10,10};
///    mdp_lattice lattice(3,box);
///    mdp_vector_field h(lattice,10);
///    mdp_site x(lattice);
///    forallsites(x)
///      h(x)=0.0+0.0*I;
/// @endverbatim
class mdp_vector_field: public mdp_field<mdp_complex> {
public:
  int rows, columns, imax;
  mdp_vector_field() {
    rows=columns=imax=0;
    mdp_field<mdp_complex>::reset_field();
  }
  mdp_vector_field(mdp_vector_field &field) {
    rows=field.rows;
    columns=field.columns;
    imax=field.imax;
    allocate_field(field.lattice(),field.imax);
  }
  mdp_vector_field(mdp_lattice &a, int i) {
    rows=i;
    columns=1;
    imax=i;
    allocate_field(a,imax);
  }
  void allocate_mdp_vector_field(mdp_lattice &a, int i) {
    deallocate_field();
    rows=i;
    columns=1;
    imax=i;
    allocate_field(a,imax);
  }
  mdp_matrix operator() (mdp_site x) {
    return mdp_matrix(address(x),rows,columns);
  }
  mdp_complex& operator() (mdp_site x, int i) {
    return address(x)[i];
  }
  const mdp_complex& operator() (mdp_site x, int i) const {
    return address(x)[i];
  }
};

