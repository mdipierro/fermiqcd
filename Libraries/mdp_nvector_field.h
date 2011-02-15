/////////////////////////////////////////////////////////////////
/// @file mdp_nvector_field.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains class mdp_nvector_field (deprecated)
///
/// Licensed under GPL2 license
/// Read attached license in file mdp_license.pdf
/// This file cannot be distributed without file mdp_license.pdf
//////////////////////////////////////////////////////////////////

/// @brief field of vectors of vectors (DEPRECATED)
class mdp_nvector_field: public mdp_field<mdp_complex> {
public:
  uint rows, columns, imax, imax2;
  mdp_nvector_field() {
    rows=columns=imax=imax2=0;
    mdp_field<mdp_complex>::reset_field();
  }
  mdp_nvector_field(mdp_nvector_field &field) {
    rows=field.rows;
    columns=field.columns;
    imax=field.imax;
    imax2=field.imax2;
    allocate_field(field.lattice(),field.imax);
  }
  mdp_nvector_field(mdp_lattice &a, int n, int i) {
    rows=i;
    columns=1;
    imax=i*n;
    imax2=i;
    allocate_field(a,imax);
  }
  void allocate_mdp_nvector_field(mdp_lattice &a, int n, int i) {
    deallocate_field();
    rows=i;
    columns=1;
    imax=i*n;
    imax2=i;
    allocate_field(a,imax);
  }
  mdp_matrix operator() (mdp_site x,int n) {
    return mdp_matrix(address(x,n*imax2),rows,columns);
  }
  mdp_complex & operator() (mdp_site x, int n, int i) {
    return address(x,n*imax2)[i];
  }
  const mdp_complex & operator() (mdp_site x, int n, int i) const {
    return address(x,n*imax2)[i];
  }
};


