/////////////////////////////////////////////////////////////////
/// @file mdp_field_test.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains a sample test (main) function
///
/// Licensed under GPL2 license
/// Read attached license in file mdp_license.pdf
/// This file cannot be distributed without file mdp_license.pdf
//////////////////////////////////////////////////////////////////

/// For debugging purposes only
bool mdp_field_test(int argc, char** argv) {
  mpi.open_wormholes(argc, argv);
  
  int box[]={4,4,4,4};
  mdp_lattice      lattice(4,box);
  mdp_matrix_field M(lattice,3,3);
  mdp_site         x(lattice);
  double counter=0;

  forallsites(x)
    M(x)=M.lattice().random(x).SU(3);
  
  forallsites(x)
    counter+=real(trace(M(x)*inv(M(x))))/3;

  mpi.add(counter);
  counter/=lattice.size();

  assert((counter-1)<mdp_precision);
  printf("lattice and field ...test passed.\n");

  mpi.close_wormholes();
  return (mdp_int) counter;
}
