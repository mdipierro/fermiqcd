//#define BLOCKSITE 4
//#define TWISTED_BOUNDARY
#include "fermiqcd.h"

void test_wilson() {
  mdp << "START TESTING CLOVER ACTIONS\n";

  int box[]={64,10,10,10}, nc=3;
  generic_lattice lattice(4,box);
  gauge_field U(lattice,nc);
  fermi_field psi(lattice, nc);
  fermi_field chi2(lattice, nc);
  coefficients coeff;

  coeff["kappa_s"]=0.1;
  coeff["kappa_t"]=0.1;

  set_hot(U);
  set_random(psi);
  double t0,t1;
  inversion_stats stats;

  default_fermi_action=FermiCloverActionFast::mul_Q;

  default_fermi_inverter=MinimumResidueInverter<fermi_field,gauge_field>;
  t0=mpi.time();
  stats=mul_invQ(chi2,psi,U,coeff);
  t1=(mpi.time()-t0)/lattice.nvol_gl/stats.steps;
  cout << "Wilson Min Res TIME=" << t1 << endl;

  default_fermi_inverter=BiConjugateGradientStabilizedInverter<fermi_field,gauge_field>;
  t0=mpi.time();
  stats=mul_invQ(chi2,psi,U,coeff);
  t1=(mpi.time()-t0)/lattice.nvol_gl/stats.steps;
  cout << "Wilson BiCGStab TIME=" << t1 << endl;

  default_fermi_action=FermiCloverActionSSE2::mul_Q;

  default_fermi_inverter=MinimumResidueInverter<fermi_field,gauge_field>;
  t0=mpi.time();
  stats=mul_invQ(chi2,psi,U,coeff);
  t1=(mpi.time()-t0)/lattice.nvol_gl/stats.steps;
  cout << "Wilson SSE Min Res TIME=" << t1 << endl;

  default_fermi_inverter=BiConjugateGradientStabilizedInverter<fermi_field,gauge_field>;
  t0=mpi.time();
  stats=mul_invQ(chi2,psi,U,coeff);
  t1=(mpi.time()-t0)/lattice.nvol_gl/stats.steps;
  cout << "Wilson SSE BiCGStab TIME=" << t1 << endl;

}

void test_clover() {
  mdp << "START TESTING CLOVER ACTIONS\n";

  int box[]={64,6,6,6}, nc=3;
  generic_lattice lattice(4,box);
  gauge_field U(lattice,nc);
  fermi_field psi(lattice, nc);
  fermi_field chi2(lattice, nc);
  coefficients coeff;

  coeff["kappa_s"]=0.1;
  coeff["kappa_t"]=0.1;
  coeff["c_{sw}"]=1.00;

  set_hot(U);
  compute_em_field(U);
  set_random(psi);
  double t0,t1;
  inversion_stats stats;

  default_fermi_action=FermiCloverActionFast::mul_Q;

  default_fermi_inverter=MinimumResidueInverter<fermi_field,gauge_field>;
  t0=mpi.time();
  stats=mul_invQ(chi2,psi,U,coeff);
  t1=(mpi.time()-t0)/lattice.nvol_gl/stats.steps;
  cout << "Clover Min Res TIME=" << t1 << endl;

  default_fermi_inverter=BiConjugateGradientStabilizedInverter<fermi_field,gauge_field>;
  t0=mpi.time();
  stats=mul_invQ(chi2,psi,U,coeff);
  t1=(mpi.time()-t0)/lattice.nvol_gl/stats.steps;
  cout << "Clover BiCGStab TIME=" << t1 << endl;

  default_fermi_action=FermiCloverActionSSE2::mul_Q;

  default_fermi_inverter=MinimumResidueInverter<fermi_field,gauge_field>;
  t0=mpi.time();
  stats=mul_invQ(chi2,psi,U,coeff);
  t1=(mpi.time()-t0)/lattice.nvol_gl/stats.steps;
  cout << "Clover SSE Min Res TIME=" << t1 << endl;

  default_fermi_inverter=BiConjugateGradientStabilizedInverter<fermi_field,gauge_field>;
  t0=mpi.time();
  stats=mul_invQ(chi2,psi,U,coeff);
  t1=(mpi.time()-t0)/lattice.nvol_gl/stats.steps;
  cout << "Clover SSE BiCGStab TIME=" << t1 << endl;

}

void test_staggered() {
  mdp << "START TESTING STAGGERED ACTIONS\n";

  int box[]={64,6,6,6}, nc=3;
  generic_lattice lattice(4,box,default_partitioning<0>,
			  torus_topology, 0, 3);
  gauge_field U(lattice,nc);
  gauge_field V(lattice,nc);
  staggered_field psi(lattice, nc);
  staggered_field chi1(lattice, nc);
  staggered_field chi2(lattice, nc);
  coefficients coeff;
  coeff["mass"]=1.0;
  double t0, t1;
  inversion_stats stats;
  set_hot(U);
  set_random(psi);

  mdp << "ATTENTION: need to adjust asqtad coefficnets\n";

  default_staggered_action=StaggeredAsqtadActionFast::mul_Q;

  default_staggered_inverter=MinimumResidueInverter<staggered_field,gauge_field>;
  t0=mpi.time();
  stats=mul_invQ(chi2,psi,U,coeff);
  t1=(mpi.time()-t0)/lattice.nvol_gl/stats.steps;
  cout << "Staggered Min Res TIME=" << t1 << endl;

  default_staggered_inverter=BiConjugateGradientStabilizedInverter<staggered_field,gauge_field>;
  t0=mpi.time();
  stats=mul_invQ(chi2,psi,U,coeff);
  t1=(mpi.time()-t0)/lattice.nvol_gl/stats.steps;
  cout << "Staggered BiCGStab TIME=" << t1 << endl;

  default_staggered_inverter=StaggeredBiCGUML::inverter;
  t0=mpi.time();
  stats=mul_invQ(chi2,psi,U,coeff);
  t1=(mpi.time()-t0)/lattice.nvol_gl/stats.steps;
  cout << "Staggered SSE BiCGStabUML TIME=" << t1 << endl;

  default_staggered_action=StaggeredAsqtadActionSSE2::mul_Q;

  default_staggered_inverter=MinimumResidueInverter<staggered_field,gauge_field>;
  t0=mpi.time();
  stats=mul_invQ(chi2,psi,U,coeff);
  t1=(mpi.time()-t0)/lattice.nvol_gl/stats.steps;
  cout << "Staggered SSE Min Res TIME=" << t1 << endl;

  default_staggered_inverter=BiConjugateGradientStabilizedInverter<staggered_field,gauge_field>;
  t0=mpi.time();
  stats=mul_invQ(chi2,psi,U,coeff);
  t1=(mpi.time()-t0)/lattice.nvol_gl/stats.steps;
  cout << "Staggered SSE BiCGStab TIME=" << t1 << endl;

  default_staggered_inverter=StaggeredBiCGUML::inverter;
  t0=mpi.time();
  stats=mul_invQ(chi2,psi,U,coeff);
  t1=(mpi.time()-t0)/lattice.nvol_gl/stats.steps;
  cout << "Staggered SSE BiCGStabUML TIME=" << t1 << endl;
}

void test_asqtad() {
  mdp << "START TESTING STAGGERED ACTIONS\n";

  int box[]={64,6,6,6}, nc=3;
  generic_lattice lattice(4,box,default_partitioning<0>,
			  torus_topology, 0, 3);
  gauge_field U(lattice,nc);
  gauge_field V(lattice,nc);
  staggered_field psi(lattice, nc);
  staggered_field chi1(lattice, nc);
  staggered_field chi2(lattice, nc);
  coefficients coeff;
  coeff["mass"]=1.0;
  double t0, t1;
  inversion_stats stats;
  set_hot(U);
  set_random(psi);

  mdp << "ATTENTION: need to adjust asqtad coefficnets\n";

  lepage_improved_links(V,U,lepage_coefficients(0.4, "Full"));

  default_staggered_action=StaggeredAsqtadActionFast::mul_Q;

  default_staggered_inverter=MinimumResidueInverter<staggered_field,gauge_field>;
  t0=mpi.time();
  stats=mul_invQ(chi2,psi,V,coeff);
  t1=(mpi.time()-t0)/lattice.nvol_gl/stats.steps;
  cout << "Asqtad Min Res TIME=" << t1 << endl;

  default_staggered_inverter=BiConjugateGradientStabilizedInverter<staggered_field,gauge_field>;
  t0=mpi.time();
  stats=mul_invQ(chi2,psi,V,coeff);
  t1=(mpi.time()-t0)/lattice.nvol_gl/stats.steps;
  cout << "Asqtad BiCGStab TIME=" << t1 << endl;

  default_staggered_inverter=StaggeredBiCGUML::inverter;
  t0=mpi.time();
  stats=mul_invQ(chi2,psi,V,coeff);
  t1=(mpi.time()-t0)/lattice.nvol_gl/stats.steps;
  cout << "Asqtad SSE BiCGStabUML TIME=" << t1 << endl;

  default_staggered_action=StaggeredAsqtadActionSSE2::mul_Q;

  default_staggered_inverter=MinimumResidueInverter<staggered_field,gauge_field>;
  t0=mpi.time();
  stats=mul_invQ(chi2,psi,V,coeff);
  t1=(mpi.time()-t0)/lattice.nvol_gl/stats.steps;
  cout << "Asqtad SSE Min Res TIME=" << t1 << endl;

  default_staggered_inverter=BiConjugateGradientStabilizedInverter<staggered_field,gauge_field>;
  t0=mpi.time();
  stats=mul_invQ(chi2,psi,V,coeff);
  t1=(mpi.time()-t0)/lattice.nvol_gl/stats.steps;
  cout << "Asqtad SSE BiCGStab TIME=" << t1 << endl;

  default_staggered_inverter=StaggeredBiCGUML::inverter;
  t0=mpi.time();
  stats=mul_invQ(chi2,psi,V,coeff);
  t1=(mpi.time()-t0)/lattice.nvol_gl/stats.steps;
  cout << "Asqtad SSE BiCGStabUML TIME=" << t1 << endl;
}

int main(int argc, char** argv) {    
  mdp.open_wormholes(argc, argv);   
  define_base_matrices("FERMILAB");    

  test_wilson();
  test_clover();
  test_staggered();
  test_asqtad();
  
  mdp.close_wormholes();             
  return 0;                         
}
                                           
