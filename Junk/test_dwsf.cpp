//#define BLOCKSITE 4
//#define TWISTED_BOUNDARY
#include "fermiqcd.h"

void test_gauge() {
  int box[]={4,4,4,4}, nc=3;
  generic_lattice lattice(4,box);
  gauge_field U(lattice,nc);
  coefficients coeff;
  coeff["beta"]=6.0;
  set_cold(U);
  set_hot(U);
  for(uint i=0; i<10; i++) {
    WilsonGaugeAction::heatbath(U,coeff);
    mpi << "plaquette = " << average_plaquette(U) << '\n';
  }
  GaugeFixing::fix(U);
  mpi << "plaquette after gauge fixing = " << average_plaquette(U) << '\n';

}

void test_gauge_improved() {
  int box[]={4,4,4,4}, nc=3;
  generic_lattice lattice(4,box,default_partitioning<0>,
			  torus_topology, 0, 3);
  gauge_field U(lattice,nc);
  coefficients coeff;
  coeff["beta"]=6.0;
  set_hot(U);
  for(uint i=0; i<10; i++) {
     ImprovedGaugeAction::heatbath(U,coeff,1,"MILC");
     mpi << "plaquette = " << average_plaquette(U) << '\n';
  }
}


void test_fermi() {
  int box[]={4,4,4,4}, nc=3;

  mpi << "\n\nTEST FERMI FIELDS\n\n";

  generic_lattice lattice(4,box);
  gauge_field U(lattice,nc);
  fermi_field psi(lattice, nc);
  fermi_field chi1(lattice, nc);
  fermi_field chi2(lattice, nc);
  coefficients coeff;

  coeff["kappa_s"]=0.1;
  coeff["kappa_t"]=0.1;
  coeff["c_{sw}"]=1.00;

  set_hot(U);
  compute_em_field(U);
  set_random(psi);

  default_fermi_inverter=MinimumResidueInverter<fermi_field,gauge_field>;
  default_fermi_action=FermiCloverActionSlow::mul_Q;
  mul_invQ(chi1,psi,U,coeff);
  
  default_fermi_inverter=MinimumResidueInverter<fermi_field,gauge_field>;
  default_fermi_action=FermiCloverActionFast::mul_Q;
  mul_invQ(chi2,psi,U,coeff);

  mpi << "\n\nCheching that CloverActionFast and CloverActionSlow agree\n\n";
  check_differences(chi1, chi2);

#if defined(SSE2) 
  default_fermi_inverter=MinimumResidueInverter<fermi_field,gauge_field>;
  default_fermi_action=FermiCloverActionSSE2::mul_Q;
  mul_invQ(chi2,psi,U,coeff);

  mpi << "\n\nCheching that CloverActionSlow and CloverActionSSE2 agree\n\n";
  check_differences(chi1, chi2);
#endif

  mpi << "\n\nCheching that inversion was correct\n\n";
  mul_Q(chi1,chi2,U,coeff);
  check_differences(psi, chi1);

  default_fermi_inverter=BiConjugateGradientStabilizedInverter<fermi_field,gauge_field>;
  default_fermi_action=FermiCloverActionSlow::mul_Q;
  mul_invQ(chi1,psi,U,coeff);

  mpi << "\n\nCheching that MinimumResidue and BiConjugateGradientStabilized agree\n\n";
  check_differences(chi1, chi2);
  }



void test_staggered() {
  int box[]={4,4,4,4}, nc=3;

  mpi << "\n\nTEST STAGGERED FIELDS\n\n";

  generic_lattice lattice(4,box,default_partitioning<0>,
			  torus_topology, 0, 3);
  gauge_field U(lattice,nc);
  gauge_field V(lattice,nc);
  staggered_field psi(lattice, nc);
  staggered_field chi1(lattice, nc);
  staggered_field chi2(lattice, nc);
  coefficients coeff;
  coeff["mass"]=1.0;

  set_hot(U);
  set_random(psi);

  mpi << "ATTENTION: need to adjust asqtad coefficnets\n";

  lepage_improved_links(V,U,lepage_coefficients(0.4, "Full"));

  default_staggered_inverter=BiConjugateGradientStabilizedInverter<staggered_field,gauge_field>;
  default_staggered_action=StaggeredAsqtadActionSlow::mul_Q;
  mul_invQ(chi2,psi,V,coeff);

  default_staggered_inverter=BiConjugateGradientStabilizedInverter<staggered_field,gauge_field>;
  default_staggered_action=StaggeredAsqtadActionFast::mul_Q;
  mul_invQ(chi1,psi,V,coeff);

  mpi << "\n\nCheching that AsqtadActionSlow and AsqtadActionFast agree\n\n";
  check_differences(chi1, chi2);

#if defined(SSE2)
  default_staggered_inverter=BiConjugateGradientStabilizedInverter<staggered_field,gauge_field>;
  default_staggered_action=StaggeredAsqtadActionSSE2::mul_Q;
  mul_invQ(chi2,psi,V,coeff);
  mpi << "\n\nCheching that AsqtadActionSlow and AsqtadActionSSE2 agree\n\n";
  check_differences(chi1, chi2);
#endif

  mpi << "\n\nCheching that inversion(s) were correct\n\n";
  mul_Q(chi1,chi2,V,coeff);
  check_differences(psi, chi1);

 
  default_staggered_inverter=StaggeredBiCGUML::inverter;
  default_staggered_action=StaggeredAsqtadActionFast::mul_Q;
  mul_invQ(chi1,psi,V,coeff);

  mpi << "\n\nCheching that MinimumResidue and BiConjugateGradientStabilized agree\n\n";
  check_differences(chi1, chi2);
}

void test_dwfermi() {
  int box[]={4,4,4,4}, nc=3, L5=5;

  mpi << "\n\nTEST DWFIELDS FIELDS\n\n";

  generic_lattice lattice(4,box);
  gauge_field U(lattice,nc);
  dwfermi_field psi(lattice, L5, nc);
  dwfermi_field chi1(lattice, L5, nc);
  dwfermi_field chi2(lattice, L5, nc);
  coefficients coeff;
  coeff["m_f"]=1.0;
  coeff["m_5"]=1.0;

  set_hot(U);
  set_random(psi);
  compute_em_field(U);
  default_dwfermi_inverter=MinimumResidueInverter<dwfermi_field,gauge_field>;
  default_dwfermi_action=DWFermiActionSlow::mul_Q;
  mul_invQ(chi1,psi,U,coeff);

  default_dwfermi_inverter=BiConjugateGradientStabilizedInverter<dwfermi_field,gauge_field>;
  default_dwfermi_action=DWFermiActionSlow::mul_Q;
  mul_invQ(chi2,psi,U,coeff);

  mpi << "\n\nCheching that MinimumResidue and BiConjugateGradientStabilized agree\n\n";
  check_differences(chi1, chi2);
  mul_Q(chi1,chi2,U,coeff);
  mpi << "\n\nCheching that inversion was correct\n\n";
  check_differences(psi, chi1);
}

/*
  WORK IN PROGRESS

void test_sdwf() {
  int box[]={4,4,4,4}, nc=3, L5=5;

  if(ME==0) printf("\n\nTEST SDWF FIELDS\n\n");

  generic_lattice lattice(4,box);
  gauge_field U(lattice,nc);
  sdwf_field psi(lattice, L5, nc);
  sdwf_field chi1(lattice, L5, nc);
  sdwf_field chi2(lattice, L5, nc);
  coefficients coeff;

  coeff["m_0"]=0.0;
  coeff["m_f"]=1.0;
  coeff["m_5"]=0.3;

  set_cold(U);
  set_random(psi);

  default_sdwf_inverter=BiConjugateGradientStabilizedInverter<sdwf_field,gauge_field>;
  default_sdwf_action=SDWFActionSlow::mul_Q;
  compute_swirls_field(U);
  mul_invQ(chi2,psi,U,coeff);

  mul_Q(chi1,chi2,U,coeff);
  if(ME==0) printf("\n\nCheching that inversion was correct\n\n");
  check_differences(psi, chi1);
}
*/

int main(int argc, char** argv) {    
  mpi.open_wormholes(argc, argv);   
  define_base_matrices("FERMILAB");    

  test_gauge();
  test_gauge_improved();
  test_fermi();
  test_staggered();
  test_dwfermi();
  // test_sdwf();
  
  mpi.close_wormholes();             
  return 0;                         
}
                                           
