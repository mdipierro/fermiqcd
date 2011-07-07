//#define BLOCKSITE 4
//#define TWISTED_BOUNDARY
#include "fermiqcd.h"

void test_gauge(int nt, int nx, int ny, int nz, int nc) {
  mdp << "START TESTING GAUGE ACTION\n";
  int box[]={nt,nx,ny,nz};
  generic_lattice lattice(4,box);
  gauge_field U(lattice,nc);
  coefficients coeff;
  coeff["beta"]=6.0;
  set_cold(U);
  for(uint i=0; i<10; i++) {
    WilsonGaugeAction::heatbath(U,coeff);
    mdp << "plaquette = " << average_plaquette(U) << '\n';
  }
  GaugeFixing::fix(U);
  mdp << "plaquette after gauge fixing = " << average_plaquette(U) << '\n';
  mdp << "END TESTING GAUGE ACTION\n";
}

void test_gauge_improved(int nt, int nx, int ny, int nz, int nc) {
  mdp << "START TESTING IMPROVED GAUGE ACTION\n";
  int box[]={nt,nx,ny,nz};
  generic_lattice lattice(4,box,default_partitioning<0>,
			  torus_topology, 0, 3);
  gauge_field U(lattice,nc);
  coefficients coeff;
  coeff["beta"]=6.0;
  set_hot(U);
  for(uint i=0; i<10; i++) {
    ImprovedGaugeAction::heatbath(U,coeff,1,"MILC");
    mdp << "plaquette = " << average_plaquette(U) << '\n';
  }
  mdp << "END TESTING IMPROVED GAUGE ACTION\n";
}


void test_fermi(int nt, int nx, int ny, int nz, int nc) {
  mdp << "START TESTING CLOVER ACTIONS\n";
  int box[]={nt,nx,ny,nz};
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

  mdp << "Using FermiCloverActionSlow\n";
  default_fermi_action=FermiCloverActionSlow::mul_Q;
  mul_Q(chi1,psi,U,coeff);  

  mdp << "Using FermiCloverActionFast\n";
  default_fermi_action=FermiCloverActionFast::mul_Q;
  mul_Q(chi2,psi,U,coeff);
  if(check_differences(chi1, chi2)>1e-5) { mdp << "FAILURE\n"; exit(1); }
#if defined(SSE2) 
  mdp << "Using FermiCloverActionSSE2\n";
  default_fermi_action=FermiCloverActionSSE2::mul_Q;
  mul_Q(chi2,psi,U,coeff);
  if(check_differences(chi1, chi2)>1e-5) { mdp << "FAILURE\n"; exit(1); }
#endif  
  mdp << "Using MinRes::inverter\n";
  default_fermi_inverter=MinRes::inverter<fermi_field,gauge_field>;
  mul_invQ(chi1,chi2,U,coeff);
  if(check_differences(chi1, psi)>1e-5) { mdp << "FAILURE\n"; exit(1); }

  mdp << "Using BiCGStab::inverter\n";
  default_fermi_inverter=BiCGStab::inverter<fermi_field,gauge_field>;
  mul_invQ(chi1,chi2,U,coeff);
  if(check_differences(chi1, psi)>1e-5) { mdp << "FAILURE\n"; exit(1); }

  mdp << "END TESTING CLOVER ACTIONS\n";
}



void test_staggered(int nt, int nx, int ny, int nz, int nc) {
  mdp << "START TESTING STAGGERED ACTIONS\n";
  int box[]={nt,nx,ny,nz};
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

  mdp << "ATTENTION: need to adjust asqtad coefficnets\n";
  lepage_improved_links(V,U,lepage_coefficients(0.4, "Full"));

  mdp << "Using StaggeredAsqtadActionSlow\n";
  default_staggered_action=StaggeredAsqtadActionSlow::mul_Q;
  mul_Q(chi1,psi,V,coeff);  

  mdp << "Using StaggeredAsqtadActionFast\n";
  default_staggered_action=StaggeredAsqtadActionFast::mul_Q;
  mul_Q(chi2,psi,V,coeff);
  if(check_differences(chi1, chi2)>1e-5) { mdp << "FAILURE\n"; exit(1); }
#if defined(SSE2) 
  mdp << "Using FermiCloverActionSSE2\n";
  default_staggered_action=StaggeredAsqtadActionSSE2::mul_Q;
  mul_Q(chi2,psi,V,coeff);
  if(check_differences(chi1, chi2)>1e-5) { mdp << "FAILURE\n"; exit(1); }
#endif  
  mdp << "Using MinRes::inverter\n";
  default_staggered_inverter=MinRes::inverter<staggered_field,gauge_field>;
  mul_invQ(chi1,chi2,V,coeff);
  if(check_differences(chi1, psi)>1e-5) { mdp << "FAILURE\n"; exit(1); }

  mdp << "Using BiCGStab::inverter\n";
  default_staggered_inverter=BiCGStab::inverter<staggered_field,gauge_field>;
  mul_invQ(chi1,chi2,V,coeff);
  if(check_differences(chi1, psi)>1e-5) { mdp << "FAILURE\n"; exit(1); }

  mdp << "Using StaggeredBiCGUML::inverter\n";
  default_staggered_inverter=StaggeredBiCGUML::inverter;
  mul_invQ(chi1,chi2,V,coeff);
  if(check_differences(chi1, psi)>1e-5) { mdp << "FAILURE\n"; exit(1); }

  mdp << "END TESTING STAGGERED ACTIONS\n";
}

void test_dwfermi(int nt, int nx, int ny, int nz, int nc) {
  mdp << "START TESTING DOMAIN WALL ACTIONS\n";
  int box[]={nt,nx,ny,nz};
  int L5=10;
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

  mdp << "Using DWFermiActionSlow\n";
  default_dwfermi_action=DWFermiActionSlow::mul_Q;
  mul_Q(chi2,psi,U,coeff);  

  mdp << "Using MinRes::inverter\n";
  default_dwfermi_inverter=MinRes::inverter<dwfermi_field,gauge_field>;
  mul_invQ(chi1,chi2,U,coeff);
  if(check_differences(chi1, psi)>1e-5) { mdp << "FAILURE\n"; exit(1); }

  mdp << "Using BiCGStab::inverter\n";
  default_dwfermi_inverter=BiCGStab::inverter<dwfermi_field,gauge_field>;
  mul_invQ(chi1,chi2,U,coeff);
  if(check_differences(chi1, psi)>1e-5) { mdp << "FAILURE\n"; exit(1); }

  mdp << "END TESTING DOMAIN WALL ACTIONS\n";
}

int main(int argc, char** argv) {    
  mdp.open_wormholes(argc, argv);   
  define_base_matrices("FERMILAB");    

  test_gauge(4,4,4,4,3);
  test_gauge_improved(4,4,4,4,3);
  test_fermi(4,4,4,4,3);
  test_staggered(4,4,4,4,3);
  test_dwfermi(4,4,4,4,3);
  
  mdp.close_wormholes();             
  return 0;                         
}
                                           
