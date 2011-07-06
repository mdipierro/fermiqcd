#include "fermiqcd.h"              

void accumulate_c2(mdp_array<mdp_complex,1>& c2, 
		   fermi_propagator& Sl, // note props are passed by reference
		   mdp_matrix G1,
		   fermi_propagator& Sh, // note props are passed by reference
		   mdp_matrix G2) {
  // I can write this much faster but this is easier to read!
  mdp_site x(Sl.lattice());
  int a,b,c,d,i,j;
  G1=Gamma5*G1;
  G2=hermitian(Gamma5*G2);
  forallsites(x) 
    for(a=0; a<4; a++)
      for(b=0; b<4; b++)
	for(c=0; c<4; c++)
	  for(d=0; d<4; d++)
	    for(i=0; i<Sl.nc; i++)
	      for(j=0;j<Sl.nc; j++) 
		c2(x(0))+=real(Sl(x,a,b,i,j)*G1(b,c)*conj(Sh(x,d,c,i,j))*G2(d,a));
  mpi.add(c2.address(),c2.size());
};

int main(int argc, char** argv) {
  mdp.open_wormholes(argc,argv);    // open communications
  define_base_matrices("FERMILAB");

  if(argc<2) {
    mdp << 
      "usage:\n"    
      "a) to make a local pion:\n"
      "   a.out -ps -lkappa=0.123 -lcsw=0.0 -hkappa=0.123 -hcsw=0.0 -input=name\n\n"
      "b) to make a local scalar:\n"
      "   a.out -s -lkappa=0.123 -lcsw=0.0 -hkappa=0.123 -hcsw=0.0 -input=name\n\n"
      "c) to make a local vector:\n"
      "   a.out -v -lkappa=0.123 -lcsw=0.0 -hkappa=0.123 -hcsw=0.0 -input=name\n\n"
      "d) to make a local pseudovector:\n"
      "   a.out -pv -lkappa=0.123 -lcsw=0.0 -hkappa=0.123 -hcsw=0.0 -input=name\n\n"
      "e) to make a local tensor 01+02+03:\n"
      "   a.out -i0 -lkappa=0.1 -lcsw=0.0 -hkappa=0.1 -hcsw=0.0 -input=name\n\n"
      "f) to make a local tensor 12+13+23:\n"
      "   a.out -ij -lkappa=0.1 -lcsw=0.0 -hkappa=0.1 -hcsw=0.0 -input=name\n\n"
      "[name] here is the name of the gauge configuration.\n\n"
      "Other optional flags:\n"
      "  -ap=1e-8  (absolute precision)\n"
      "  -rp=1e-6  (relative precision)\n"
      "  -minres   (minimum residue, default inverter)\n"
      "  -bicgstab (stabilized biconjugate gradient)\n"
      "  -verbose  (output all steps of the inverter)\n"
      "  -clean    (recompute propagators)\n"
      "\nAttention:\n"
      "The program creates temporary propagator files in the working directory.\n"
      "These files are reused if not deleted.\n";
    exit(1);
  }
  
  mdp_field_file_header header;
  int ndim=4,nc,t,mu,nu;
  int *L;
  int verbose=false, clean=false;
  char input[1024]="";
  char output[1024]="";
  mdp_real lkappa=0.10,lcsw=0.0,hkappa=0.10,hcsw=0.0;
  mdp_real ap=1e-10,rp=1e-8;
  enum {s,ps,v,pv,i0,ij}  meson=ps;
  enum {minres, bicgstab} inverter=minres;
  
  // //////////////////////////////
  // Parsing command line arguments
  // //////////////////////////////
  for(int i=1; i<argc; i++) {
    if(strncmp(argv[i],"-verbose",8)==0)     verbose=true;
    else if(strncmp(argv[i],"-input",6)==0)  sscanf(argv[i],"-input=%s",input);
#ifndef USE_DOUBLE_PRECISION
    else if(strncmp(argv[i],"-lkappa",7)==0) sscanf(argv[i],"-lkappa=%f",&lkappa);  
    else if(strncmp(argv[i],"-lcsw",5)==0)   sscanf(argv[i],"-lcsw=%f",&lcsw);
    else if(strncmp(argv[i],"-hkappa",7)==0) sscanf(argv[i],"-hkappa=%f",&hkappa);  
    else if(strncmp(argv[i],"-hcsw",5)==0)   sscanf(argv[i],"-hcsw=%f",&hcsw);
    else if(strncmp(argv[i],"-ap",3)==0)     sscanf(argv[i],"-ap=%f",&ap);  
    else if(strncmp(argv[i],"-rp",3)==0)     sscanf(argv[i],"-rp=%f",&rp);  
#else
    else if(strncmp(argv[i],"-lkappa",7)==0) sscanf(argv[i],"-lkappa=%lf",&lkappa);  
    else if(strncmp(argv[i],"-lcsw",5)==0)   sscanf(argv[i],"-lcsw=%lf",&lcsw);
    else if(strncmp(argv[i],"-hkappa",7)==0) sscanf(argv[i],"-hkappa=%lf",&hkappa);  
    else if(strncmp(argv[i],"-hcsw",5)==0)   sscanf(argv[i],"-hcsw=%lf",&hcsw);
    else if(strncmp(argv[i],"-ap",3)==0)     sscanf(argv[i],"-ap=%lf",&ap);  
    else if(strncmp(argv[i],"-rp",3)==0)     sscanf(argv[i],"-rp=%lf",&rp);  
#endif
    else if(strncmp(argv[i],"-minres",7)==0)   inverter=minres;
    else if(strncmp(argv[i],"-bicgstab",9)==0) inverter=bicgstab;
    else if(strncmp(argv[i],"-s",2)==0)  meson=s;    
    else if(strncmp(argv[i],"-ps",3)==0) meson=ps;    
    else if(strncmp(argv[i],"-v",2)==0)  meson=v;    
    else if(strncmp(argv[i],"-pv",3)==0) meson=pv;    
    else if(strncmp(argv[i],"-i0",3)==0) meson=i0;    
    else if(strncmp(argv[i],"-ij",3)==0) meson=ij;    
    else if(strncmp(argv[i],"-clean",3)==0) clean=true;    
    else error("Wrong command line option");
  }

  // if gauge field exists determine its size (issues with single/double precision)
  if(is_file(input)) header=get_info(input);
  else error("Unable to access input gauge configuration\n");
  if(header.ndim!=4) error("Sorry, mesons only in 4D");
  nc=(int) sqrt((double) header.bytes_per_site/(4*sizeof(mdp_complex)));
  L=header.box;

  // //////////////////////////////
  // Output parameters
  // //////////////////////////////
  mdp << "=============================================\n";
  mdp << "Reading gauge configuration = " << input << '\n';
  mdp << "Number of colors = " << nc << '\n';
  mdp << "Lattice size = " << L[0] << "x" << L[1] << "x" << L[2] << "x" << L[3] << '\n'; 
  switch(meson) {
  case  s: mdp << "Computing a scalar propagator\n"; break;
  case ps: mdp << "Computing a speudo-scalar propagator\n"; break;
  case  v: mdp << "Computing a vector propagator\n"; break;
  case pv: mdp << "Computing a speudo-vector propagator\n"; break;
  case i0: mdp << "Computing a i0-tensor propagator\n"; break;
  case ij: mdp << "Computing a ij-tensor propagator\n"; break;
  }
  switch(inverter) {
  case minres: mdp << "Using minimum residue inverter\n"; break;
  case bicgstab: mdp << "Using bicg stabilized inverter\n"; break;
  }
  mdp << "Required absolute precision = " << ap << '\n';
  mdp << "Required relative precision = " << rp << '\n';
  mdp << "=============================================\n";

  
  // //////////////////////////////
  // Computation
  // //////////////////////////////

  if(!verbose)  mdp.print=false; // eventualy print off
  mdp_lattice   lattice(ndim,L); // declare lattice
  mdp_site      x(lattice);      // declare site variable
  gauge_field	U(lattice,nc);   // declare SU(3) field

  fermi_propagator Sl(lattice,nc);
  fermi_propagator Sh(lattice,nc);
  coefficients	light_quark;
  coefficients	heavy_quark;
  mdp_array<mdp_complex,1> c2(lattice.size(0));

  // load gauge field
  U.load(input);
  mdp << "average_plaquette = " << average_plaquette(U) << '\n';
  if(lcsw!=0.0 || hcsw!=0) compute_em_field(U);

  // choose inverter
  if(inverter==minres)   
    default_fermi_inverter=MinimumResidueInverter<fermi_field,gauge_field>;
  else if(inverter==bicgstab) 
    default_fermi_inverter=BiConjugateGradientStabilizedInverter<fermi_field,gauge_field>;

  // choose action implementation
  default_fermi_action=FermiCloverActionFast::mul_Q;
#ifdef SSE2
  if(nc==3) default_fermi_action=FermiCloverActionSSE2::mul_Q;
#endif

  // creating and saving light quark propagator 
  // load it if it exists already
  light_quark["kappa"]=lkappa;
  light_quark["c_{sw}"]=lcsw;
  sprintf(output,"fp_%s_l_kappa%.4f_csw%.4f",input,lkappa,lcsw);
  if(is_file(output) && !clean) 
    Sl.load(output);
  else {
    generate(Sl,U,light_quark,ap,rp);
    Sl.save(output);
  }

  // creating and saving heavy quark propagator 
  // load it if it exists already
  heavy_quark["kappa"]=hkappa;
  heavy_quark["c_{sw}"]=hcsw;
  sprintf(output,"fp_%s_l_kappa%.4f_csw%.4f",input,hkappa,hcsw);
  if(is_file(output) && !clean) 
    Sh.load(output);
  else {
    generate(Sh,U,heavy_quark,ap,rp);
    Sh.save(output);
  }

  // contract propagators with gamma matrices and make meson
  for(t=0;t<lattice.size(0);t++) c2(t)=0;
  switch(meson) {
  case s: 
    accumulate_c2(c2,Sl,Gamma1,Sh,Gamma1);
    break;
  case ps: 
    accumulate_c2(c2,Sl,Gamma5,Sh,Gamma5);
    break;
  case v: 
    for(mu=1;mu<4;mu++)  
      accumulate_c2(c2,Sl,Gamma[mu],Sh,Gamma[mu]);
    break;  
  case pv: 
    for(mu=1;mu<4;mu++) 
      accumulate_c2(c2,Sl,Gamma5*Gamma[mu],Sh,Gamma5*Gamma[mu]);
    break;  
  case i0: 
    for(mu=1;mu<4;mu++) 
      accumulate_c2(c2,Sl,Sigma[0][mu],Sh,Sigma[0][mu]);
    break;  
  case ij: 
    for(mu=1;mu<4;mu++) 
      for(nu=1;nu<4;nu++) 
	if(nu>mu) 
	  accumulate_c2(c2,Sl,Sigma[mu][nu],Sh,Sigma[mu][nu]);
    break;  
  }
  mdp.add(c2.address(),c2.size());

  // print result
  if(mdp.me()==0) mdp.print=true;
  mdp << "t, Real(c2(t)), Imag(c2(t))\n";
  for(t=0; t<lattice.size(0); t++) {
    mdp << t << ", " << real(c2(t)) << ", " << imag(c2(t)) << '\n';
  }

  mdp.close_wormholes();
  return 0;
}

