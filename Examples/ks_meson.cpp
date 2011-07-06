#include "fermiqcd.h"              

int main(int argc, char** argv) {
  mdp.open_wormholes(argc,argv);    // open communications
  define_base_matrices("FERMILAB");

  if(argc<2) {
    mdp << 
      "usage:\n"    
      "a) to make a goldstone pion (wall-local) using plain ks fermions:\n"
      "   a.out -meson=5x5 -lmass=0.1 -hmass=0.1 -input=name\n\n"
      "b) to make a goldstone pion (wall-local) using asqtad ks fermions:\n"
      "   a.out -meson=5x5 -lmass=0.1 -hmass=0.1 -u0=1.0 -input=name\n\n"
      "[name] here is the name of the gauge configuration.\n\n"
      "To make other meson replace 5x5 with any of the following:\n"
      "  5x5, 05x05, 5x35, 05x12, 5x03, 05x3, 5x0\n\n"
      "Other optional flags:\n"
      "  -ap=1e-8  (absolute precision)\n"
      "  -bicguml  (uml stabilized biconjugate gradient, defualt)\n"
      "  -minres   (minimum residue)\n"
      "  -bicgstab (stabilized biconjugate gradient)\n"
      "  -verbose  (output all steps of the inverter)\n";
    exit(1);
  }
  
  mdp_field_file_header header;
  int ndim=4,nc,t,i,j;
  int *L;
  int verbose=false;
  char meson[1024]="5x5";
  char input[1024]="";
  mdp_real lmass=0.10,hmass=0.10,u0=-1.0;
  mdp_real ap=1e-10;
  enum {minres, bicgstab, bicguml} inverter=bicguml;
  
  // //////////////////////////////
  // Parsing command line arguments
  // //////////////////////////////
  for(int i=1; i<argc; i++) {
    if(strncmp(argv[i],"-verbose",8)==0)     verbose=true;
    else if(strncmp(argv[i],"-input",6)==0)  sscanf(argv[i],"-input=%s",input);
    else if(strncmp(argv[i],"-meson",6)==0)  sscanf(argv[i],"-meson=%s",meson);
#ifndef USE_DOUBLE_PRECISION
    else if(strncmp(argv[i],"-lmass",6)==0)  sscanf(argv[i],"-lmass=%f",&lmass); 
    else if(strncmp(argv[i],"-hmass",6)==0)  sscanf(argv[i],"-hmass=%f",&hmass); 
    else if(strncmp(argv[i],"-u0",3)==0)     sscanf(argv[i],"-u0=%f",&u0);  
    else if(strncmp(argv[i],"-ap",3)==0)     sscanf(argv[i],"-ap=%f",&ap);  
#else
    else if(strncmp(argv[i],"-lmass",6)==0)  sscanf(argv[i],"-lmass=%lf",&lmass); 
    else if(strncmp(argv[i],"-hmass",6)==0)  sscanf(argv[i],"-hmass=%lf",&hmass); 
    else if(strncmp(argv[i],"-u0",3)==0)     sscanf(argv[i],"-u0=%lf",&u0);  
    else if(strncmp(argv[i],"-ap",3)==0)     sscanf(argv[i],"-ap=%lf",&ap);  
#endif
    else if(strncmp(argv[i],"-bicguml",7)==0)  inverter=bicguml;
    else if(strncmp(argv[i],"-minres",7)==0)   inverter=minres;
    else if(strncmp(argv[i],"-bicgstab",9)==0) inverter=bicgstab;
    else error("Wrong command line option");
  }

  //
  // Note is u0<=0 ks fermions, if u0>0 asqtad ks
  //

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
  if(u0<0) mdp << "Creating a meson " << meson << " using ks fermions\n";
  else     mdp << "Creating a meson " << meson << " using asqtad ks fermions (u0=" <<u0<< ")\n";
  mdp << "light staggered quark mass = " << lmass << '\n';
  mdp << "heavy staggered quark mass = " << hmass << '\n';
  switch(inverter) {
  case bicguml: mdp << "Using bicg uml (milc) inverter\n"; break;
  case minres: mdp << "Using minimum residue inverter\n"; break;
  case bicgstab: mdp << "Using bicg stabilized inverter\n"; break;
  }
  mdp << "Required absolute precision = " << ap << '\n';
  mdp << "=============================================\n";
    
  // //////////////////////////////
  // Computation
  // //////////////////////////////

  if(!verbose)  mdp.print=false; // eventualy print off
  // declare lattice (note the next_next=3 for long links)
  if(u0<0) i=1; else i=3;
  mdp_lattice   lattice(ndim,L,default_partitioning0,torus_topology,0,i); 
  mdp_site      x(lattice);      // declare site variable
  gauge_field	U(lattice,nc);   // declare SU(3) field
  gauge_field	V(lattice,nc);   // declare SU(3) field (for asqtad action)
  coefficients  light_quark;
  coefficients  heavy_quark;
  mdp_matrix    prop;

  // load gauge field
  U.load(input);
  mdp << "average plaquette = " << average_plaquette(U) << '\n';
  if(u0>0) 
    lepage_improved_links(V,U,lepage_coefficients(pow(u0,4), "Full"));
  else
    V=U;

  // choose inverter
  if(inverter==bicguml)   
    default_staggered_inverter=StaggeredBiCGUML::inverter;
  else if(inverter==minres)   
    default_staggered_inverter=MinimumResidueInverter<staggered_field,gauge_field>;
  else if(inverter==bicgstab) 
    default_staggered_inverter=BiConjugateGradientStabilizedInverter<staggered_field,gauge_field>;

  // choose action implementation
  default_staggered_action=StaggeredAsqtadActionFast::mul_Q;
#ifdef SSE2
  if(nc==3) default_staggered_action=StaggeredAsqtadActionSSE2::mul_Q;
#endif

  mdp_matrix G1; // spin structure of the meson
  mdp_matrix G2; // flavour structure ot the meson

  if(string(meson)=="5x5") {
    // this is the Goldston pion !!!
    G1=Gamma5;
    G2=Gamma5;
  } else if(string(meson)=="05x05") {
    G1=Gamma[0]*Gamma5;
    G2=Gamma[0]*Gamma5;
  } else if(string(meson)=="5x35") {
    G1=Gamma5;
    G2=Gamma[3]*Gamma5;
  } else if(string(meson)=="05x12") {
    G1=Gamma[0]*Gamma5;
    G2=Gamma[1]*Gamma[2];
  } else if(string(meson)=="5x03") {
    G1=Gamma5;
    G2=Gamma[0]*Gamma[3];
  } else if(string(meson)=="05x3") {
    G1=Gamma[0]*Gamma5;
    G2=Gamma[3];
  } else if(string(meson)=="5x0") {
    G1=Gamma5;
    G2=Gamma[0];
  } else error("meson non recognized!");

  light_quark["mass"]=lmass;
  heavy_quark["mass"]=hmass;
  // look inside tthe function below in file fermiqcd_staggered_mesons.h
  prop=make_meson(U,V,G1,G2,light_quark,heavy_quark,wall_source,local_source,ap);

  // print result
  if(mdp.me()==0) mdp.print=true;
  mdp << "t, Real(c2(t)), Imag(c2(t))\n";
  for(t=0; t<lattice.size(0); t+=2) {
    mdp << t << ", " << real(prop(0,t)) << ", " << imag(prop(0,t)) << '\n';
  }

  mdp.close_wormholes();
  return 0;
}

