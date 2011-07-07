#include "fermiqcd.h"

void plaquette_vtk(gauge_field& U, string filename) {
  mdp_field<mdp_real> q(U.lattice());
  mdp_site x(U.lattice());
  forallsites(x) if(x(0)==0) {
    q(x)=0;
    for(int mu=0; mu<4; mu++)
      for(int nu=mu+1; nu<4; nu++)
        q(x)+=real(trace(plaquette(U,x,mu,nu)));
  }  
  q.save_vtk(filename,-1);
}

void polyakov_vtk(gauge_field& U, string filename) {
  int L[3];
  L[0]=U.lattice().size(1);
  L[1]=U.lattice().size(2);
  L[2]=U.lattice().size(3);
  mdp_lattice space(3,L,
                        default_partitioning<1>,
                        torus_topology,
                        0, 1,false);
  mdp_matrix_field V(space,U.nc,U.nc);
  mdp_field<mdp_real> q(space,2);
  mdp_site x(U.lattice());
  mdp_site y(space);

  int k,mu=0,nu=1;
  mdp_complex s=0;

  forallsites(y)
    V(y)=1;
  for(int t=0; t<L[0]; t++) {
    forallsites(y) {
      x.set(t,y(0),y(1),y(2));
      V(y)=V(y)*U(x,0);
    }
  }

  forallsites(y) {
    mdp_complex z=trace(V(y));
    q(y,0)=real(z);
    q(y,1)=imag(z);
  }
  q.save_vtk(filename,-1,0,0,false);
}

float topcharge_vtk(gauge_field& U, string filename) {
  return topological_charge_vtk(U,filename,0);
}


void usage() {
  mdp << 
    "For help:\n"
    "  fermiqcd\n"
    "\nExamples:\n:"
    "  fermiqcd -cold:nt=16:nx=4\n"
    "  fermiqcd -hot:nt=16:nx=4\n"
    "  fermiqcd -load:f=cold.mdp\n"
    "  fermiqcd -load:f=cold.mdp -heatbath:steps=10:beta=5.7\n"
    "  fermiqcd -load:f=cold.mdp -heatbath:steps=10:beta=5.7\n"
    "  fermiqcd -load:f=*.mdp -plaquette\n"
    "  fermiqcd -load:f=*.mdp -plaquette-vtk\n"
    "  fermiqcd -load:f=*.mdp -polyaov-vtk\n"
    "  fermiqcd -load:f=*.mdp -cool:steps=20 topcharge-vtk\n"
    "  fermiqcd -load:f=*.mdp -quark:kappa=0.12:alg=minres-vtk\n"
    "  fermiqcd -load:f=*.mdp -quark:kappa=0.12 -pion\n"
    "  fermiqcd -load:f=*.mdp -quark:kappa=0.12 -pion-vtk\n"
    "\nOptions (key = default):\n"
    "  -cold\n"
    "       nt = 16\n"
    "       nx = 4\n"
    "       ny = nx\n"
    "       nz = ny\n"
    "  -cool\n"
    "       alg = ape\n"
    "       alpha = 0.7\n"
    "       steps = 20\n"
    "       cooling = 10\n"
    "  -gauge\n"
    "       steps = 1\n"
    "       therm = 10\n"
    "       beta = 0\n"
    "       zeta = 1.0\n"
    "       u_t = 1.0\n"
    "       u_s = 1.0\n"
    "       prefix = \n"
    "       action = wilson or wilson_improved or wilson_sse2\n"
    "  -hot\n"
    "       nt = 16\n"
    "       nx = 4\n"
    "       ny = nx\n"
    "       nz = ny\n"
    "  -load\n"
    "       f = *.mdp\n"
    "  -pion\n"
    "  -pion-vtk\n"
    "  -plaquette\n"
    "  -plaquette-vtk\n"
    "  -polyakov-vtk\n"
    "  -quark\n"
    "       action = clover_fast or clover_losw or clover_sse2\n"
    "       alg = bicgstab or minres or bictstab-vtk ot minres-vtk\n"
    "       abs_precision = 1e-12\n"
    "       rel_precision = 1e-12\n"
    "       matrices = FERMILAB or MILC or UKQCD\n"
    "                         or Minkowsy-Dirac or Minkowsy-Chiral\n"
    "       kappa = 0.12\n"
    "       kappa_t = quark.kappa\n"
    "       kappa_s = quark.kappa\n"
    "       r_t = 1.0\n"
    "       r_s = 1.0\n"
    "       c_sw = 0.0\n"
    "       c_E = 0.0\n"
    "       c_B = 0.0\n"
    "  -topcharge-vtk\n";
    exit(0);
}

void pretty_print(string prefix, vector<mdp_real> data) {
  for(int t=0; t<data.size(); t++) {
    cout << prefix << "[" << t << "] = " << data[t] << endl;
  }
}

void make_quark(gauge_field &U, coefficients &gauge, coefficients &quark1,
		mdp_args& arguments, string newfilename) {

  string quark_action = arguments.get("-quark","action","clover_fast");
  string inverter = arguments.get("-quark","alg","bicgstab");
  float abs_precision = arguments.get("-quark","abs_precision",1e-12);
  float rel_precision = arguments.get("-quark","rel_precision",1e-12);
  
  if(gauge["c_{SW}"]!=0)
    compute_em_field(U);
  if (quark_action == "clover_fast")
    default_fermi_action=FermiCloverActionFast::mul_Q;
  else if (quark_action == "clover_slow")
    default_fermi_action=FermiCloverActionSlow::mul_Q;
#if defined(SSE2)
  else if (quark_action == "clover_sse2")
    default_fermi_action=FermiCloverActionSSE2::mul_Q;
#endif	  
  else
    mdp.error_message("quark action not supported");
  if (inverter == "minres")
    default_fermi_inverter=MinRes::inverter<fermi_field,gauge_field>;
  else if (inverter == "bicgstab")
    default_fermi_inverter=BiCGStab::inverter<fermi_field,gauge_field>;
  else if (inverter == "minres-vtk")
    default_fermi_inverter=MinResVtk::inverter<fermi_field,gauge_field>;
  else if (inverter == "bicgstab-vtk")
    default_fermi_inverter=BiCGStabVtk::inverter<fermi_field,gauge_field>;
  else
    mdp.error_message("quark inverter not supported");  
  
  int nc = U.nc;
  fermi_field psi(U.lattice(),nc);
  fermi_field phi(U.lattice(),nc);
  mdp_site x(U.lattice());
  vector<mdp_real> pion(U.lattice().size(0));

  int L[3];
  L[0]=U.lattice().size(1);
  L[1]=U.lattice().size(2);
  L[2]=U.lattice().size(3);
  mdp_lattice space(3,L,
		    default_partitioning<1>,
		    torus_topology,
		    0, 1,false);
  mdp_field<float> Q(space);
  mdp_site y(space);

  for(int a=0; a<4; a++)
    for(int i=0; i<nc; i++) {	
      psi = 0;
      if (on_which_process(U.lattice(),0,0,0,0)==ME) x.set(0,0,0,0);
      psi(x,a,i)=1;
      psi.update();
      
      inversion_vtk_prefix=newfilename+".s"+tostring(a,1)+".c"+tostring(i,1);
      mul_invQ(phi,psi,U,quark1,abs_precision,rel_precision);
      psi.save(inversion_vtk_prefix+".quark");
      if (arguments.have("-pion") || arguments.have("-pion-vtk")) {
	if (a==0 && i==0)
	  for(int t=0; t<U.lattice().size(0); t++)
	    pion[t]=0.0;
	forallsites(x) {
	  y.set(x(1),x(2),x(3));
	  for(int b=0; b<4; b++)
	    for(int j=0; j<nc; j++)
	      pion[x(0)]+=Q(y)=real(phi(x,b,j)*conj(phi(x,b,j)));
	}
	mpi.add(&pion[0],U.lattice().size(0));	
      }
    }
  if (arguments.have("-pion"))
    pretty_print("pion",pion);
  if (arguments.have("-pion-vtk"))
    Q.save_vtk(inversion_vtk_prefix+".pion.vtk");
}

int main(int argc, char** argv) {
  mdp.open_wormholes(argc,argv);
  mdp_args arguments(argc,argv);
  if(!arguments.length()) usage();
  define_base_matrices(arguments.get("-quark","matrices","FERMILAB"));
  coefficients gauge;
  coefficients quark1;
  int ndim = 4, nc=3;
  int size[4];
  string filename, newfilename, vtkfilename;
  vector<string> filenames;
  if(arguments.have("-cold")) {
    int nt = arguments.get("-cold","nt",16);
    int nx = arguments.get("-cold","nx",4);
    int ny = arguments.get("-cold","ny",nx);
    int nz = arguments.get("-cold","nz",ny);
    size[0]=nt;
    size[1]=nx;
    size[2]=ny;
    size[3]=nz;
    filenames.push_back("cold.mdp");
  } else if(arguments.have("-hot")) {
    int nt = arguments.get("-hot","nt",16);
    int nx = arguments.get("-hot","nx",4);
    int ny = arguments.get("-hot","ny",nx);
    int nz = arguments.get("-hot","nz",ny);
    size[0]=nt;
    size[1]=nx;
    size[2]=ny;
    size[3]=nz;
    filenames.push_back("hot.mdp");
  } else if(arguments.have("-load")) {
    string pattern = arguments.get("-load","f","*.mdp");
    cout << pattern << endl;
    filenames = glob(pattern);    
    if (filenames.size()==0)
      mdp.error_message("No files to read");
    mdp_field_file_header header = get_info(filenames[0].c_str());
    assert(header.ndim==4);
    size[0]=header.box[0];
    size[1]=header.box[1];
    size[2]=header.box[2];
    size[3]=header.box[3];
  } else {
    mdp.error_message("no input specified");
  }
  int nconfigs = arguments.get("-gauge","n",arguments.have("-gauge")?10:0);
  int nsteps = arguments.get("-gauge","steps",1);
  int ntherm = arguments.get("-gauge","therm",10);
  gauge["beta"] = arguments.get("-gauge","beta",0);
  gauge["zeta"] = arguments.get("-gauge","zeta",1.0);
  gauge["u_t"] = arguments.get("-gauge","u_t",1.0);
  gauge["u_s"] = arguments.get("-gauge","u_s",1.0);
  string prefix = arguments.get("-gauge","prefix","");
  string gauge_action = arguments.get("-gauge","action","wlson");

  quark1["kappa"] = arguments.get("-quark","kappa",0.12);
  quark1["kappa_t"] = arguments.get("-quark","kappa_t",quark1["kappa"]);
  quark1["kappa_s"] = arguments.get("-quark","kappa_s",quark1["kappa"]);
  quark1["r_t"] = arguments.get("-quark","r_t",1.0);
  quark1["r_s"] = arguments.get("-quark","r_s",1.0);  
  quark1["c_{sw}"] = arguments.get("-quark","c_sw",0.0);
  quark1["c_E"] = arguments.get("-quark","c_E",0.0);
  quark1["c_B"] = arguments.get("-quark","c_B",0.0);

  /*
  quark2["kappa"] = arguments.get("-quark2","kappa",quark1["kappa"]);
  quark2["kappa_t"] = arguments.get("-quark2","kappa_t",quark1["happa_t"]);
  quark2["kappa_s"] = arguments.get("-quark2","kappa_s",quark1["kappa_s"]);
  quark2["r_t"] = arguments.get("-quark2","r_t",quark1["r_t"]);
  quark2["r_s"] = arguments.get("-quark2","r_s",quark1["r_s"]);  
  quark2["c_{sw}"] = arguments.get("-quark2","c_sw",quark1["c_{sw}"]);
  quark2["c_E"] = arguments.get("-quark2","c_E",quark1["c_E"]);
  quark2["c_B"] = arguments.get("-quark2","c_B",quark1["c_B"]);
  */

  mdp_lattice lattice(ndim,size);
  gauge_field U(lattice,nc);
  for(int f=0; f<filenames.size(); f++) {
    filename = filenames[f];
    if (arguments.have("-cold")) {
      set_cold(U);
      U.save(filename);
    } else if (arguments.have("-hot")) { 
      U.save(filename);
    } else if (arguments.have("-load")) { 
      U.load(filename);
    }
    for(int n=-1; n<nconfigs; n++) {
      if(n>=0) {
	int niter =(n==0)?ntherm:nsteps;
	if (gauge_action=="wilson")
	  WilsonGaugeAction::heatbath(U,gauge,niter);
	else if (gauge_action=="wilson_improved")
	  ImprovedGaugeAction::heatbath(U,gauge,niter);
#if defined(SS2)
	else if (gauge_action=="wilson_sse2")
	  ImprovedGaugeActionSSE2::heatbath(U,gauge,niter);
#endif
	else
	  mdp.error_message("gauge action not supported");	
	if (prefix!="")
	  newfilename = prefix + "." + tostring(n) +".mdp";
	else if(filename.substr(filename.size()-4,4)==".mdp")
	  newfilename = filename.substr(0,filename.size()-4) + "." + tostring(n) +".mdp";
	else
	  newfilename = filename + "." + tostring(n) +".mdp";
	U.save(newfilename);
      } else {
	newfilename = filename;
      }
      if (arguments.have("-plaquette")) {
	mdp << "plaquette = " << average_plaquette(U) << endl;
      }
      if (arguments.have("-cool")) {
	if (arguments.get("-cool","alg","ape")=="ape")
	  ApeSmearing::smear(U,
			     arguments.get("-cool","alpha",0.7),
			     arguments.get("-cool","steps",20),
			     arguments.get("-cool","cooling",10));
	else
	  mdp.error_message("cooling algorithm not supported");
      }      
      if (arguments.have("-plaquette-vtk")) {
	plaquette_vtk(U,newfilename+".plaquette.vtk");
      }      
      if (arguments.have("-polyakov-vtk")) {
	polyakov_vtk(U,newfilename+".polyakov.vtk");
      }      
      if (arguments.have("-topcharge-vtk")) {
	float tc = topcharge_vtk(U,newfilename+".topcharge.vtk");
	mdp << "topcharge = " << tc << endl;
      }      
      if (arguments.have("-quark")) {
	make_quark(U,gauge,quark1,arguments,newfilename);
      }
    }
  }
  mdp.close_wormholes();
  return 0;
}
