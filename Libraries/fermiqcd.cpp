#include "fermiqcd.h"

void usage() {
  mdp << 
    "For help:\n"
    "  fermiqcd\n\n"
    "Examples:\n"
    "  fermiqcd -gauge:start=cold:nt=16:nx=4\n"
    "  fermiqcd -gauge:start=hot:nt=16:nx=4\n"
    "  fermiqcd -gauge:load=cold.mdp\n"
    "  fermiqcd -gauge:load=cold.mdp:steps=10:beta=5.7\n"
    "  fermiqcd -gauge:load=*.mdp -plaquette\n"
    "  fermiqcd -gauge:load=*.mdp -plaquette-vtk\n"
    "  fermiqcd -gauge:load=*.mdp -polyaov-vtk\n"
    "  fermiqcd -gauge:load=*.mdp -cool:steps=20 topcharge-vtk\n"
    "  fermiqcd -gauge:load=*.mdp -cool-vtk:steps=20\n"
    "  fermiqcd -gauge:load=*.mdp -quark:kappa=0.12:alg=minres-vtk\n"
    "  fermiqcd -gauge:load=*.mdp -quark:kappa=0.12 -pion\n"
    "  fermiqcd -gauge:load=*.mdp -quark:kappa=0.12 -pion-vtk\n"
    "Options (key = default):\n\n"
    "  -gauge\n"
    "       start = hot or cold or load\n"
    "       load = *.mdp\n"
    "       save = true\n"
    "       nt = 16\n"
    "       nx = 4\n"
    "       ny = nx\n"
    "       nz = ny\n"
    "       steps = 1\n"
    "       therm = 10\n"
    "       beta = 0\n"
    "       zeta = 1.0\n"
    "       u_t = 1.0\n"
    "       u_s = 1.0\n"
    "       prefix = \n"
    "       action = wilson or wilson_improved or wilson_sse2\n"
    "  -cool\n"
    "       alg = ape\n"
    "       alpha = 0.7\n"
    "       steps = 20\n"
    "       cooling = 10\n"
    "  -plaquette\n"
    "  -plaquette-vtk\n"
    "  -polyakov-vtk\n"
    "  -topcharge-vtk\n"
    "  -quark\n"
    "       load = false\n"
    "       save = true\n"
    "       action = clover_fast or clover_slow or clover_sse2\n"
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
    "  -pion\n"
    "  -pion-vtk\n";
    exit(0);
}

void cool(gauge_field& U, mdp_args& arguments) {
  if (arguments.get("-cool","alg","ape")=="ape")
    ApeSmearing::smear(U,
		       arguments.get("-cool","alpha",0.7),
		       arguments.get("-cool","steps",20),
		       arguments.get("-cool","cooling",10));
  else
    mdp.error_message("cooling algorithm not supported");
}

void cool_vtk(gauge_field& U, mdp_args& arguments, string filename) {
  if (arguments.get("-cool","alg","ape")=="ape")
    for(int k=0; k<arguments.get("-cool-vtk","steps",20); k++) {
      ApeSmearing::smear(U,
			 arguments.get("-cool-vtk","alpha",0.7),
			 arguments.get("-cool-vtk","steps",1),
			 arguments.get("-cool-vtk","cooling",10));
      topological_charge_vtk(U,filename+".cool"+tostring(k,2)+".vtk",0);
    }
  else
    mdp.error_message("cooling algorithm not supported");
}

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


void pretty_print(string prefix, vector<mdp_real> data) {
  if (ME==0) {
    for(int t=0; t<data.size(); t++) {
      cout << prefix << "[" << t << "] = " << data[t] << endl;
    }
  }
}

void make_quark(gauge_field &U, coefficients &gauge, coefficients &quark,
		mdp_args& arguments, string newfilename) {

  string quark_action = arguments.get("-quark","action","clover_fast|clover_slow|clover_sse2");
  string inverter = arguments.get("-quark","alg","bicgstab|minres|bicgstab-vtk|minres-vtk");
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
  int NT = U.lattice().size(0);
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
  string quarkfilename;
  mdp_complex open_prop[4][4][3][3][256];
  
  for(int a=0; a<4; a++)
    for(int i=0; i<nc; i++) {	
      psi = 0;
      if (on_which_process(U.lattice(),0,0,0,0)==ME) x.set(0,0,0,0);
      psi(x,a,i)=1;
      psi.update();
      
      inversion_vtk_prefix=newfilename+".s"+tostring(a,1)+".c"+tostring(i,1);
      quarkfilename = inversion_vtk_prefix+".quark";
      if (arguments.get("-quark","load","false|true")=="true") {
	psi.save(inversion_vtk_prefix+".quark");
      } else {
	mul_invQ(phi,psi,U,quark,abs_precision,rel_precision);
	if (arguments.get("-quark","save","true|false")=="true")
	  psi.save(inversion_vtk_prefix+".quark");
      }
      if (arguments.have("-pion") || arguments.have("-pion-vtk")) {
	if (a==0 && i==0)
	  for(int t=0; t<NT; t++)
	    pion[t]=0.0;
	forallsites(x) {
	  y.set(x(1),x(2),x(3));
	  for(int b=0; b<4; b++)
	    for(int j=0; j<nc; j++)
	      pion[x(0)]+=Q(y)=real(phi(x,b,j)*conj(phi(x,b,j)));	  
	}
	mpi.add(&pion[0],NT);	
      }
      if(arguments.have("-4quarks")) {
	for(int t=0; t<NT; t++)
	  for(int b=0; b<4; b++)
	    for(int j=0; j<nc; j++)
	      open_prop[a][b][i][j][t] = 0.0;
	forallsites(x)
	  for(int b=0; b<4; b++)
	    for(int j=0; j<nc; j++)
	      open_prop[a][b][i][j][x(0)]+=phi(x,b,j)*conj(phi(x,b,j));
	mpi.add((mdp_complex*) &open_prop[0],NT*4*4*3*3);	      	
      }
    }  
    if (arguments.have("-pion"))
      pretty_print("C2",pion);      
    if (arguments.have("-pion-vtk"))
      Q.save_vtk(inversion_vtk_prefix+".pion.vtk");
    if(arguments.have("-4quarks")) {
      mdp_matrix G1, G2;
      if(arguments.get("-4quarks","operator","5Ix5I")=="0Ix0I") {      
	G1=Gamma[0]*Gamma5;
	G2=Gamma[0]*Gamma5;
      }
      if(arguments.get("-4quarks","operator","5Ix5I")=="1Ix1I") {      
	G1=Gamma[1]*Gamma5;
	G2=Gamma[1]*Gamma5;
      }
      if(arguments.get("-4quarks","operator","5Ix5I")=="5Ix5I") {      
	G1=Gamma5*Gamma5;
	G2=Gamma5*Gamma5;
      }
      // others operators may be 0Tx0T,1Tx1T,5Tx5T,etc.
      for(int t1=0;t1<NT;t1++)
	for(int t2=0;t2<NT;t2++) {
	  mdp_real c3 = 0;
	  // manually add other contractions....
	  for(int a=0; a<4; a++)
	    for(int b=0; b<4; b++)
	      for(int c=0; c<4; c++)
		for(int d=0; d<4; d++) {
		  mdp_complex g1 = G1(b,a);
		  mdp_complex g2 = G2(d,c);
		  if(g1!=0 && g2!=0) 
		    for(int i=0; i<3; i++)
		      for(int j=0; j<3; j++)
			c3+=real(open_prop[a][b][i][i][NT-1-t1]*g1*
				 open_prop[c][d][j][j][t2]*g2);
		}
	  mdp << "C3[" << t1 << "]["<< t2 << "] = " << c3 << endl;
	}
    }
}

int main(int argc, char** argv) {
  mdp.open_wormholes(argc,argv);
  mdp_args arguments(argc,argv);
  if(!arguments.length())
    usage();
  define_base_matrices(arguments.get("-quark","matrices",
    "FERMILAB|MILC|UKQCD|Minkowsy-Dirac|Minkowsy-Chiral"));
  coefficients gauge;
  coefficients quark;
  int ndim = 4, nc=3;
  int size[4];
  string filename, newfilename, vtkfilename;
  vector<string> filenames;
  if (arguments.get("-gauge","start","load|cold|hot")=="cold") {
    int nt = arguments.get("-gauge","nt",16);
    int nx = arguments.get("-gauge","nx",4);
    int ny = arguments.get("-gauge","ny",nx);
    int nz = arguments.get("-gauge","nz",ny);
    size[0]=nt;
    size[1]=nx;
    size[2]=ny;
    size[3]=nz;
    filenames.push_back("cold.mdp");
  } else if (arguments.get("-gauge","start","load|cold|hot")=="hot") {
    int nt = arguments.get("-gauge","nt",16);
    int nx = arguments.get("-gauge","nx",4);
    int ny = arguments.get("-gauge","ny",nx);
    int nz = arguments.get("-gauge","nz",ny);
    size[0]=nt;
    size[1]=nx;
    size[2]=ny;
    size[3]=nz;
    filenames.push_back("hot.mdp");
  } else if(arguments.get("-gauge","start","load|cold|hot")=="load") {
    string pattern = arguments.get("-gauge","load","demo.mdp");
    cout << "pattern=" << pattern << endl;
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
    mdp << "no input specified\n";
    exit(0);
  }
  int nconfigs = arguments.get("-gauge","n",0);
  int nsteps = arguments.get("-gauge","steps",1);
  int ntherm = arguments.get("-gauge","therm",10);
  gauge["beta"] = arguments.get("-gauge","beta",0);
  gauge["zeta"] = arguments.get("-gauge","zeta",1.0);
  gauge["u_t"] = arguments.get("-gauge","u_t",1.0);
  gauge["u_s"] = arguments.get("-gauge","u_s",1.0);
  string prefix = arguments.get("-gauge","prefix","");
  string gauge_action = arguments.get("-gauge","action",
     "wilson|wilson_improved|wilson_sse2");
  quark["kappa"] = arguments.get("-quark","kappa",0.12);
  quark["kappa_t"] = arguments.get("-quark","kappa_t",quark["kappa"]);
  quark["kappa_s"] = arguments.get("-quark","kappa_s",quark["kappa"]);
  quark["r_t"] = arguments.get("-quark","r_t",1.0);
  quark["r_s"] = arguments.get("-quark","r_s",1.0);  
  quark["c_{SW}"] = arguments.get("-quark","c_sw",0.0);
  quark["c_E"] = arguments.get("-quark","c_E",0.0);
  quark["c_B"] = arguments.get("-quark","c_B",0.0);

  mdp_lattice lattice(ndim,size);
  gauge_field U(lattice,nc);
  for(int f=0; f<filenames.size(); f++) {
    filename = filenames[f];
    if (arguments.get("-gauge","start","load")=="cold") {
      set_cold(U);
      if (arguments.get("-gauge","save","true")=="true")
	U.save(filename);
    } else if (arguments.get("-gauge","start","load")=="hot") {
      set_hot(U);
      if (arguments.get("-gauge","save","true")=="true")
	U.save(filename);
    } else if (arguments.get("-gauge","load","")!="") {
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
	if(filename.substr(filename.size()-4,4)==".mdp")
	  newfilename = filename.substr(0,filename.size()-4) + "." + tostring(n) +".mdp";
	else if (prefix!="")
	  newfilename = prefix + "." + tostring(n) +".mdp";
	else
	  newfilename = filename + "." + tostring(n) +".mdp";
	if (arguments.get("-gauge","save","true")=="true")
	  U.save(newfilename);
      } else {
	newfilename = filename;
      }
      if (arguments.have("-plaquette")) {
	mdp << "plaquette = " << average_plaquette(U) << endl;
      }
      if (arguments.have("-cool")) {
	cool(U,arguments);
      }      
      if (arguments.have("-cool-vtk")) {
	cool_vtk(U,arguments,newfilename);
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
	make_quark(U,gauge,quark,arguments,newfilename);	
      }
    }
  }
  mdp.close_wormholes();
  return 0;
}

