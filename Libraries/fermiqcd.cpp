#include "fermiqcd.h"
/*
./a.out -gauge:start=cold:nt=32:nx=8 -quark:kappa=0.115:source_point=center:load=true -pion:vtk=true -current-static:source_gamma=1:sink_gamma=1:vtk=true
python qcdutils_vtk.py -u 0.5 -l 0.005 cold.mdp.point.at00016.00004.00004.00004.s3.c2.current-static.vtk
open cold.mdp.point.at00016.00004.00004.00004.s3.c2.current-static.vtk.html
 */

void usage() {
  mdp << "file file should be called by qcdutils_run.py (googlecode)";
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
    for(int k=0; k<arguments.get("-cool_vtk","n",20); k++) {
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
  forallsites(x) if(x(TIME)==0) {
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
                        default_partitioning<0>,
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
      x.set(t,y(TIME),y(1),y(2));
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

mdp_matrix parse_gamma(string g) {
  if(g=="I") return Gamma1;
  if(g=="0") return Gamma[0];
  if(g=="1") return Gamma[1];
  if(g=="2") return Gamma[2];
  if(g=="3") return Gamma[3];
  if(g=="5") return Gamma5;
  if(g=="05") return Gamma[0]*Gamma5;
  if(g=="15") return Gamma[1]*Gamma5;
  if(g=="25") return Gamma[2]*Gamma5;
  if(g=="35") return Gamma[3]*Gamma5;
  if(g=="01") return Gamma[0]*Gamma[1];
  if(g=="02") return Gamma[0]*Gamma[2];
  if(g=="03") return Gamma[0]*Gamma[3];
  if(g=="12") return Gamma[1]*Gamma[2];
  if(g=="13") return Gamma[1]*Gamma[3];
  if(g=="23") return Gamma[2]*Gamma[3];
  throw string("undefined gamma structure");
}

void choose_action_and_inverter(mdp_args& arguments) {
  string quark_action = 
    arguments.get("-quark","action","clover_fast|clover_slow|clover_sse2");
  string inverter = 
    arguments.get("-quark","alg","bicgstab|minres|bicgstab-vtk|minres-vtk");
  mdp << "using action=" << quark_action << " inverter=" << inverter << endl;
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
}

void make_quark(gauge_field &U, coefficients &gauge, coefficients &quark,
		mdp_args& arguments, string newfilename) {

  float abs_precision = arguments.get("-quark","abs_precision",1e-12);
  float rel_precision = arguments.get("-quark","rel_precision",1e-12);

  choose_action_and_inverter(arguments);

  if(gauge["c_{SW}"]!=0) compute_em_field(U);
  
  int nc = U.nc;
  fermi_field psi(U.lattice(),nc);
  fermi_field phi(U.lattice(),nc);
  mdp_site x(U.lattice());
  mdp_site z(U.lattice());
  vector<mdp_real> pion(U.lattice().size(TIME));
  vector<mdp_real> meson(U.lattice().size(TIME));
  vector<mdp_real> current(U.lattice().size(TIME));
  int NT = U.lattice().size(TIME);
  int L[3];
  L[0]=U.lattice().size(1);
  L[1]=U.lattice().size(2);
  L[2]=U.lattice().size(3);
  mdp_lattice space(3,L,
		    default_partitioning<0>,
		    torus_topology,
		    0,1,false);
  mdp_field<float> Q(space);
  mdp_site y(space);
  string prefix;
  string quarkfilename;
  mdp_real tmp;
  mdp_complex s1,s2;
  mdp_matrix G1,G2,G3;
 
  // this is conditional because we need all S in some cases
  fermi_propagator S;
  mdp_complex open_prop[4][4][10][10][512];
  if(arguments.have("-current-static")||
     arguments.have("-meson")||
     arguments.have("-baryon"))
    S.allocate_fermi_propagator(U.lattice(),U.nc);

  int t0 = arguments.get("-quark","source_t",0);
  int x0 = arguments.get("-quark","source_x",0);
  int y0 = arguments.get("-quark","source_y",0);
  int z0 = arguments.get("-quark","source_z",0);
  if(arguments.get("-quark","source_point","zero|center")=="center") {
    t0 = NT/2; x0 = L[0]/2; y0 = L[1]/2; z0 = L[2]/2;    
  }

  string source_type = arguments.get("-quark","source_type","point|wall");
  for(int a=0; a<4; a++)
    for(int i=0; i<nc; i++) {	
      mdp << "quark source spin="<<a<<" color="<<i<<endl;
      if(source_type=="point") {
	psi = 0;      
	if (on_which_process(U.lattice(),t0,x0,y0,z0)==ME) {
	  x.set(t0,x0,y0,z0);
	  psi(x,a,i)=1;
	}
      } else if(source_type=="wall") {
	forallsites(x) {
	  forspincolor(b,j,U.nc) {
	    psi(x,b,j)=(x(TIME)==t0 && b==a && j==i)?1:0;
	  }
	}
      }
      // optional ... smer source here
      psi.update();
      prefix = newfilename+"."+source_type+".k"+tostring(quark["kappa"]);
      if (t0*t0+x0*x0+y0*y0+z0*z0 > 0)
	prefix = prefix + ".at"+tostring(t0) + "." + tostring(x0) + \
	  "." + tostring(y0) + "." + tostring(z0);
      inversion_vtk_prefix = prefix + ".s"+tostring(a,1)+".c"+tostring(i,1);
      quarkfilename = inversion_vtk_prefix + ".quark";
      if (arguments.get("-quark","load","false|true")=="true") {
	phi.load(quarkfilename);
      } else {
	mul_invQ(phi,psi,U,quark,abs_precision,rel_precision);
	// optional ... insert smearing here
	phi.save(quarkfilename);
      }
      if (arguments.have("-pion")) {
	if (a==0 && i==0) {
	  for(int t=0; t<NT; t++)
	    pion[(t-t0+NT)%NT]=0.0;
	  Q=0;
	}
	forallsitesandcopies(x) {
	  y.set(x(1),x(2),x(3));
	  for(int b=0; b<4; b++)
	    for(int j=0; j<nc; j++) {
	      tmp = real(phi(x,b,j)*conj(phi(x,b,j)));	  
	      pion[(x(TIME)-t0+NT)%NT] += tmp;
	      Q(y) += tmp;
	    }
	}
      }
      if(arguments.have("-4quarks")||arguments.have("-meson")||arguments.have("-current-static"))
	forallsites(x)
	  for(int b=0; b<4; b++)
	    for(int j=0; j<nc; j++)
	      S(x,a,b,i,j) = phi(x,b,j);
    }
  if(arguments.have("-pion")) {  
    mpi.add(&pion[0],NT);	
    pretty_print("C2",pion);      
    if (arguments.get("-pion","vtk","false|true")=="true")
      Q.save_vtk(prefix+".pion.vtk");
  }
  /// mind - before here Q ony to be used for pion
  if(arguments.have("-meson")) {
    Q=0;
    for(int t=0; t<NT; t++) meson[t]=0;
    G1 = Gamma5*parse_gamma(arguments.get("-meson","source_gamma","5"));
    G2 = parse_gamma(arguments.get("-meson","sink_gamma","5"))*Gamma5;
    forspincolor(a,i,U.nc) {
      forspincolor(b,j,U.nc) {
	forallsites(x) {
	  s1=s2=0;
	  for(int c=0;c<4;c++) {
	    s1 += S(x,a,c,i,j)*G2(c,b);
	    s2 += conj(S(x,c,b,i,j))*G1(c,a);
	  }
	  tmp = abs(s1*s2);
	  meson[(x(TIME)-t0+NT)%NT] += tmp;
	}
      }
    }
    mpi.add(&meson[0],NT);	
    pretty_print("C2_meson",meson);      
    if (arguments.get("-meson","vtk","false|true")=="true")
      Q.save_vtk(prefix+".meson.vtk");
  }
  if(arguments.have("-current-static")) {
    /// this part does not work in parallel (yet)
    Q = 0;
    G1 = parse_gamma(arguments.get("-current-static","source_gamma","5"))*Gamma5;
    G2 = parse_gamma(arguments.get("-current-static","sink_gamma","5"));
    G3 = Gamma5*parse_gamma(arguments.get("-current-static","current_gamma","I"));
    mdp_matrix_field Sh(U.lattice(),U.nc,U.nc);
    for(int t=0; t<NT; t++) meson[t]=0;
    for(int t=0; t<U.lattice().size(TIME)/2;t++)
      if(t==0) {	
	forallsites(x)
	  if(x(TIME)==(NT+t+t0)%NT) Sh(x)=1;   /// FIX PROBLEM WITH CENTER
      } else {
	forallsites(x)
	  if(x(TIME)==(NT+t+t0)%NT) {
	    z.set((NT+t0-t)%NT,x(1),x(2),x(3));	  
	    Sh(x)=U(z,0)*Sh(x-0)*U(x-0,0);
	  }
      }			     
    forallsites(x) 
      if(x(TIME)>=0) {
	z.set((NT+2*t0-x(TIME))%NT,x(1),x(2),x(3));
	y.set(x(1),x(2),x(3));
	forspincolor(a,i,U.nc) { 
	  forspincolor(b,j,U.nc) {
	    s1 = s2 = 0;
	    for(int c=0; c<4; c++) {
	      s1 += conj(S(z,c,a,j,i))*G3(c,b);
	      for(int k=0; k<U.nc; k++)
		s2 += S(x,b,c,j,k)*G2(c,a)*conj(Sh(x,k,i));
	    }
	    if(a>2) { // non-zero components of (1-gamma[0])/2
	      tmp = abs(s1*s2);
	      current[(x(TIME)-t0+NT)%NT] += tmp;
	      Q(y) += tmp;
	    }
	  }
	}
      }
    mpi.add(&current[0],NT);	
    pretty_print("C2_current",current);      
    if (arguments.get("-current-static","vtk","false|true")=="true")
      Q.save_vtk(prefix+".current-static.vtk");
  }
  if(arguments.have("-4quarks")) {
    mdp_matrix G = 
      parse_gamma(arguments.get("-4quark","source","5|I|0|1|2|3|05|15|25|35|01|02|03|12|13|23"))*Gamma5;
    forspincolor(a,i,U.nc) {
      forspincolor(b,j,U.nc) {
	for(int t=0; t<U.lattice().size(TIME); t++)
	  open_prop[a][b][i][j][t] = 0.0;	
	for(int c=0; c<4; c++) 
	  for(int d=0; d<4; d++) 
	    if(G(c,d)!=0)
	      forallsites(x)
		for(int k=0; k<U.nc; k++)
		  open_prop[a][b][i][j][(x(TIME)-t0+NT)%NT] +=
		    S(x,a,c,i,k)*conj(S(x,b,d,j,k))*G(c,d);
	mpi.add((mdp_complex*) &open_prop[0],NT*4*4*3*3);	      	
      }
    }
    mdp_matrix G1, G2;
    string op4q = arguments.get("-4quarks","operator","5Ix5I|0Ix0I|1Ix1I|2Ix2I|3Ix3I|05Ix05I|15Ix15I|25Ix25I|35Ix35I|01Ix01I|02Ix02I|03Ix03I|12Ix12I|13Ix13I|23Ix23I|5Tx5T|0Tx0T|1Tx1T|2Tx2T|3Tx3T|05Tx05T|15Tx15T|25Tx25T|35Tx35T|01Tx01T|02Tx02T|03Tx03T|12Tx12T|13Tx13T|23Tx23T");
    bool rotate=false;
    if(op4q[op4q.size()-1]=='T') rotate = true;
    G1 = G2 = Gamma5*parse_gamma(op4q.substr(0,op4q.find("x")-1));
    // others operators may be 0Tx0T,1Tx1T,5Tx5T,etc.
    for(int t1=0;t1<NT/2;t1++)
      for(int t2=0;t2<NT/2;t2++) {
	int t1s = (NT+t0-t1) % NT;
	int t2s = (NT+t0+t2) % NT;
	mdp_real c3a = 0;
	mdp_real c3b = 0;
	// manually add other contractions....
	for(int a=0; a<4; a++)
	  for(int b=0; b<4; b++)
	    for(int c=0; c<4; c++)
	      for(int d=0; d<4; d++) {
		mdp_complex g1 = G1(b,a);
		mdp_complex g2 = G2(d,c);
		if(g1!=0 && g2!=0) 
		  for(int i=0; i<U.nc; i++)
		    for(int j=0; j<U.nc; j++)
		      if(!rotate) {
			c3a+=abs(open_prop[a][b][i][i][t1s]*g1*
				 open_prop[c][d][j][j][t2s]*g2);
			c3b+=abs(open_prop[c][b][j][i][t1s]*g1*
				 open_prop[a][d][i][j][t2s]*g2);
		      } else
			for(int z=1; z<9; z++)
			  for(int k1=0; k1<U.nc; k1++)
			    for(int k2=0; k2<U.nc; k2++) {
			      c3a+=abs(open_prop[a][b][i][k1][t1s]*g1*
				       Lambda[z](k1,i)*
				       open_prop[c][d][j][k2][t2s]*g2*
				       Lambda[z](k2,j))/4;
			      c3b+=abs(open_prop[c][b][j][k1][t1s]*g1*
				       Lambda[z](k1,i)*
				       open_prop[a][d][i][k2][t2s]*g2*
				       Lambda[z](k2,j))/4;
			    }
	      }
	mdp << "C3a[" << t1 << "]["<< t2 << "] = " << c3a << endl;
	mdp << "C3b[" << t1 << "]["<< t2 << "] = " << c3b << endl;
      }
  }
}

int main(int argc, char** argv) {
  mdp.open_wormholes(argc,argv);
  mdp_args arguments(argc,argv);
  if(!arguments.length())
    usage();
  define_base_matrices(arguments.get("-quark","matrices","FERMILAB|MILC|UKQCD|Minkowsy-Dirac|Minkowsy-Chiral"));
  coefficients gauge;
  coefficients quark;
  int ndim = 4;
  int size[4];
  string filename, newfilename, vtkfilename;
  vector<string> filenames;
  int nt = arguments.get("-gauge","nt",16);
  int nx = arguments.get("-gauge","nx",4);
  int ny = arguments.get("-gauge","ny",nx);
  int nz = arguments.get("-gauge","nz",ny);
  int nc = arguments.get("-gauge","nc",3);
  size[0]=nt; size[1]=nx; size[2]=ny; size[3]=nz;
  if (arguments.get("-gauge","start","load|cold|hot|instantons")=="cold") {
    filenames.push_back("cold.mdp"); 
  } else if (arguments.get("-gauge","start","load|cold|hot|instantons")=="hot") {
    filenames.push_back("hot.mdp");
  } else if (arguments.get("-gauge","start","load|cold|hot|instantons")=="instantons") {
    filenames.push_back("custom.mdp");
  } else if(arguments.get("-gauge","start","load|cold|hot|instantons")=="load") {
    string pattern = arguments.get("-gauge","load","demo.mdp");
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
    exit(TIME);
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

  mdp_lattice lattice(ndim,size,
		      default_partitioning<1>,
		      torus_topology,
		      0,1,(arguments.get("-gauge","start","load|cold|hot|instantons")!="load" || nconfigs>0));
  gauge_field U(lattice,nc);
  for(int f=0; f<filenames.size(); f++) {
    filename = filenames[f];
    if (arguments.get("-gauge","start","load|cold|hot|instantons")=="cold") {
      set_cold(U);
      if (arguments.get("-gauge","save","true")=="true")
	U.save(filename);
    } else if (arguments.get("-gauge","start","load|cold|hot|instantons")=="hot") {
      set_hot(U);
      if (arguments.get("-gauge","save","true")=="true")
	U.save(filename);
    } else if (arguments.get("-gauge","start","load|cold|hot|instantons")=="instantons") {
      float t0 = arguments.get("-gauge","t0",0);
      float x0 = arguments.get("-gauge","x0",0);
      float y0 = arguments.get("-gauge","y0",0);
      float z0 = arguments.get("-gauge","z0",0);
      float r0 = arguments.get("-gauge","r0",1.0);
      float t1 = arguments.get("-gauge","t1",1);
      float x1 = arguments.get("-gauge","x1",1);
      float y1 = arguments.get("-gauge","y1",1);
      float z1 = arguments.get("-gauge","z1",1);
      float r1 = arguments.get("-gauge","r1",0.0);
      InstantonGenerator4D generator;      
      vector<SingleInstanton4D> instantons;    
      instantons.push_back(SingleInstanton4D(t0,x0,y0,z0,abs(r0),(r0>0)?+1:-1));
      if(r1!=0) 
	instantons.push_back(SingleInstanton4D(t1,x1,y1,z1,abs(r1),(r1>0)?+1:-1));
      generator.generate(U,instantons);      
    } else if (arguments.get("-gauge","start","load|cold|hot|instantons")=="load") {
      cout << filename << "\n";
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
      if (arguments.have("-cool-vtk")) {
	cool_vtk(U,arguments,newfilename);
      } else if (arguments.have("-cool")) {
	cool(U,arguments);
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

