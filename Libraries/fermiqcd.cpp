#include "fermiqcd.h"
/*
./a.out -gauge:start=cold:nt=32:nx=8 -quark:kappa=0.115:source_point=center:load=true -pion:vtk=true -current_static:source_gamma=1:sink_gamma=1:vtk=true
python qcdutils_vtk.py -u 0.5 -l 0.005 cold.mdp.point.at00016.00004.00004.00004.s3.c2.current_static.vtk
open cold.mdp.point.at00016.00004.00004.00004.s3.c2.current_static.vtk.html
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
			 arguments.get("-cool_vtk","alpha",0.7),
			 arguments.get("-cool_vtk","steps",1),
			 arguments.get("-cool_vtk","cooling",10));
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

void pretty_print(string prefix, vector<mdp_real> data) {
  if (ME==0) {
    for(int t=0; t<data.size(); t++) {
      cout << prefix << "[" << t << "] = " << data[t] << endl;
    }
  }
}

void make_quark(gauge_field &U, coefficients &gauge, coefficients &quark,
		mdp_args& arguments, string newfilename) {

  float abs_precision = arguments.get("-quark","abs_precision",1e-12);
  float rel_precision = arguments.get("-quark","rel_precision",1e-12);
  string quark_action = 
    arguments.get("-quark","action","clover_fast|clover_slow|clover_sse2");
  string inverter = 
    arguments.get("-quark","alg","bicgstab|minres|bicgstabvtk|minresvtk");
  mdp << "using action=" << quark_action << " inverter=" << inverter << endl;

  select_action_and_inverter(quark_action, inverter);

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
  mdp_field<float> Q(U.lattice());
  string prefix;
  string quarkfilename;
  mdp_real tmp;
  mdp_complex s1,s2;
  mdp_matrix G1,G2,G3,G4;
 
  // this is conditional because we need all S in some cases
  fermi_propagator S;
  mdp_complex open_prop[4][4][10][10][512];
  bool use_propagator = 
    arguments.have("-current_static")||
    arguments.have("-meson")||
    arguments.have("-wave_static")||
    arguments.have("-baryon");
  if(use_propagator) S.allocate_fermi_propagator(U.lattice(),U.nc);

  int t0 = arguments.get("-quark","source_t",0);
  int x0 = arguments.get("-quark","source_x",0);
  int y0 = arguments.get("-quark","source_y",0);
  int z0 = arguments.get("-quark","source_z",0);
  if(arguments.get("-quark","source_point","zero|center")=="center") {
    t0 = 0; x0 = L[0]/2; y0 = L[1]/2; z0 = L[2]/2;    
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
      if(arguments.get("-quark","vtk","false|true")=="true") {
	forspincolor(b,j,nc) {
	  forallsites(x) Q(x)=abs(phi(x,b,j));
	  Q.save_vtk(inversion_vtk_prefix+".quark"+tostring(b,1)+tostring(j,1)+".vtk",-1);
	}
      }
      if (arguments.have("-pion")) {
	if (a==0 && i==0) {
	  for(int t=0; t<NT; t++)
	    pion[(t-t0+NT)%NT]=0.0;
	  Q=0;
	}
	forallsitesandcopies(x) {
	  forspincolor(b,j,nc) {
	    tmp = real(phi(x,b,j)*conj(phi(x,b,j)));	  
	    pion[(x(TIME)-t0+NT)%NT] += tmp;
	    Q(x) += tmp;
	  }
	}
      }
      if(use_propagator) {
	forallsites(x) {
	  forspincolor(b,j,nc) {
	    S(x,a,b,i,j) = phi(x,b,j);
	  }
	}
      }
    }
  if(use_propagator) {
    int smear_steps = arguments.get("-quark","smear_steps",0);
    float smear_alpha = arguments.get("-quark","smear_alpha",1.0);
    smear_propagator(S,U,smear_steps,smear_alpha);
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
    G1 = Gamma5*parse_gamma(arguments.get("-meson","source_gamma","5|0|1|2|3|01|02|03|12|13|05|15|25|35|I"));
    G2 = parse_gamma(arguments.get("-meson","sink_gamma","5|0|1|2|3|01|02|03|12|13|05|15|25|35|I"))*Gamma5;
    forspincolor(a,i,U.nc) {
      forspincolor(b,j,U.nc) {
	forallsites(x) {
	  s1=s2=0;
	  for(int c=0;c<4;c++) {
	    s1 += S(x,a,c,i,j)*G2(c,b);
	    s2 += conj(S(x,c,b,i,j))*G1(c,a);
	  }
	  tmp = real(s1*s2);
	  meson[(x(TIME)-t0+NT)%NT] += tmp;
	  Q(x) += tmp;
	}
      }
    }
    mpi.add(&meson[0],NT);	
    pretty_print("C2_meson",meson);      
    if (arguments.get("-meson","vtk","false|true")=="true")
      Q.save_vtk(prefix+".meson.vtk");
  }
  if(arguments.have("-current_static")) {
    /// this part does not work in parallel (yet)
    Q = 0;
    G1 = parse_gamma(arguments.get("-current_static","source_gamma","5|0|1|2|3|01|02|03|12|13|05|15|25|35|I"))*Gamma5;
    G2 = parse_gamma(arguments.get("-current_static","sink_gamma","5|0|1|2|3|01|02|03|12|13|05|15|25|35|I"));
    G3 = Gamma5*parse_gamma(arguments.get("-current_static","current_gamma","I|0|1|2|3|5|01|02|03|12|13|05|15|25|35"));
    G4 = G2*(1-Gamma[0])/2*G1;
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
    forallsites(x) // FIX G1 !!!!!!!!
      if(x(TIME)>=0) {
	z.set((NT+2*t0-x(TIME))%NT,x(1),x(2),x(3));
	forspincolor(a,i,U.nc) { 
	  forspincolor(b,j,U.nc) {
	    s1 = s2 = 0;
	    for(int c=0; c<4; c++) {
	      s1 += conj(S(z,c,a,j,i))*G3(c,b);
	      for(int k=0; k<U.nc; k++)
		s2 += S(x,b,c,j,k)*G4(c,a)*conj(Sh(x,i,k));
	    }
	    tmp = real(s1*s2);
	    current[(x(TIME)-t0+NT)%NT] += tmp;
	    Q(x) += tmp;
	  }
	}
      }
    mpi.add(&current[0],NT);	
    pretty_print("C2_current",current);      
    if (arguments.get("-current_static","vtk","false|true")=="true")
      Q.save_vtk(prefix+".current_static.vtk");
  }
  if(arguments.have("-4quark")) {
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
    string op4q = arguments.get("-4quark","operator","5Ix5I|0Ix0I|1Ix1I|2Ix2I|3Ix3I|05Ix05I|15Ix15I|25Ix25I|35Ix35I|01Ix01I|02Ix02I|03Ix03I|12Ix12I|13Ix13I|23Ix23I|5Tx5T|0Tx0T|1Tx1T|2Tx2T|3Tx3T|05Tx05T|15Tx15T|25Tx25T|35Tx35T|01Tx01T|02Tx02T|03Tx03T|12Tx12T|13Tx13T|23Tx23T");
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
			c3a+=real(open_prop[a][b][i][i][t1s]*g1*
				  open_prop[c][d][j][j][t2s]*g2);
			c3b+=real(open_prop[c][b][j][i][t1s]*g1*
				  open_prop[a][d][i][j][t2s]*g2);
		      } else
			for(int z=1; z<9; z++)
			  for(int k1=0; k1<U.nc; k1++)
			    for(int k2=0; k2<U.nc; k2++) {
			      c3a+=real(open_prop[a][b][i][k1][t1s]*g1*
					Lambda[z](k1,i)*
					open_prop[c][d][j][k2][t2s]*g2*
					Lambda[z](k2,j))/4;
			      c3b+=real(open_prop[c][b][j][k1][t1s]*g1*
					Lambda[z](k1,i)*
					open_prop[a][d][i][k2][t2s]*g2*
					Lambda[z](k2,j))/4;
			    }
	      }
	mdp << "C3a[" << t1 << "]["<< t2 << "] = " << c3a << endl;
	mdp << "C3b[" << t1 << "]["<< t2 << "] = " << c3b << endl;
      }
  }
  if(arguments.have("-wave_static")) {
    throw string("NotImplemented");
    // WORK IN PROGRESS - only works on cold - no paths
    string source_gamma = arguments.get("-wave_static","source_gamma","5|0|1|2|3|01|02|03|12|13|05|15|25|35|I");
    int smear_steps = arguments.get("-wave_static","smear_steps",10);
    G1 = parse_gamma(source_gamma)*(1+Gamma[0])/2;
    Q=0;
    forallsites(x)
      for(int a=0; a<4; a++)
	for(int b=0; b<4; b++)
	  if(G1(a,b)!=0)
	    for(int i=0; i<U.nc; i++)
	      Q(x) += pow(abs(S(x,b,a,i,i)*G1(a,b)),2);
    // smear_propagator(S,U,smear_steps);
    Q.save_vtk(prefix+string(".")+source_gamma+".wave.vtk");
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
  string gauge_start = arguments.get("-gauge","start","load|cold|hot|instantons");
  if (gauge_start=="cold")
    filenames.push_back("cold.mdp"); 
  else if (gauge_start=="hot")
    filenames.push_back("hot.mdp");
  else if (gauge_start=="instantons")
    filenames.push_back("custom.mdp");
  else if(gauge_start=="load") {
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
		      0,1,(gauge_start!="load" || nconfigs>0));
  gauge_field U(lattice,nc);
  for(int f=0; f<filenames.size(); f++) {
    filename = filenames[f];
    if (gauge_start=="cold") {
      set_cold(U);
      if (arguments.get("-gauge","save","true")=="true")
	U.save(filename);
    } else if (gauge_start=="hot") {
      set_hot(U);
      if (arguments.get("-gauge","save","true")=="true")
	U.save(filename);
    } else if (gauge_start=="instantons") {
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
    } else if (gauge_start=="load") {
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
      if (arguments.have("-cool_vtk")) {
	cool_vtk(U,arguments,newfilename);
      } else if (arguments.have("-cool")) {
	cool(U,arguments);
      }
      if (arguments.have("-plaquette_vtk")) {
	plaquette_vtk(U,newfilename+".plaquette.vtk");
      }      
      if (arguments.have("-polyakov_vtk")) {
	polyakov_vtk(U,newfilename+".polyakov.vtk");
      }      
      if (arguments.have("-topcharge_vtk")) {
	float tc = topological_charge_vtk(U,newfilename+".topcharge.vtk",-1);
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

