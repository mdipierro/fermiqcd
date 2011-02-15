#include "fermiqcd.h"



class FermiQCD {
public:
  int nc;
  int L[4];
  mdp_lattice* plattice;
  gauge_field* pU;
  string prefix;
  int counter;
  ofstream os;
  string make_prefix() {
    time_t time0;
    time(&time0);
    string s=ctime(&time0);
    while(1) {
      int i=s.find(" ");
      if(i>=0) s=s.replace(i,1,"_");
      else return s.substr(0,s.size()-1)+"/";
    }
    return s.substr(0,s.size()-1)+"/";
  }  
  FermiQCD() {
    plattice=0;
    pU=0;
    prefix=make_prefix();    
    if(ME==0) system((string("mkdir ")+prefix).c_str());    
    if(ME==0) os.open((prefix+"README.log").c_str());
    os << "prefix=" << prefix << endl;
    os << "initialization completed\n";
    os << "running on npocs=" << Nproc << " parallel processes\n";      
    os << endl;
    counter=0;
  }
  template<class T>
  FermiQCD& operator<<(const T& obj) {
    cout << ME << ":" << obj;
    if(ME==0) {
      os << obj;
      os.flush();
    }
    mdp << obj;
    return *this;
  }
  ~FermiQCD() {
    if(!plattice) delete plattice;
    if(!pU) delete pU;
    os.close();
  }
  void init_cold(int LT,int LX, int LY,int LZ, int nc) {
    this->L[0]=LT;
    this->L[1]=LX;
    this->L[2]=LY;
    this->L[3]=LZ;
    this->nc=nc;
    os << "making an lattice TxXxYxZ=" << LT << "x" << LX << "x" << LY << "x" << LZ << endl;
    this->plattice=new mdp_lattice(4,L);    
    os << "making a cold gauge configuration U with nc=" << nc << endl; 
    this->pU=new gauge_field(*plattice,nc);    
    cout << ME << endl;
    pU->update();
    // set_cold(*pU);  
    //prefix=prefix+string("C")+tostring((int)mdp.time());
    counter=0;
  }
  void init_hot(int LT,int LX, int LY,int LZ, int nc) {
    this->L[0]=LT;
    this->L[1]=LX;
    this->L[2]=LY;
    this->L[3]=LZ;
    this->nc=nc;
    os << "making an lattice TxXxYxZ=" << LT << "x" << LX << "x" << LY << "x" << LZ << endl;
    this->plattice=new mdp_lattice(4,L);
    os << "making a hot gauge configuration U with nc=" << nc << endl; 
    this->pU=new gauge_field(*plattice,nc);
    set_hot(*pU);
    prefix=string("C")+tostring((int)mdp.time());
    counter=0;
  }
  void init_load(string filename) {
    mdp_field_file_header header=get_info(filename);
    this->L[0]=header.box[0];
    this->L[1]=header.box[1];
    this->L[2]=header.box[2];
    this->L[3]=header.box[3];
    // (4-8)*4*(1-4-9-25-36-49-64-81-100)
    int precision=4;
    switch(header.bytes_per_site) {
      case 4*4*1: precision=4; this->nc=1; break; 
      case 8*4*1:  precision=8; this->nc=1; break; 
      case 4*4*4:  precision=4; this->nc=2; break; 
      case 8*4*4:  precision=8; this->nc=2; break; 
      case 4*4*9:  precision=4; this->nc=3; break; 
      case 8*4*9:  precision=8; this->nc=3; break; 
      case 4*4*16:  precision=4; this->nc=4; break; 
      case 8*4*16:  precision=8; this->nc=4; break; 
      case 4*4*25:  precision=4; this->nc=5; break; 
      case 8*4*25:  precision=8; this->nc=5; break; 
      case 4*4*36:  precision=4; this->nc=6; break; 
      case 8*4*36:  precision=8; this->nc=6; break; 
      case 4*4*49:  precision=4; this->nc=7; break; 
      case 8*4*49:  precision=8; this->nc=7; break; 
      case 4*4*64:  precision=4; this->nc=8; break; 
      case 8*4*64:  precision=8; this->nc=8; break; 
      case 4*4*81:  precision=4; this->nc=9; break; 
      case 8*4*81:  precision=8; this->nc=9; break; 
      case 4*4*100:  precision=4; this->nc=10; break; 
      case 8*4*100:  precision=8; this->nc=10; break;
    }    
    os << "making an lattice TxXxYxZ=" << this->L[0] << "x" << this->L[1] << "x" << this->L[2] << "x" << this->L[3] << endl;
    this->plattice=new mdp_lattice(4,L);
    os << "loading gauge configuration " << filename << " into U with nc=" << nc << endl; 
    this->pU=new gauge_field(*plattice,nc);
    if(precision==4) pU->load_as_float(filename);    
    if(precision==8) pU->load_as_double(filename);    
    prefix=filename;
    counter=0;
  }
  void wilson_heatbath(float beta, int steps=1) {
    coefficients gauge;
    gauge["beta"]=beta;
    for(int step=0; step<steps; step++) {
      (*this) << "WilsonGaugeAction::heatbath(beta=" << beta << ") step=" << step << endl;
      WilsonGaugeAction::heatbath(*pU,gauge,1);      
    }
    (*this) << "average_plaquette=" << average_plaquette(*pU);    
    counter=counter+steps;
  }
  void ape_smearing(float alpha=0.7, int iterations=20, int cooling_steps=10) {    
    (*this) << "starting Ape smering...\n";
    ApeSmearing::smear(*pU,alpha,iterations,cooling_steps);
    (*this) << "done!\n";
  }
  void coulomb_gauge_fix(int iterations=100) {
    (*this) << "starting coulomb gauge fixing...\n";
    float precision=1e-6;
    float boost=1;
    gaugefixing_stats stats=GaugeFixing::fix(*pU,GaugeFixing::Coulomb,iterations,precision,boost,false);
    (*this) << "done!\n";  
  }
  void coulomb_z3_gauge_fix(int iterations=100) {
    (*this) << "starting coulomb gauge fixing...\n";
    float precision=1e-6;
    float boost=1;
    gaugefixing_stats stats=GaugeFixing::fix(*pU,GaugeFixing::Coulomb,iterations,precision,boost,true);
    (*this) << "done!\n";  
  }
  void landau_gauge_fix(int iterations=100) {
    (*this) << "starting coulomb gauge fixing...\n";
    float precision=1e-6;
    float boost=1;
    gaugefixing_stats stats=GaugeFixing::fix(*pU,GaugeFixing::Landau,iterations,precision,boost,false);
    (*this) << "done!\n";  
  }
  void compute_topological_charge() {
    mdp_field<float> Q(*plattice);
    topological_charge(Q,*pU);
    Q.save_vtk(prefix+"topological_charge");
  }  
  void compute_partitioning() {
    mdp_site x(*plattice);
    mdp_field<float> Q(*plattice);
    forallsites(x)
      p(x)=on_which_process(*plattice,x(0),x(1),x(2),x(3));
    Q.save_vtk(prefix+"partitioning");
  }
};

int main(int argc, char** argv) {
  mdp.open_wormholes(argc,argv);

  
  FermiQCD fermiqcd;
  fermiqcd.init_cold(10,10,4,4,3);
  fermiqcd.compute_partitioning();
  fermiqcd.compute_topological_charge();
 

  mdp.close_wormholes();
  return 0;
}
