/////////////////////////////////////////////////////////////////
/// @file mdp_jackboot.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains class mdp_jackboot
///
/// Licensed under GPL2 license
/// Read attached license in file mdp_license.pdf
/// This file cannot be distributed without file mdp_license.pdf
//////////////////////////////////////////////////////////////////

float mdp_jackboot_plain(float *x, void *a);

/// @brief coniatiner class for jackknife and boostrap analysis
///
/// Example:
/// @verbatim
///    mdp_jackboot jb(10,2);
///    for(int k=0; k<10; k++) {       
///       jb(k,0)=mdp_random.plain();
///       jb(k,1)=mdp_random.plain();
///    }
///    jb.plain(0);
///    cout << "mean of jb(k,0) =" << mean() << endl;
///    cout << "jackknife of mean jb(k,0) =" << j_err() << endl;
///    cout << "boostrap of mean jb(k,0) =" << b_err() << endl;
///    jb.plain(1);
///    cout << "mean of jb(k,1) =" << mean() << endl;
///    cout << "jackknife of mean jb(k,1) =" << j_err() << endl;
///    cout << "boostrap of mean jb(k,1) =" << b_err() << endl;
/// @endverbatim
class mdp_jackboot {
 private:
  friend float mdp_jackboot_plain(float *x, void *a) {
    int *i_ptr=(int*) a;
    return x[*i_ptr];
  }
  int mdp_jackboot_plain_int;
 public:
  float (*f)(float*, void *);
  int nconf,narg,conf;
  float *m;
  void* handle;
  mdp_jackboot() {
    m=0;
    dimension(0,0);
  }
  /// allocate container for nconf_ datasets of nargs_ numbers each
  mdp_jackboot(int nconf_, int narg_=1) {
    m=0;
    dimension(nconf_, narg_);
  }
  void dimension(int nconf_, int narg_=1) {
    nconf=nconf_;
    narg=narg_;
    conf=0;
    handle=0;
    if(m!=0) {
      delete[] m;
      m=0;
    }
    if(m==0) m=new float[nconf*narg];
    int i,j;
    for(i=0; i<nconf; i++)
      for(j=0; j<narg; j++)
	m[i*narg+j]=0;
  }
  virtual ~mdp_jackboot() {
    delete[] m;
  }
  float* address(int conf) {
    return m+conf*narg;
  }
  float &operator() (int present_conf, int arg) {
    conf=present_conf;
    return m[conf*narg+arg];
  }
  float &operator() (int arg) {
    return m[conf*narg+arg];
  }
  void plain(int i) {
    if((i<0) || (i>=narg)) error("incorrect mdp_jackboot.plain() argument");
    mdp_jackboot_plain_int=i;
    handle=(void*) &mdp_jackboot_plain_int;
    f=mdp_jackboot_plain;
  }
  void makesample(int *p, int nboot) {
    int boot, j;
    for(boot=0; boot<nboot; boot++)
      for(j=0; j<=conf; j++)
	p[j+(conf+1)*boot]=(int) ((conf+1)*mdp_random.plain());
  }

  float mean() {
    int i,j;
    float tmp;
    float *average= new float[narg];
    if(conf==0) return (*f)(&m[0],handle);
    for(i=0; i<narg; i++) {
      average[i]=0;
      for(j=0; j<=conf; j++)
	average[i]+=m[j*narg+i]/(conf+1);
    }
    tmp=(*f)(average,handle);
    delete[] average;
    return tmp;
  }
  float j_err() {
    int i,j,k;
    float mean0=this->mean();
    float *average= new float[(conf+1)*narg];
    float stddev=0;
    if(conf<2) return 0;
    for(i=0; i<narg; i++)
      for(j=0; j<=conf; j++) {
	average[j*narg+i]=0;
	for(k=0; k<=conf; k++)
	  // typo noticed by Chris Schroeder was fixed!
	  if(j!=k) average[j*narg+i]+=(m[k*narg+i]/(nconf-1));
      }
    stddev=0;
    for(j=0; j<=conf; j++) 
      stddev+=(float) pow((double) (*f)(&(average[j*narg]),handle)-mean0,2.0);
    delete[] average;
    return (float) sqrt(stddev*conf/(conf+1));
  }  
  float b_err(int nboot=100) {
    int i,j,boot,x,y;
    float vx, vy, tmp;
    float *average=new float[nboot*narg];
    int *p=new int[(conf+1)*nboot];
    makesample(p,nboot);
    if(conf==0) return 0;
    for(i=0; i<narg; i++) 
      for(boot=0; boot<nboot; boot++) {
	average[boot*narg+i]=0;
	for(j=0; j<=conf; j++) 
	  average[boot*narg+i]+=m[p[j+(conf+1)*boot]*narg+i];
	average[boot*narg+i]/=conf+1;
      }
    for(x=1; x<nboot; x++) {
      vx=(*f)(&(average[x*narg]),handle);
      for(y=x-1; y>=0; y--) {
	vy=(*f)(&(average[y*narg]),handle);
	if(vy>vx) 
	  for(i=0; i<narg; i++) {
	    tmp=average[(y+1)*narg+i];
	    average[(y+1)*narg+i]=average[y*narg+i];
	    average[y*narg+i]=tmp;
	  } else y=-1;
      }
    }
    vx=(*f)(&(average[((int)((float) nboot*0.16))*narg]),handle);
    vy=(*f)(&(average[((int)((float) nboot*0.84))*narg]),handle);
    delete[] average;
    delete[] p;
    return (float) (vy-vx)/2.0; 
  }
 };
