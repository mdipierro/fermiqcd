/////////////////////////////////////////////////////////////////
/// @file mdp_site.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains delcaration of class mdp_site for complex numbers
///
/// Licensed under GPL2 license
/// Read attached license in file mdp_license.pdf
/// This file cannot be distributed without file mdp_license.pdf
//////////////////////////////////////////////////////////////////

/// checks which process of the lattice a stores locally the site of 
/// coordinates x0,x1,x2,...,x9
/// to be used before calling mdp_site::set()
/// (note: prototyping of friend functions is required by some compilers)
int on_which_process(mdp_lattice &a, 
                     int x0=0, int x1=0, int x2=0, int x3=0, int x4=0, 
                     int x5=0, int x6=0, int x7=0, int x8=0, int x9=0);


/// @brief site object to loop on a lattice
///
/// Example:
/// @verbatim
///   int box[]={10,10,10};
///   mdp_lattice lattice(3,box);
///   mdp_site x(lattice);
///   forallsites(x) cout << x << endl;
///   if(on_which_process(lattice,1,1,1)==ME) {
///      x.set(1,1,1);
///      cout << lattice.random(x).plain() << endl;
///   }
/// @endverbatim
class mdp_site {
 private:
  mdp_lattice *ptr;  /// this points to the lattice for this field 
 public:
  mdp_int idx;          /// value of the mdp_site in local coordinate     
#ifdef BLOCKSITE
  int block[BLOCKSITE];
#endif
  mdp_site() {
    // one whould not use this!
    ptr=0;
    idx=0;
  }
  /// declares object of class mdp_site living on the 
  /// lattice passed by reference
  mdp_site(const mdp_lattice &a) {
    on(a);
  }
  void on(const mdp_lattice &a) {
    ptr=(mdp_lattice*) &a;
    idx=(*ptr).start[ME][0];
#ifdef BLOCKSITE
    for(int k=0; k<BLOCKSITE; k++) block[k]=0;
#endif
  }
  /// returns by reference the lattice the site lives on
  inline mdp_lattice &lattice() {
    return *ptr;
  }
  inline mdp_site(mdp_int i, mdp_lattice* ptr2) {
    idx=i;
    ptr=ptr2;
#ifdef BLOCKSITE
    for(int k=0; k<BLOCKSITE; k++) block[k]=0;
#endif
  }
#ifdef BLOCKSITE
  inline mdp_site(mdp_int i, mdp_lattice* ptr2, int b[], int sign=0, int mu=0) {
    idx=i;
    ptr=ptr2;
    for(int k=0; k<BLOCKSITE; k++) block[k]=b[k];
    block[mu]+=sign;
  }
#endif
  mdp_site(const mdp_site &x) {
    idx=x.idx;
    ptr=x.ptr;
#ifdef BLOCKSITE
    for(int k=0; k<BLOCKSITE; k++) block[k]=x.block[k];
#endif
  }
  inline mdp_site operator=(mdp_int i) {
    idx=lattice().start[ME][0]+i;
    return mdp_site(idx, ptr);
  }
  inline mdp_site operator=(mdp_site x) {
    if(ptr==x.ptr) idx=x.idx;
    else           set_global(x.global_index());
#ifdef BLOCKSITE
    for(int k=0; k<BLOCKSITE; k++) block[k]=x.block[k];
    return mdp_site(idx, ptr, block);
#endif
    return mdp_site(idx, ptr);
  }
  inline int operator==(mdp_site x) {
    if((idx==NOWHERE) || (x.idx==NOWHERE)) return false;
    return (global_index()==x.global_index());
  }
  inline int operator!=(mdp_site x) {
    if((idx==NOWHERE) || (x.idx==NOWHERE)) return true;
    return (global_index()!=x.global_index());
  }
  inline void start(int np=0) {
#ifdef BLOCKSITE
  for(int k=0; k<BLOCKSITE; k++) block[k]=0;
#endif
    idx=lattice().start[ME][np];
  }
  inline void next() {
    idx++;
  }
  /// checks if the site is inside the portion of the lattice stored by
  /// the current process
  inline int is_in() {
    if((idx>=lattice().start[ME][0]) && 
       (idx<lattice().stop[ME][1])) return 1;
    else return 0;
  }
  /// checks if the site is inside the portion of the lattice stored by
  /// the current process or if the site is in a local copy of a remote
  /// site 
  inline int is_here() {
    if((idx>=0) && (idx<lattice().nvol)) return 1;
    else return 0;
  }
  /// returns the parity EVEN or ODD of the site
  inline int parity() {
    return lattice().parity[idx];
  }
  /// true if the site is stored locally as a copy of a site local 
  /// in another process
  inline int is_in_boundary() {
    return (lattice().wh[idx]!=ME);
  }
  /// returns the local index of the site
  /// local index is assigned by the process to the local sites and copies
  /// of remote sites. local index is not unique thoughout the lattice.
  inline mdp_int local_index() {
    return idx;
  }
  /// returns the global (unique) index of the site
  inline mdp_int global_index() {
    return lattice().gl[idx];
  }
  /// sets the site by its local index (dangerous)
  inline void set_local(mdp_int idx2) {
    idx=idx2;
  }
  /// sets the site by its global index
  inline void set_global(mdp_int idx_gl) {
    mdp_int idx2=lattice().local(idx_gl);
    if(idx2==NOWHERE) 
      error("set_global() trying to access a site that is not here");
    idx=idx2;
  }
  /// returns the site shifted forward in direction mu=(0...ndim-1)
  inline mdp_site operator+ (int mu) {
    mdp_int idx2=lattice().up[idx][mu];
    if(idx2==NOWHERE) {
      cout << ME << " " << (*this)(0) << (*this)(1) 
	   << (*this)(2) << (*this)(3) << " " << mu << endl; 
      error("You cannot exit from your portion of the lattice");
    }
#ifdef BLOCKSITE
    if((mu<BLOCKSITE) && (lattice().co[idx][mu]==lattice().nx[mu]-1))
      return mdp_site(idx2, ptr, block,1,mu);
    return mdp_site(idx2, ptr, block);
#endif
    return mdp_site(idx2, ptr);
  }
  /// returns the site shifted backwards in direction mu=(0...ndim-1)
  inline mdp_site operator- (int mu) {
    mdp_int idx2=lattice().dw[idx][mu];
    if(idx2==NOWHERE)
      error("You cannot exit from your portion of the lattice");
#ifdef BLOCKSITE
    if((mu<BLOCKSITE) && (lattice().co[idx][mu]==0))
      return mdp_site(idx2, ptr, block, -1,mu);
    return mdp_site(idx2, ptr, block);
#endif
    return mdp_site(idx2, ptr);
  }
  /// returns a site shifted i position (backwards if i<0 or forward if i>0)
  /// in direction mu=(0...mdim-1) 
  inline mdp_site hop(int i, int mu) {
    mdp_site y(lattice());
    y=(*this);
    while(i!=0) {
      if(i<0) {y=y-mu; i++;}
      if(i>0) {y=y+mu; i--;}
    }
    return y;
  }
  /// sets the site to the coordinates stored in vector v
  inline mdp_site operator= (mdp_vector v) {
    set(v.x[0],v.x[1],v.x[2],v.x[3],v.x[4],
	v.x[5],v.x[6],v.x[7],v.x[8],v.x[9]);
    return *this;
  }
  /// retruns a site similar to the present but
  /// each coordinates mu of the site shifted according to v[mu]
  inline mdp_site operator+ (mdp_vector v) {
    int mu,step;
    mdp_site y=*this;
    for(mu=0; mu<lattice().ndim; mu++) {
      if(v.x[mu]>0) for(step=0; step<v.x[mu]; step++) y=y+mu;
      else          for(step=0; step<-v.x[mu]; step++) y=y-mu;
    }
    return y;
  }
  /// retruns a site similar to the present but
  /// each coordinates mu of the site shifted according to -v[mu]
  inline mdp_site operator- (mdp_vector v) {
    int mu,step;
    mdp_site y=*this;
    for(mu=0; mu<lattice().ndim; mu++) {
      if(v.x[mu]>0) for(step=0; step<v.x[mu]; step++) y=y-mu;
      else          for(step=0; step<-v.x[mu]; step++) y=y+mu;
    }
    return y;
  }
  /// returns mu coordinate of the site
  inline int operator() (int mu) {
    return lattice().co[idx][mu];
  }
  void operator= (int *x) {
    int ndim=lattice().ndim;
    idx=x[0];
    for(int mu=1; mu<ndim; mu++) idx=idx*lattice().nx[mu]+x[mu];
    idx=lattice().local(idx);
    if(idx==NOWHERE) {
      bool print=mdp.print;
      mdp.print=true;
      mdp << "Warning message from ME=" << ME << ":\n";
      mdp << "You assigned a site that is not here!\n";
      mdp.print=print;
    }
  }
  /// sets the site to a the location spacified by the coordinates
  /// and assumes the site is local (or at least a copy)
  /// @see on_which_process()
  void set(int x0, int x1=0, int x2=0, int x3=0, int x4=0, 
	   int x5=0, int x6=0, int x7=0, int x8=0, int x9=0) {
    int ndim=lattice().ndim;
    idx=x0;
    if(ndim>1) idx=idx*lattice().nx[1]+x1;
    if(ndim>2) idx=idx*lattice().nx[2]+x2;
    if(ndim>3) idx=idx*lattice().nx[3]+x3;
    if(ndim>4) idx=idx*lattice().nx[4]+x4;
    if(ndim>5) idx=idx*lattice().nx[5]+x5;
    if(ndim>6) idx=idx*lattice().nx[6]+x6;
    if(ndim>7) idx=idx*lattice().nx[7]+x7;
    if(ndim>8) idx=idx*lattice().nx[8]+x8;
    if(ndim>9) idx=idx*lattice().nx[9]+x9;
    idx=lattice().local(idx);
    if(idx==NOWHERE) {
      bool print=mdp.print;
      mdp.print=true;
      mdp << "Warning message from ME=" << ME << ":\n";
      mdp << "You assigned a site that is not here!\n";
      mdp.print=print;
    }
  }
  int operator==(int *x) {
    int ndim=lattice().ndim;
    int is_it=1;
    for(int mu=0; mu<ndim; mu++) 
      if(x[mu]!=lattice().co[idx][mu]) is_it=0;
    return is_it;
  }
  int operator!=(int *x) {
    return !(*this == x);
  }
  /// checks the site coordinates vs the coordinates passed as args
  int is_equal(int x0, int x1=0, int x2=0, int x3=0, int x4=0, 
	       int x5=0, int x6=0, int x7=0, int x8=0, int x9=0) {
    int ndim=lattice().ndim;
    int is_it=1;
    if((ndim>0) && (x0!=lattice().co[idx][0])) is_it=0;
    if((ndim>1) && (x1!=lattice().co[idx][1])) is_it=0;
    if((ndim>2) && (x2!=lattice().co[idx][2])) is_it=0;
    if((ndim>3) && (x3!=lattice().co[idx][3])) is_it=0;
    if((ndim>4) && (x4!=lattice().co[idx][4])) is_it=0;
    if((ndim>5) && (x5!=lattice().co[idx][5])) is_it=0;
    if((ndim>6) && (x6!=lattice().co[idx][6])) is_it=0;
    if((ndim>7) && (x7!=lattice().co[idx][7])) is_it=0;
    if((ndim>8) && (x8!=lattice().co[idx][8])) is_it=0;
    if((ndim>9) && (x9!=lattice().co[idx][9])) is_it=0;
    return is_it;
  }
  /// converts a site into a binary number 
  /// to be used only if the site is a vertex of an hypercube centered 
  /// at the origin. this is used to make staggered mesons
  inline friend mdp_int site2binary(mdp_site x) {
    int mu, a=0;
    for(mu=0; mu<x.lattice().ndim; mu++) {
#ifdef CHECK_ALL
      if(fabs(0.5-x(mu))>1) error("site2binary");
#endif
      a+=(0x1 << mu)*x(mu);
    }
    return a;
  }
  /// checks which process of the lattice a stores locally the site of 
  /// coordinates x0,x1,x2,...,x9
  /// to be used before calling mdp_site::set()
  friend int on_which_process(mdp_lattice &a, int x0, int x1, int x2, int x3, 
                              int x4, int x5, int x6, int x7, int x8, int x9) {
    int x[10];
    x[0]=x0;
    if(a.ndim>1) x[1]=x1;
    if(a.ndim>2) x[2]=x2;
    if(a.ndim>3) x[3]=x3;
    if(a.ndim>4) x[4]=x4;
    if(a.ndim>5) x[5]=x5;
    if(a.ndim>6) x[6]=x6;
    if(a.ndim>7) x[7]=x7;
    if(a.ndim>8) x[8]=x8;
    if(a.ndim>9) x[9]=x9;
    return (*(a.where))(x, a.ndim, a.nx);
  }
};

#ifdef MDP_LATTICE

/// Returns the local object mdp_prng at site x of the lattice 
inline mdp_prng &mdp_lattice::random(mdp_site x) {
  if(local_random_generator) {
    if(!x.is_in()) error("request the random generator of a non local site");
    return random_obj[x.idx-start[ME][0]];
  }
  return mdp_random;
}

#endif

/// When compiled with TWISTED_BOUNDARY the mdp_site class keeps track of 
/// sites that moved around the boundary of the torus topology. this function 
/// Returns false if this is one such site, true otherwise.
inline int in_block(mdp_site x) {
#ifdef TWISTED_BOUNDARY
  static int mu;
  for(mu=0; mu<x.lattice().ndim; mu++)
    if(x.block[mu]!=0) return false;
#endif
  return true;
};

ostream& operator<<(ostream& os, mdp_site &x) {
  for(int i=0; i<x.lattice().ndim; i++) {
    if(i==0) os << "(" << x(i);
    else os << "," << x(i);
  }
  os << ")";
  return os;
}
