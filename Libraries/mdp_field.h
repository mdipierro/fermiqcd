/////////////////////////////////////////////////////////////////
/// @file mdp_field.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains declaration of class mdp_field
///
/// Licensed under GPL2 license
/// Read attached license in file mdp_license.pdf
/// This file cannot be distributed without file mdp_license.pdf
//////////////////////////////////////////////////////////////////

/// @brief header for field file IO
///
/// Used to store the binary file haeader that precedes the data
/// When storing an object of class mdp_field<> in a file
class mdp_field_file_header {
 public:
  char  file_id[60];
  char  program_version[60];
  char  creation_date[60];
  uint32_t endianess;
  int32_t   ndim;
  int32_t   box[10];
  int32_t  bytes_per_site;
  int32_t  sites;
  mdp_field_file_header() {
    reset();
  }
  void reset() {
    strcpy(file_id, "File Type: MDP FIELD\n");
    strcpy(program_version,mdp_program_name);
    endianess=mdp_local_endianess;
  }
  void set_time() {
    int i;
    time_t time_and_date;
    time(&time_and_date);
    strcpy(creation_date, ctime(&time_and_date));
    for(i = strlen(creation_date)+1; i<sizeof(creation_date); i++)
      creation_date[i] = '\0';
  }
  /// tries to swicth the endianess of the numerical members of the header
  friend bool switch_header_endianess(mdp_field_file_header &header) {
    if(header.endianess==mdp_local_endianess) return false;
    switch_endianess_byte4(header.endianess);
    if(header.endianess==mdp_local_endianess) {
      switch_endianess_byte4(header.endianess);
      switch_endianess_byte4(header.ndim);
      switch_endianess_byte4(header.box[0]);
      switch_endianess_byte4(header.box[1]);
      switch_endianess_byte4(header.box[2]);
      switch_endianess_byte4(header.box[3]);
      switch_endianess_byte4(header.box[4]);
      switch_endianess_byte4(header.box[5]);
      switch_endianess_byte4(header.box[6]);
      switch_endianess_byte4(header.box[7]);
      switch_endianess_byte4(header.box[8]);
      switch_endianess_byte4(header.box[9]);
      switch_endianess_byte4(header.bytes_per_site);
      switch_endianess_byte4(header.sites);
      return true;
    } else {
      switch_endianess_byte4(header.endianess);
      return false;
    }
  }
};

/// @brief most generic field object
/// 
/// Example:
/// @verbatim
///    int box[]={10,10,10};
///    mdp_lattice lattice(3,box);
///    mdp_field<float> psi(lattice,10);
///    mdp_site x(lattice);
///    forallsites(x)
///      for(int i=0; i<10; i++)
///         psi(x,i)=0.0;
///    psi.update(); // synchronization
///    psi.save("myfield");
///    psi.load("myfield");
/// @endverbatim
template<class T>
class mdp_field  {
 protected:
  mdp_lattice* ptr; /* this points to the lattice for this field  */
  T* m;                 /* this is to store the main field            */
  mdp_int Tsize;
  mdp_int size;            /* this is the size of the field in sizeof(T) */
  int  field_components; /* this is the size of the structure per site */
 public:
  /// the field file header, contains data only if field was read from file
  mdp_field_file_header header;
  /// declare empty field (zero size)
  mdp_field() {
    m=0;
    Tsize=sizeof(T);
    size=field_components=0;
  }
  /// declares a field on lattice a and allocates a vector of n T at each site
  mdp_field(mdp_lattice &a, int n=1) {
    m=0;
    allocate_field(a, n);
  }
  mdp_field(const mdp_field &field) {
    m=0;
    allocate_field(field.lattice(), field.field_components);
    mdp_int i_min=physical_local_start(EVENODD);
    mdp_int i_max=physical_local_stop(EVENODD);
    for(mdp_int i=i_min; i<i_max; i++) m[i]=field.m[i];
  }
  /// checks if a field is allocated of has zero-size
  bool allocated() {
    if(m==0) return false;
    return true;
  }  
  /// Allows dynamical allocation of a field that is not allocated
  void allocate_field(mdp_lattice &a, int n=0) {   
    deallocate_field();
    if(n==0) n=field_components;
    else     field_components=n;
    if(field_components==0) 
      error("You cannot have a field of zero size!");
    size=a.nvol*field_components;
    Tsize=sizeof(T);
    m=new T[size];
    if(m==0) error("OUT OF MEMORY !!!!");
    ptr=&a;
    fill_header();
  }
  /*
  void allocate_field(T* memory_address, mdp_lattice &a, int n=0) {
    m=0;
    if(n==0) n=field_components;
    else     field_components=n;
    if(field_components==0) 
      error("You cannot have a field of zero size!");
    Tsize=sizeof(T);
    size=a.nvol*field_components;
    m=memory_address;
    ptr=&a;
    fill_header();
  }
  */
  void fill_header() {
    int i;
    header.bytes_per_site=Tsize*field_components;
    header.sites=lattice().size();
    header.ndim=lattice().ndim;
    for(i=0; i<lattice().ndim; i++) header.box[i]=lattice().size(i);
    for(; i<10; i++) header.box[i]=0;    
  }
  void deallocate_memory() {
    if(m!=0) delete[] m;
    m=0;
    size=field_components=0;
  }
  /// do not use, may cause memory leaks
  void reset_field() {
    m=0;
    size=field_components=0;
  }
  /// dynamically deallocate field
  void deallocate_field() {
    deallocate_memory();
  }
  virtual ~mdp_field() {
    deallocate_memory();
  }
  /// returns component i of the vector of objects T stored at site x
  inline T& operator() (mdp_site x, int i=0) {
    #ifdef CHECK_ALL
    if(!(x.is_here())) {
      error("You are looking for a site that is not here");
    }
    #endif
    return m[x.idx*field_components+i];
  }

  inline T& operator() (int idx, int i=0) {
    return m[idx*field_components+i];
  }
  /// retruns the address of the vector of objects T stored at site x
  inline T* operator[] (mdp_site x) {
    return address(x,0);
  }
  inline T& operator[] (mdp_int i) {
    return m[i];
  }

  inline T* address(mdp_site x, int i=0) const {
#ifdef CHECK_ALL
    if(!(x.is_here())) {
      error("You are looking for a site that is not here");
    }
#endif
    return m+x.idx*field_components+i;
  }
  /// shifts the entire fields in direction mu of i steps 
  /// (i can be positive or negative)
  /// note that if i=1, field(x-mu) is assigned to field(x)
  /// function requires communication
  void shift(int i, int mu) {
    mdp_field tmp(lattice(),field_components);
    mdp_site x(lattice());
    int k;
    while(i!=0) {
      update();
      if(i==+1) {
	forallsites(x)
	  for(k=0; k<field_components; k++)
	    tmp(x,k)=(*this)(x-mu,k); // mind here
	i--;
      } else if(i==-1) {
	forallsites(x)
	  for(k=0; k<field_components; k++)
	    tmp(x,k)=(*this)(x+mu,k); // mind here
	i++;
      }
      (*this)=tmp;
    }	      
  }

  void operator= (const mdp_field &a) {
    if(&lattice()!=&a.lattice() || 
       size!=a.size || 
       field_components!=a.field_components) 
      error("mdp_field: operator=() incompatible fields");
    mdp_int i=0;
    for(; i<size; i++) m[i]=a.m[i];
  }
  void operator= (const T a) {
    for(mdp_int i=0; i<size; i++) m[i]=a;
  }
  void operator+=(const mdp_field &a) {
    for(mdp_int i=0; i<size; i++) m[i]+=a.m[i];
  }
  void operator-=(const mdp_field &a) {
    for(mdp_int i=0; i<size; i++) m[i]-=a.m[i];
  }
  template <class T2> 
    void operator*=(const T2 a) {
    for(mdp_int i=0; i<size; i++) m[i]*=a;
  }
  template <class T2>
    void operator/=(const T2 a) {
    for(mdp_int i=0; i<size; i++) m[i]/=a;
  }
  /// returns by reference the lattice this field is defined on
  inline mdp_lattice &lattice() const {
    return *ptr;
  }
  /// returns the total memory in bytes occupied by the field
  mdp_int field_size() {
    return 
      lattice().size()*field_components*Tsize;
  }
  /// returns the total space in bytes required to store the field
  mdp_int file_size() {
    return 
      sizeof(mdp_field_file_header)+
      field_size();
  }
  /// only used by mdp_field::load() and mdp_field::save()
  int where_global(mdp_int i) {
    int x[10];
    lattice().global_coordinate(i, x);
    return (*(lattice().where))(x, lattice().ndim, lattice().nx);
  }
  void switch_endianess_4bytes() {
    // I am not sure if this workd for complex<double>
    int32_t *p;
    uint i;

    if(Tsize*field_components % 4 !=0) error("Field not % 4");
    mdp_site x(lattice());
    forallsitesandcopies(x) {
      p=(int32_t*) address(x);
      for(i=0; i<Tsize*field_components/4; i++) {	
	switch_endianess_byte4(*(p+i));
      }
    }
  }
  void switch_endianess_8bytes() {
    // I am not sure if this workd for complex<double>
    int64_t *p;
    uint i;

    if(Tsize*field_components % 8 !=0) error("Field not % 8");
    mdp_site x(lattice());
    forallsitesandcopies(x) {
      p=(int64_t*) address(x);
      for(i=0; i<Tsize*field_components/8; i++) {
	cout << '.';
	switch_endianess_byte8(*(p+i));
      }
    }
  }

  // ////////////////////////////////////
  // functions to access member variables
  // ////////////////////////////////////

  /// lattice size in units of sizeof(T)
  inline mdp_int global_size() {
    return field_components*lattice().global_volume();
  }
  inline mdp_int physical_size() {
    return size;
  }
  inline mdp_int size_per_site() {
    return field_components;
  }
  inline mdp_int physical_local_start(int i=2) {
    if(i==2) i=0;
    return field_components*lattice().start[ME][i];
  }
  inline mdp_int physical_local_stop(int i=2) {
    if(i==2) i=1;
    return field_components*lattice().stop[ME][i];
  }
  inline T* physical_address(mdp_int i=0) {
    return m+i; 
  }

  /// the most important communication function in MDP.
  /// it must be called after each field variables are modified.
  /// it restores the syncronization between parallel processes.
  void update(int np=2, int d=-1, int size=1);

  // ////////////////////////////////////////////////////////////
  // IO funcyions: load(filename, processIO, buffersize)
  //               save(filename, processIO, buffersize)   
  // filename should possibly include the path.
  // processIO is the process that physically perform the IO.
  // buffersize if the size of the buffer associated to the
  //   communication to each process. buffersize*Nproc 
  //   should fit in the memory of processIO.
  //   By default buffersize=1024 and it works reasonably fast.
  // ///////////////////////////////////////////////////////////

  bool load(string filename, 
	    int processIO=0, 
	    mdp_int max_buffer_size=1024, 
	    bool load_header=true,
	    mdp_int skip_bytes=0,
	    bool (*user_read)(FILE*, void*, mdp_int, mdp_int, mdp_int, const mdp_lattice&)=0,
	    bool try_switch_endianess=true);
  
  
  bool save(string filename, 
	    int processIO=0, 
	    mdp_int max_buffer_size=1024, 
	    bool load_header=true,
	    mdp_int skip_bytes=0,
	    bool (*user_write)(FILE*, void*, mdp_int, mdp_int, mdp_int, const mdp_lattice&)=0);

  bool save_vtk(string filename, 
		int t=-1,
		int component=-1,
		int processIO=0,
		bool ASCII=false); 
  
#ifdef INCLUDE_DEPRECATED_IO
  void load(char filename[], 
	    int processIO, 
	    mdp_int max_buffer_size, 
	    char *header, 
	    mdp_int header_size=0,
	    mdp_int (*sort_x)(mdp_lattice&,mdp_int)=0,
	    int auto_switch_endianess=true);
  
  void save(char filename[], 
	    int processIO, 
	    mdp_int max_buffer_size,
	    char *header,
	    mdp_int header_size=0,
	    mdp_int (*sort_x)(mdp_lattice&,mdp_int)=0,
	    char *mode="w");
#endif
  
};



