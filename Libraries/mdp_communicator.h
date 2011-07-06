/////////////////////////////////////////////////////////////////
/// @file mdp_communicator.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains declaration of class mdp_array
///
/// Licensed under GPL2 license
/// Read attached license in file mdp_license.pdf
/// This file cannot be distributed without file mdp_license.pdf
//////////////////////////////////////////////////////////////////

#ifdef PARALLEL
#include "mpi.h"
     
typedef MPI_Request mdp_request;

class mdp_communicator : public mdp_log {
 private:  
  MPI_Comm communicator;
  double   mytime; // total time
  int      wormholes_open;
  int      my_id;
  int      my_nproc;
 public:
  double   comm_time;  // time spent in communications
  void open_wormholes(int argc, char** argv, 
		      MPI_Comm communicator_=MPI_COMM_WORLD) {
    if(wormholes_open) return;
    communicator=communicator_;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(communicator, &(my_id));
    MPI_Comm_size(communicator, &(my_nproc));
    print=true;
    for(int i=0; i<nproc(); i++) {
      if(me()==i) 
	(*this) << "MPI PROCESS " << me() 
		<< " of " << nproc() 
		<< " STARTING ON " << getname() << "\n";
      barrier();
    }
    if(me()==0) print=true; else print=false;

    begin_function("PROGRAM");
    begin_function("open_wormholes");
    (*this) << 
      "<head>\n"
      "*************************************************\n"
      "* Starting [Matrix Distributed Processing]      *\n"
      "* created by Massimo Di Pierro                  *\n"
      "* copyrighted by www.metacryption.com           *\n"
      "*************************************************\n"
      "</head>\n";


    reset_time();
    wormholes_open=true; 

    double a,b;
    getcpuusage(a,b);
    end_function("open_wormholes");
  }
  int tag(int i, int j) {
    return i*nproc()+j;
  }
  inline const int me() {
    return my_id;
  }
  inline const int nproc() {
    return my_nproc;
  }
  void print_stats() {
#ifndef NO_POSIX
    int i;
    double *a=new double[nproc()];
    double *b=new double[nproc()];
    double *c=new double[nproc()];
    for(i=0; i<nproc(); i++) a[i]=b[i]=c[i]=0; 
    getcpuusage(a[me()],b[me()]);
    c[me()]=100.0*comm_time/(time());
    add(a,nproc());
    add(b,nproc());
    add(c,nproc());
    char buffer[256];
    for(i=0; i<nproc(); i++) {
      sprintf(buffer,
	      "* Process %i stats: CPU=%.2f%% PROCESS=%.2f%% COMM=%.2f%%\n",
	      i, a[i],b[i],c[i]);
      (*this) << buffer;
    }
    (*this) << "* (above numbers make no sense under windows)\n";
    delete[] a;
    delete[] b;
    delete[] c;
#endif
  }
  void close_wormholes() {
    begin_function("close_wormholes");
    (*this) << 
      "<foot>\n"
      "*************************************************\n"
      "* Ending [Matrix Distributed Processing]        *\n";
    print_stats();
    (*this) << 
      "*************************************************\n"
      "</foot>\n";
    (*this) << "PROCESS " << me() << " ENDING AFTER " << time() << " sec.\n";
    wormholes_open=false;
    MPI_Finalize();
    end_function("close_wormholes");
    end_function("PROGRAM");
  }
  void abort() {
    MPI_Abort(communicator,1);
  }
  mdp_communicator() {
    wormholes_open=false;
  }
  virtual ~mdp_communicator() {
  }
  template<class T>
    void put(T &obj, int destination) {
    mdp_request r;
    MPI_Isend(&obj, sizeof(T)/sizeof(char), MPI_CHAR, destination, 
	     tag(me(),destination), communicator, &r);
    wait(r);
  }
  template<class T>
    void put(T &obj, int destination, mdp_request &r) {
    MPI_Isend(&obj, sizeof(T)/sizeof(char), MPI_CHAR, destination, 
	     tag(me(),destination), communicator, &r);
  }
  template<class T>
    void get(T &obj, int source) {
    MPI_Status  status;
    MPI_Recv(&obj, sizeof(T)/sizeof(char), MPI_CHAR, source, 
	     tag(source,me()), communicator, &status);     
  }
  template<class T>
    void put(T* objptr, mdp_int length, int destination) {
    mdp_request r;
    MPI_Isend(objptr, length*sizeof(T)/sizeof(char), MPI_CHAR, destination, 
	      tag(me(),destination), communicator, &r);
    wait(r);
  }
  template<class T>
    void put(T* objptr, mdp_int length, int destination, mdp_request &r) {
    if(MPI_Isend(objptr, length*sizeof(T)/sizeof(char), MPI_CHAR, destination, 
		 tag(me(),destination), communicator, &r)!=MPI_SUCCESS)
      error("Failure to send");
  }
  template<class T>
    void get(T* objptr, mdp_int length, int source) {
    MPI_Status  status;
    if(MPI_Recv(objptr, length*sizeof(T)/sizeof(char), MPI_CHAR, source, 
	     tag(source,me()), communicator, &status)!=MPI_SUCCESS) 
      error("Failure to receive");
  }
  void wait(mdp_request &r) {
    MPI_Status status;
    MPI_Wait(&r, &status);
  }
  void wait(mdp_request *r, int length) {
    MPI_Status status;
    MPI_Waitall(length, r, &status);
  }
  // Compatibility functions
  // obj1 is input, obj2 is output
  void add(float &obj1, float &obj2) {
    MPI_Allreduce(&obj1, &obj2, 1, MPI_FLOAT, MPI_SUM, communicator);
  }
  void add(float* obj1, float *obj2, mdp_int length) {
    MPI_Allreduce(obj1, obj2, length, MPI_FLOAT, MPI_SUM, communicator);
  }
  void add(double &obj1, double &obj2) {
    MPI_Allreduce(&obj1, &obj2, 1, MPI_DOUBLE, MPI_SUM, communicator);
  }
  void add(double* obj1, double *obj2, mdp_int length) {
    MPI_Allreduce(obj1, obj2, length, MPI_DOUBLE, MPI_SUM, communicator);
  }
  void add(mdp_int &obj1) {
    mdp_int obj2=0;
    MPI_Allreduce(&obj1, &obj2, 1, MPI_LONG, MPI_SUM, communicator);
    obj1=obj2;
  }
  void add(float &obj1) {
    float obj2=0;
    MPI_Allreduce(&obj1, &obj2, 1, MPI_FLOAT, MPI_SUM, communicator);
    obj1=obj2;
  }
  void add(double &obj1) {
    double obj2=0;
    MPI_Allreduce(&obj1, &obj2, 1, MPI_DOUBLE, MPI_SUM, communicator);
    obj1=obj2;
  }
  void add(mdp_int *obj1, mdp_int length) {
    mdp_int i;
    mdp_int *obj2=new mdp_int[length];
    for(i=0; i<length; i++) obj2[i]=0;
    MPI_Allreduce(obj1, obj2, length, MPI_LONG, MPI_SUM, communicator);
    for(i=0; i<length; i++) obj1[i]=obj2[i];
    delete[] obj2;
  }
  void add(float *obj1, mdp_int length) {
    mdp_int i;
    float *obj2=new float[length];
    for(i=0; i<length; i++) obj2[i]=0;
    MPI_Allreduce(obj1, obj2, length, MPI_FLOAT, MPI_SUM, communicator);
    for(i=0; i<length; i++) obj1[i]=obj2[i];
    delete[] obj2;
  }
  void add(double *obj1, mdp_int length) {
    mdp_int i;
    double *obj2=new double[length];
    for(i=0; i<length; i++) obj2[i]=0;
    MPI_Allreduce(obj1, obj2, length, MPI_DOUBLE, MPI_SUM, communicator);
    for(i=0; i<length; i++) obj1[i]=obj2[i];
    delete[] obj2;
  }
  void add(mdp_complex &obj1) {
    mdp_complex obj2=0;
#if !defined(USE_DOUBLE_PRECISION)
    MPI_Allreduce(&obj1, &obj2, 2, MPI_FLOAT, MPI_SUM, communicator);
#else
    MPI_Allreduce(&obj1, &obj2, 2, MPI_DOUBLE, MPI_SUM, communicator);
#endif
    obj1=obj2;
  }
  void add(mdp_complex *obj1, mdp_int length) {
    mdp_int i;
    mdp_complex *obj2=new mdp_complex[length];
    for(i=0; i<length; i++) obj2[i]=0;
#if !defined(USE_DOUBLE_PRECISION)
    MPI_Allreduce(obj1, obj2, 2*length, MPI_FLOAT, MPI_SUM, communicator);
#else
    MPI_Allreduce(obj1, obj2, 2*length, MPI_DOUBLE, MPI_SUM, communicator);
#endif
    for(i=0; i<length; i++) obj1[i]=obj2[i];
    delete[] obj2;
  }
  void add(mdp_matrix &a) {
    add(a.address(), a.rowmax()*a.colmax());
  }
  void add(mdp_matrix *a, mdp_int length) {
    mdp_int i;
    for(i=0; i<length; i++) add(a[i]);
  }
  template<class T>
    void add(vector<T> &a) {
    add(&a[0],a.size());
  }

  void barrier() {
    MPI_Barrier(communicator);
  }
  template<class T>
    void broadcast(T &obj, int p) {
    MPI_Bcast(&obj, sizeof(T)/sizeof(char), MPI_CHAR, p, communicator); 
  }
  template<class T>
    void broadcast(T* obj, mdp_int length, int p) {
    MPI_Bcast(obj, length*sizeof(T)/sizeof(char), MPI_CHAR, p,communicator); 
  }
  void reset_time() {
    mytime=MPI_Wtime();
    comm_time=0;
  }
  double time() {
    return MPI_Wtime()-mytime;
  }
};

#else
#include "time.h"
typedef int mdp_request;

/// @brief DO NOT INSTANTIATE use object mdp instead
///
/// Example:
/// @verbatim
/// int main(int argc, char**argv) {
///    mdp.open_wormholes(argc,argv);
///    // your code here
///    mdp << 3.14 << endl;  // only process 0 prints
///    mdp.close_wormholes();
///    return 0;
/// }
/// @endverbatim
class mdp_communicator : public mdp_log {
 private:
#ifndef NO_POSIX
  mdp_psim *nodes;
#endif
  int    communicator;
  int    wormholes_open;
  double  mytime; // total time
  double MPI_Wtime() {
    return walltime();
  }
  int    my_id;
  int    my_nproc;
 public:
  double  comm_time;  /// time spent in communications
  mdp_communicator() {
    wormholes_open=false;
  }
  template<class T>
    void put(T &obj, int destination) {
#ifndef NO_POSIX
    nodes->send(destination,"",obj);
#endif
  }
  template<class T>
    void put(T &obj, int destination, mdp_request &r) {
#ifndef NO_POSIX
    nodes->send(destination,"",obj);
#endif
  }
  template<class T>
    void get(T &obj, int source) {
#ifndef NO_POSIX
    nodes->recv(source,"",obj);
#endif
  }
  template<class T>
    void put(T* objptr, mdp_int length, int destination) {
#ifndef NO_POSIX
    nodes->send(destination,"",objptr, length);
#endif
  }
  template<class T>
    void put(T* objptr, mdp_int length, int destination, mdp_request &r) {
#ifndef NO_POSIX
    nodes->send(destination,"",objptr, length);
#endif
  }
  template<class T>
    void get(T* objptr, mdp_int length, int source) {
#ifndef NO_POSIX
    nodes->recv(source,"",objptr,length);
#endif
  }
  // ////////////////////////////////
  // Compatibility functions
  // obj1 is input, obj2 is output
  void add(float &obj1, float &obj2) {
#ifndef NO_POSIX
    obj2=nodes->add(obj1);
#endif
  }
  void add(float* obj1, float *obj2, mdp_int length) {
#ifndef NO_POSIX
    for(int i=0; i<length; i++) {
      obj2[i]=nodes->add(obj1[i]);
    }
#endif
  }
  void add(double &obj1, double &obj2) {
#ifndef NO_POSIX
    obj2=nodes->add(obj1);
#endif
  }
  void add(double* obj1, double *obj2, mdp_int length) {
#ifndef NO_POSIX
    for(int i=0; i<length; i++) {
      obj2[i]=nodes->add(obj1[i]);
    }
#endif
  }
  // end compatilibility functions
  // /////////////////////////////////
  void add(mdp_int &obj1) {
#ifndef NO_POSIX
    obj1=nodes->add(obj1);
#endif
  }
  void add(float &obj1) {
#ifndef NO_POSIX
    obj1=nodes->add(obj1);
#endif
  }
  void add(double &obj1) {
#ifndef NO_POSIX
    obj1=nodes->add(obj1);
#endif
  }
  void add(mdp_int *obj1, mdp_int length) {
#ifndef NO_POSIX
    for(int i=0; i<length; i++) {
      obj1[i]=nodes->add(obj1[i]);
    }
#endif
  }
  /*
  void add(mdp_int *obj1, mdp_int length) {
#ifndef NO_POSIX
    for(mdp_int i=0; i<length; i++) obj1[i]=nodes->add(obj1[i]);
#endif
  }
  */
  void add(float *obj1, mdp_int length) {
#ifndef NO_POSIX
    for(mdp_int i=0; i<length; i++) obj1[i]=nodes->add(obj1[i]);
#endif
  }
  void add(double *obj1, mdp_int length) {
#ifndef NO_POSIX
    for(mdp_int i=0; i<length; i++) obj1[i]=nodes->add(obj1[i]);
#endif
  }
  void add(mdp_complex &obj1) {
#ifndef NO_POSIX
    obj1=nodes->add(obj1);
#endif
  }
  void add(mdp_complex *obj1, mdp_int length) {
#ifndef NO_POSIX
    for(mdp_int i=0; i<length; i++) obj1[i]=nodes->add(obj1[i]);
#endif
  }
  void add(mdp_matrix &a) {
#ifndef NO_POSIX
    for(mdp_int i=0; i<a.size(); i++) 
      a.address()[i]=nodes->add(a.address()[i]);
#endif
  }
  void add(mdp_matrix *a, mdp_int length) {
#ifndef NO_POSIX
    for(mdp_int j=0; j<length; j++)
      for(mdp_int i=0; i<a[j].size(); i++) 
	a[j].address()[i]=nodes->add(a[j].address()[i]);
#endif
  }
  template<class T>
    void add(vector<T> &a) {
    add(&a[0],a.size());
  }
  template<class T>
    void broadcast(T &obj, int p) {
#ifndef NO_POSIX
    nodes->broadcast(p,obj);
#endif
  }
  template<class T>
    void broadcast(T* obj, mdp_int length, int p) {
#ifndef NO_POSIX
    nodes->broadcast(p,obj,length);
#endif
  }
  void wait(mdp_request &r) {}
  void wait(mdp_request *r, int length) {}
  inline const int me() {
    return my_id;
  }
  inline const int nproc() {
    return my_nproc;
  }
  void barrier() {}
  int tag(int i, int j) {
    return i*nproc()+j;
  }
  void reset_time() {
    mytime=MPI_Wtime();
    comm_time=0;
  }
  /// returns the time in seconds since call to mdp_communicator::open_wormholes
  double time() { 
    return MPI_Wtime()-mytime;
  }
  /// starts communications 
  /// parses command line argument for MPI or PSIM parameters
  void open_wormholes(int argc, char** argv) {
    if(wormholes_open) return;

#ifndef NO_POSIX      
    nodes=new mdp_psim(argc, argv);
    my_id=nodes->id();
    my_nproc=nodes->nprocs();
#else
    my_id=0;
    my_nproc=1;
#endif
    if(me()==0) print=true; else print=false;
    begin_function("PROGRAM");
    begin_function("open_wormholes");
    (*this) << 
      "<head>\n"
      "*************************************************\n"
      "* Starting [Matrix Distributed Processing]      *\n"
      "* created by Massimo Di Pierro                  *\n"
      "* copyrighted by www.metacryption.com           *\n"
      "*************************************************\n"
      "</head>\n";
    reset_time();
    double a,b;
    getcpuusage(a,b);
    wormholes_open=true;
    end_function("open_wormholes");
  }
  /// prints statistics about parallel processes
  void print_stats() {
#ifndef NO_POSIX
    int i;
    double *a=new double[nproc()];
    double *b=new double[nproc()];
    double *c=new double[nproc()];
    for(i=0; i<nproc(); i++) a[i]=b[i]=c[i]=0; 
    getcpuusage(a[me()],b[me()]);
    c[me()]=100.0*comm_time/(time());
    add(a,nproc());
    add(b,nproc());
    add(c,nproc());
    char buffer[256];
    for(i=0; i<nproc(); i++) {
      sprintf(buffer,
	      "* Process %i stats: CPU=%.2f%% PROCESS=%.2f%% COMM=%.2f%%\n",
	      i, a[i],b[i],c[i]);
      (*this) << buffer;
    }
    (*this) << "* (above numbers make no sense under windows)\n";
    delete[] a;
    delete[] b;
    delete[] c;
#endif
  }
  /// closes parallel communications
  void close_wormholes() {
    begin_function("close_wormholes");
    (*this) << 
      "<foot>\n"
      "*************************************************\n"
      "* Ending [Matrix Distributed Processing]        *\n";
    print_stats();
    (*this) << 
      "*************************************************\n"
      "</foot>\n";
    (*this) << "PROCESS " << me() << " ENDING AFTER " << time() << " sec.\n";
    wormholes_open=false;
    end_function("close_wormholes");
    end_function("PROGRAM");
#ifndef NO_POSIX
    if(nodes) delete nodes;
#endif
  }
  /// forces the process to exit(-1)
  void abort() {
    (*this) << "calling abort...";
    exit(-1);
  }
};

#endif

/// the only communicator object
mdp_communicator mdp; 

/// alias for mdp
mdp_communicator& mpi=mdp; 

void _mpi_error_message(string a, string b="unkown", int c=0) {
  mpi.error_message(a,b,c);
}

/// Logs in xml the start of a function with message s
inline void begin_function(string s) {
  mpi.begin_function(s);
}

/// Logs in xml the end of a function with message s
inline void end_function(string s) {
  mpi.end_function(s);
}

