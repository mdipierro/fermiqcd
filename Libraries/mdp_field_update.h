/////////////////////////////////////////////////////////////////
/// @file mdp_field_update.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains mdp_field::update()
///
/// Licensed under GPL2 license
/// Read attached license in file mdp_license.pdf
/// This file cannot be distributed without file mdp_license.pdf
//////////////////////////////////////////////////////////////////

/// The only communication function for a field object
/// To be invoked every time field variables are assigned and
/// Need to be sinchronized between the parallel processes
template<class T>
void mdp_field<T>::update(int np, int d, int ncomp) {
    T *dynamic_buffer=0;
    T *where_to=0;
    mpi.comm_time-=mpi.time();
    mdp_request request;
    mdp_int start_to_send,dp,process,length,idx;
    int k, ni, nf;
    if(d==-1) {d=0; ncomp=field_components;}
    if((ncomp==field_components) && (d!=0)) 
      error("update(): packet is too big");
    if(np<2) ni=nf=np;
    else     {ni=0; nf=1; }
    for(dp=1; dp<Nproc; dp++) {       
      process=(ME+dp) % Nproc;
      if(np<2)  {
	length=lattice().len_to_send[process][np];
      } else {
	length=lattice().len_to_send[process][0]+lattice().len_to_send[process][1];
      }
      if(np==1) start_to_send=lattice().len_to_send[process][0];
      else      start_to_send=0;
      if(length>0) {
	dynamic_buffer=new T[length*ncomp];
	for(idx=0; idx<length; idx++)
	  for(k=0; k<ncomp; k++)
	    dynamic_buffer[idx*ncomp+k]=
	      *(m+lattice().to_send[process][start_to_send+idx]*
		field_components+d*ncomp+k);	
	mpi.put(dynamic_buffer, length*ncomp, process, request);
	cout.flush();
      } else {
	dynamic_buffer=0;
      }

      process=(ME-dp+Nproc) % Nproc;
      length=lattice().stop[process][nf]-lattice().start[process][ni];
      if(length>0) {
	if(ncomp==field_components) {
	  where_to=m+lattice().start[process][ni]*field_components;	  
	  mpi.get(where_to, length*field_components, process);
	  where_to=0;
	} else {
	  where_to=new T[length*ncomp];
	  mpi.get(where_to, length*ncomp, process);
	  for(idx=0; idx<length; idx++)
	    for(k=0; k<ncomp; k++)
	      *(m+(lattice().start[process][ni]+idx)*field_components+
		d*ncomp+k)=where_to[idx*ncomp+k]; 
	  delete[] where_to;
	}
      }

      process=(ME+dp) % Nproc;      
      if(dynamic_buffer!=0) {
	mpi.wait(request);
	delete[] dynamic_buffer;
      }
    }
    mpi.comm_time+=mpi.time();
  }
