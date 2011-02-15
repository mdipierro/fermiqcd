/////////////////////////////////////////////////////////////////
/// @file mdp_deprecatedIO.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Old functions for file IO now deprecated
///
/// Licensed under GPL2 license
/// Read attached license in file mdp_license.pdf
/// This file cannot be distributed without file mdp_license.pdf
//////////////////////////////////////////////////////////////////

template<class T>
void mdp_field<T>::load(char filename[], 
			int processIO, 
			mdp_int max_buffer_size, 
			char *header, 
			mdp_int header_size,
			mdp_int (*sort_x)(mdp_lattice&,mdp_int),
			int auto_switch_endianess) {
  mdp_int idx_gl, nvol_gl=lattice().nvol_gl, k;
  double mytime=mpi.time();
  int try_switch_endianess=false;
  if(ME==processIO) {
    mdp_int *buffer_size=new mdp_int[Nproc];
      mdp_array<T,3> large_buffer(Nproc,max_buffer_size,field_components);
      T *short_buffer=new T[field_components];
      int process;
      mdp_request request;

      for(process=0; process<Nproc; process++) buffer_size[process]=0;
      printf("Loading file %s from process %i (buffer = %li sites)\n", 
	     filename, processIO, max_buffer_size); 
      fflush(stdout);
      FILE *fp=fopen(filename, "rb");
      if(fp==0) error("Unable to open file");

      // #ifdef INSERT_HEADER
      int i;
      if(strcmp(header,"NATIVE")==0) {
	error("NATIVE HEADER IN DEPRECATED FUNCTION NOT SUPPORTED ANY MORE");
      } else if(header!=0 && strcmp(header, "NOHEADER")!=0) {
	if(fread(header, sizeof(char), header_size, fp)!=
	   header_size) error("Unable to read file header");
      }
      // #endif
      for(idx_gl=0; idx_gl<nvol_gl; idx_gl++) {
	process=where_global(idx_gl);
	if(process!=NOWHERE) {
	  if(sort_x!=0)
	    if(fseek(fp, sort_x(lattice(),idx_gl)*
		     Tsize*field_components+header_size, SEEK_SET)<0) 
	      error("unexpected end of file");
	  if((fread(short_buffer, Tsize, field_components, fp)-
	      field_components)!=0)
	    error("unexpected end of file");
	}
	if((process!=NOWHERE) && (process!=processIO)) {
	  for(k=0; k<field_components; k++)
	    large_buffer(process,buffer_size[process],k)=short_buffer[k];
	  buffer_size[process]++;
	  if(buffer_size[process]==max_buffer_size) {
	    mpi.put(&(large_buffer(process,0,0)), 
		    max_buffer_size*field_components, process, request); 
	    mpi.wait(request);
	    buffer_size[process]=0;
	  }
	  if(idx_gl==nvol_gl-1) 
	    for(process=0; process<Nproc; process++)
	      if((process!=ME) && 
		 (buffer_size[process]!=max_buffer_size) &&
		 (buffer_size[process]>0)) {
		mpi.put(&(large_buffer(process,0,0)), 
			buffer_size[process]*field_components, 
			process, request); 
		mpi.wait(request);
	      }
	}
	if(process==processIO) {
	  for(k=0; k<field_components; k++)
	    *(m+lattice().lg[idx_gl]*field_components+k)=short_buffer[k];
	}
      }
      delete[] buffer_size;
      delete[] short_buffer;
      fclose(fp);
    } else {
      int  process;
      mdp_int buffer_size=0, idx;
      mdp_int *local_index=new mdp_int[max_buffer_size];
      mdp_array<T,2> local_buffer(max_buffer_size,field_components);
      for(idx_gl=0; idx_gl<nvol_gl; idx_gl++) {
	process=where_global(idx_gl);
	if(process==ME) {
	  local_index[buffer_size]=lattice().lg[idx_gl];
	  buffer_size++;
	}
	if((buffer_size==max_buffer_size) || 
	   ((idx_gl==nvol_gl-1) && (buffer_size>0))) {
	  mpi.get(&(local_buffer(0,0)), buffer_size*field_components, processIO);
	  for(idx=0; idx<buffer_size; idx++)
	    for(k=0; k<field_components; k++)
	      *(m+local_index[idx]*field_components+k)=local_buffer(idx,k);
	  buffer_size=0;
	}
      }
      delete[] local_index;
    }
    update();
    if(ME==0 && mdp_shutup==false) {
      printf("... Loading time: %f (sec)\n", mpi.time()-mytime);
      fflush(stdout);
    }
    if(try_switch_endianess==true && auto_switch_endianess==true) {
      if(ME==0) printf("Trying to switch endianess.");
      switch_endianess_4bytes();
    }
}

template<class T>
void mdp_field<T>::save(char filename[], 
			int processIO,
			mdp_int max_buffer_size,
			char *header, 
			mdp_int header_size,
			mdp_int (*sort_x)(mdp_lattice&,mdp_int),
			char *mode) {
  mdp_int idx_gl, nvol_gl=lattice().nvol_gl, k;
  double mytime=mpi.time();
    if(ME==processIO) {
      mdp_int *buffer_size=new mdp_int[Nproc];
      mdp_int *buffer_ptr =new mdp_int[Nproc];
      mdp_array<T,3> large_buffer(Nproc,max_buffer_size,field_components);
      T *short_buffer=new T[field_components];
      int process;
      for(process=0; process<Nproc; process++) buffer_ptr[process]=0;
      printf("Saving file %s from process %i (buffer = %li sites)\n", 
	     filename, processIO, max_buffer_size); 
      fflush(stdout);
      FILE *fp=fopen(filename, mode);
      if(fp==0) error("Unable to open file");      


      // #ifdef INSERT_HEADER
      int i;
      if(strcmp(header,"NATIVE")==0) {
	error("NATIVE HEADER IN DEPRECATED FUNCTION NOT SUPPORTED ANY MORE");
      } else if(header!=0 && strcmp(header, "NOHEADER")!=0) {
	if(fwrite(header, sizeof(char), header_size, fp)!=
	   header_size) error("Unable to write file header");
      }
      // #endif
      for(idx_gl=0; idx_gl<nvol_gl; idx_gl++) {
	process=where_global(idx_gl);
	if((process!=NOWHERE) && (process!=processIO)) {
	  if(buffer_ptr[process]==0) {
	    mpi.get(buffer_size[process], process);
	    mpi.get(&(large_buffer(process,0,0)), 
		    buffer_size[process]*field_components, process);
	  }
	  for(k=0; k<field_components; k++)
	    short_buffer[k]=large_buffer(process,buffer_ptr[process],k);
	  buffer_ptr[process]++;
	  if(buffer_ptr[process]==buffer_size[process]) buffer_ptr[process]=0; 
	}
	if(process==processIO) {
	  for(k=0; k<field_components; k++)
	    short_buffer[k]=*(m+lattice().lg[idx_gl]*field_components+k);
	}
	if(process!=NOWHERE) {
	  if(sort_x!=0)
	    if(fseek(fp, sort_x(lattice(),idx_gl)*
		     Tsize*field_components+header_size, SEEK_SET)<0) 
	      error("unexpected end of file");
	  if((fwrite(short_buffer, Tsize, field_components, fp)-
	     field_components)!=0)
	    error("I cannot write on the file. I am confused !?!?");
	}
      }
      if(strcmp(header,"NATIVE")==0) 
	fprintf(fp, "\n\n [ MDP Standard File Format ]\n");
      delete[] buffer_size;
      delete[] buffer_ptr;
      delete[] short_buffer;
      fclose(fp);
    } else {
      int  process;
      mdp_int buffer_size=0, idx, idx_gl;
      mdp_int *local_index=new mdp_int[max_buffer_size];
      mdp_array<T,2> local_buffer(max_buffer_size,field_components);
      mdp_request request;
      for(idx_gl=0; idx_gl<nvol_gl; idx_gl++) {
	process=where_global(idx_gl);
	if(process==ME) {
	  local_index[buffer_size]=lattice().lg[idx_gl];
	  buffer_size++;
	}
	if((buffer_size==max_buffer_size) || 
	   ((idx_gl==nvol_gl-1) && (buffer_size>0))) {
	  for(idx=0; idx<buffer_size; idx++)
	    for(k=0; k<field_components; k++)
	      local_buffer(idx,k)=*(m+local_index[idx]*field_components+k);
	  mpi.put(buffer_size, processIO, request);
	  mpi.wait(request);
	  mpi.put(&(local_buffer(0,0)), buffer_size*field_components, 
		  processIO, request);
	  mpi.wait(request);
	  buffer_size=0;
	}
      }
      delete[] local_index;
    }
    if(ME==0 && mdp_shutup==false) {
      printf("... Saving time: %f (sec)\n", mpi.time()-mytime);
      fflush(stdout);
    }
}


