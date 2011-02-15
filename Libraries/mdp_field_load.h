/////////////////////////////////////////////////////////////////
/// @file mdp_field_load.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains file IO operations for class mdp_field
///
/// Licensed under GPL2 license
/// Read attached license in file mdp_license.pdf
/// This file cannot be distributed without file mdp_license.pdf
//////////////////////////////////////////////////////////////////

/// Auxiliary function
bool mdp_default_user_read(FILE *fp,
			   void* p,
			   mdp_int psize,
			   mdp_int header_size,
			   mdp_int position,
			   const mdp_lattice &lattice) {
  if(fseek(fp, (size_t) position*psize+header_size, SEEK_SET) ||
     fread(p, psize,1, fp)!=1) return false;
  return true;
}


/// Best way to load a field
template<class T>
bool mdp_field<T>::load(string filename, 
			int processIO, 
			mdp_int max_buffer_size, 
			bool load_header,
			mdp_int skip_bytes,
			bool (*user_read)(FILE*, void*, mdp_int, mdp_int, mdp_int, const mdp_lattice&),
			bool try_swicth_endianess) {

  filename=latest_file(filename);
  if(filename=="?") return false;
  mdp_int header_size=0;
  size_t idx_gl, nvol_gl=lattice().nvol_gl, k;
  size_t psize=field_components*Tsize;
  double mytime=mpi.time();
  bool reversed_header_endianess=false;
  struct stat statbuf;
  if(ME==processIO) {
    mdp_int *buffer_size=new mdp_int[Nproc];
    mdp_array<T,3> large_buffer(Nproc,max_buffer_size,field_components);
    T *short_buffer=new T[field_components];
    int process;
    mdp_request request;

    for(process=0; process<Nproc; process++) buffer_size[process]=0;
    cout << "Loading file " << filename 
	 << " from process " << processIO 
	 << " (buffer = " << max_buffer_size << " sites)" << '\n'; 
    fflush(stdout);
    stat(filename.c_str(),&statbuf);
    int total_size = statbuf.st_size;
    FILE *fp=fopen(filename.c_str(), "rb");
    if(fp==0) error("Unable to open file");
    
    int i;
    if(load_header) {
      mdp_field_file_header tmp_header;
      header_size=sizeof(mdp_field_file_header);

      if(fseek(fp, skip_bytes, SEEK_SET) || 
	   fread(&tmp_header, header_size, 1, fp)!=1) {
	fprintf(stderr,"mdp_field.load(): Unable to load file header\n");
	return false;
      }
      
      reversed_header_endianess=switch_header_endianess(tmp_header);

      cout << "reverse: " << reversed_header_endianess << endl;

      if(tmp_header.endianess!=header.endianess) 
	fprintf(stderr, "Unrecognized endianess... trying to read anyway\n");

      // UGLY BUT FIXES INCOMPATIBIITY
      int actual_size = tmp_header.box[0];
      for(int d=1; d<tmp_header.ndim; d++) actual_size*=tmp_header.box[d];
      tmp_header.sites = actual_size;
      header_size += total_size-(tmp_header.bytes_per_site*actual_size+header_size);

      if(tmp_header.ndim!=header.ndim) {
	fprintf(stderr,"mdp_field.load(): wrong ndim\n");
	return false;
      }
      for(i=0; i<lattice().ndim; i++)  
	if(tmp_header.box[i]!=header.box[i]) {
	  fprintf(stderr,"mdp_file.load(): wrong lattice size\n");
	  return false;
	}
      if(tmp_header.bytes_per_site!=header.bytes_per_site) {
	fprintf(stderr, "mdp_file.load(): wrong type of field (%i bytes per site?)\n", tmp_header.bytes_per_site );
	return false;
      }
      if(tmp_header.sites!=header.sites) {
	fprintf(stderr,"mdp_field.load(): wrong number of sites\n");
	return false;
      }
      header=tmp_header;
      
    }
    
    skip_bytes+=header_size;

    bool exception = false;

    fseek(fp, skip_bytes, SEEK_SET);
    for(idx_gl=0; idx_gl<nvol_gl; idx_gl++) {
      process=where_global(idx_gl);
      if(process!=NOWHERE) {
	if(user_read) {
	  if(!user_read(fp, short_buffer, 
			field_components*Tsize, 
			skip_bytes, 
			idx_gl, lattice()))
	    error("unexpected end of file");
	} else {
	  if(exception && fseek(fp, idx_gl*psize+skip_bytes, SEEK_SET)) {
	    cout << "debug info: " << idx_gl*psize+skip_bytes << " " << psize << endl;
	    error("unexpected end of file");
	  }
	  if(fread(short_buffer, psize, 1, fp)!=1)  {
	    cout << "failure to read" << endl;
	    error("unexpected end of file");
	  }
	}
      }  else {
	exception = true;	
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
	    *(m+lattice().local(idx_gl)*field_components+k)=short_buffer[k];
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
	local_index[buffer_size]=lattice().local(idx_gl);
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
  mpi.broadcast(reversed_header_endianess,processIO);
  if(try_swicth_endianess && reversed_header_endianess) {
    mpi << "swithing endiness...\n";
#ifdef USE_DOUBLE_PRECISION
    switch_endianess_8bytes();
#else
    switch_endianess_4bytes();
#endif
  }
  if(ME==0 && mdp_shutup==false) {
    printf("... Loading time: %f (sec)\n", mpi.time()-mytime);
    fflush(stdout);
  }
  return true;
}

