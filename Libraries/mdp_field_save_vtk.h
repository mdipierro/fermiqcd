/////////////////////////////////////////////////////////////////
/// @file mdp_field_save.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains file IO operations for class mdp_field
///
/// Licensed under GPL2 license
/// Read attached license in file mdp_license.pdf
/// This file cannot be distributed without file mdp_license.pdf
//////////////////////////////////////////////////////////////////

/**
 * Best way to save a field to a VTK file.
 * @param filename
 *   String giving name of file to create
 * @param t
 *   Int that specifies which timeslice to save.  Value of -1 saves all timeslices.
 * @param component
 *   Int specifying which component to save.  Value of -1 saves all components.
 * @processIO
 *   Int naming the ID of processor that should perform the save
 * @param ASCII
 *   Bool flag to save data in ASCII when true, binary when false
 * @return
 *   Returns true on success, throws an error for unexpected conditions.
 */
template<class T>
bool mdp_field<T>::save_vtk(string filename, 
			    int t,
			    int component,
			    int processIO,
			    bool ASCII) {

  filename=next_to_latest_file(filename);
  string filename_tmp=filename+".tmp";

  int max_buffer_size=1024;
  int timeslice;
  mdp_int header_size=0;
  mdp_int psize=field_components*Tsize;
  mdp_int idx_gl, nvol_gl=lattice().nvol_gl, k;
  double mytime=mpi.time();
  header.reset();
  if(lattice().ndim<3 || lattice().ndim>4) error("mdp_field::save_vtk only works for ndim=3 and 4");
  if(ME==processIO) {
    mdp_int *buffer_size=new mdp_int[Nproc];
    mdp_int *buffer_ptr =new mdp_int[Nproc];
    mdp_array<T,3> large_buffer(Nproc,max_buffer_size,field_components);
    T *short_buffer=new T[field_components];
    int process;
    for(process=0; process<Nproc; process++) buffer_ptr[process]=0;
    cout << "Saving file " << filename 
	 << " from process " << processIO 
	 << " (buffer = " << max_buffer_size << " sites)" << '\n'; 
    fflush(stdout);
    FILE *fp=0;    

    float fval;
    char tmp[1024];
    char header[1024];
    int skip_bytes=0;
        
    int space_volume=lattice().size()/lattice().size(0);
    int LZ=lattice().size(1);
    int LY=lattice().size(2);
    int LX=lattice().size(3);    
    if(lattice().ndim==3) {
      space_volume=lattice().size();
      LZ=lattice().size(0);
      LY=lattice().size(1);
      LX=lattice().size(2);
      t=-1;
    }
    int fc,fc2;
    if(component==-1) fc2=fc=field_components; else fc2=fc=component;

    sprintf(header,
	    "# vtk DataFile Version 2.0\n"
	    "%s\n"
	    "%s\n"
	    "DATASET STRUCTURED_POINTS\n"
	    "DIMENSIONS %i %i %i\n"
	    "ORIGIN     0   0   0\n"
	    "SPACING    1   1   1\n"
	    "POINT_DATA %i",
	    filename.c_str(),
	    ((ASCII)?"ASCII":"BINARY"),
	    LX,LY,LZ,LX*LY*LZ);

    fp=fopen(filename_tmp.c_str(),"wb");
    fwrite(header,sizeof(char),strlen(header),fp);
   
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
	  short_buffer[k]=*(m+lattice().local(idx_gl)*field_components+k);
      }
      if(process!=NOWHERE) {
	timeslice=idx_gl/space_volume;
	if(idx_gl % space_volume==0 && (t<0 || timeslice==t)) {	  
	  sprintf(header,"\nSCALARS scalars_t%i float\nLOOKUP_TABLE default\n",timeslice);
	  fwrite(header,sizeof(char),strlen(header),fp);
	}
	if(t<0 || timeslice==t || lattice().ndim==3) {	  
	  for(fc=0; fc<field_components; fc++) if(component==-1 || fc==component) {
	    fval=(float) short_buffer[fc];
	    if(!ASCII) {
	      switch_endianess_byte4(fval);
	      if(fwrite(&fval,sizeof(fval),1,fp)!=1) 
		error("probably out of disk space");
	    } else {
	      sprintf(tmp,"%e ",fval);
	      if(fwrite(tmp,strlen(tmp), 1,fp)!=1)
		error("probably out of disk space");
	    }
	  }
	}
      }
    }
    if(fp) {
      fclose(fp);
      remove(filename.c_str());
      rename(filename_tmp.c_str(),filename.c_str());
    }
    delete[] buffer_size;
    delete[] buffer_ptr;
    delete[] short_buffer;
  } else {
    int  process;
    mdp_int buffer_size=0, idx, idx_gl;
    mdp_int *local_index=new mdp_int[max_buffer_size];
    mdp_array<T,2> local_buffer(max_buffer_size,field_components);
    mdp_request request;
    for(idx_gl=0; idx_gl<nvol_gl; idx_gl++) {
      process=where_global(idx_gl);
      if(process==ME) {
	local_index[buffer_size]=lattice().local(idx_gl);
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
  return true;
}

mdp_field<float>& cumulate_field(mdp_field<float>& field, string filename) {
  static map<string,int> counter;
  static map<string,mdp_field<float>*> fields;
  mdp_site p(field.lattice());
  int k=0;
  if(counter.find(filename)==counter.end()) {
    fields[filename]=new mdp_field<float>(field.lattice(),field.size_per_site());
    forallsites(p) for(int i=0; i<field.size_per_site(); i++) (*fields[filename])(p,i)=field(p,i);
    counter[filename]=1;
  } else {
    k=counter[filename];
    forallsites(p) for(int i=0; i<field.size_per_site(); i++)
      (*fields[filename])(p,i)=(k*(*fields[filename])(p,i)+field(p,i))/(k+1);
    counter[filename]++;
  }
  return (*fields[filename]);
}

/*
void save_vtk(mdp_field<float> &field,
	      string filename, 
	      int t=-1,
	      int component=-1,
	      int processIO=0,
	      bool ASCII=false) {
  static map<string,int> counter;
  static map<string,mdp_field<float>*> fields;
  mdp_site p(field.lattice());
  int k=0;
  if(counter.find(filename)==counter.end()) {
    fields[filename]=new mdp_field<float>(field.lattice(),field.size_per_site());
    forallsites(p) for(int i=0; i<field.size_per_site(); i++) (*fields[filename])(p,i)=field(p,i);
    counter[filename]=1;    
  } else {
    k=counter[filename];
    forallsites(p) for(int i=0; i<field.size_per_site(); i++) 
      (*fields[filename])(p,i)=(k*(*fields[filename])(p,i)+field(p,i))/(k+1);
    counter[filename]++;
  }
  if(filename[filename.size()-1]!='*') {
    filename=filename.substr(0,filename.size()-4);
    field.save_vtk(filename+".vtk",t,component,processIO,ASCII);
    fields[filename]->save_vtk(filename+"_average.vtk",t,component,processIO,ASCII);
  } else {
    string filename2=filename.substr(0,filename.size()-1);
    field.save_vtk(filename2+"_"+tostring(k)+".vtk",t,component,processIO,ASCII);
    fields[filename]->save_vtk(filename2+"_average_"+tostring(k)+".vtk",t,component,processIO,ASCII);
  }
}
*/
