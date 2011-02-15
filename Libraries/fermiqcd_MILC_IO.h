/////////////////////////////////////////////////////////////////
/// @file fermiqcd_MILC_IO.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Functions to read a MILC gauge configuration without conversion
///
/// Distributed under GPL2 License
/// 
/// Created with support from the US Department of Energy
//////////////////////////////////////////////////////////////////


bool milc_read_as_float_noswitch(FILE *fp,
				 void* data,
				 mdp_int psize,
				 mdp_int header_size,
				 mdp_int position,
				 const mdp_lattice &lattice) {
  cout << position << '\n';
  double *p=(double*) data;
  site x(lattice);
  x.set_global(position);
#ifdef USE_DOUBLE_PRECISION
  float  *q=(float*) malloc(psize/2);
  position=(((x(0)*lattice.size(3)+x(3))*lattice.size(2)+x(2))*lattice.size(1)+x(1));
  if(fseek(fp, position*psize/2+header_size, SEEK_SET) ||
     fread(q, psize/2,1, fp)!=1) {
    return false;
  }
  for(uint i=0; i<psize/sizeof(double); i++) {
    p[i]=q[i];
  }
  free(q);
#else 
  float  *q=(float*) malloc(psize);
  position=(((x(0)*lattice.size(3)+x(3))*lattice.size(2)+x(2))*lattice.size(1)+x(1));
  if(fseek(fp, position*psize+header_size, SEEK_SET) ||
     fread(q, psize,1, fp)!=1) {
    return false;
  }
  for(uint i=0; i<psize/sizeof(float); i++) {
    p[i]=q[i];
  }
  free(q);
#endif
  return true;
}

bool milc_read_as_float_switch(FILE *fp,
			       void* data,
			       mdp_int psize,
			       mdp_int header_size,
			       mdp_int position,
			       const mdp_lattice &lattice) {
  cout << "switch" << position << '\n';
  double *p=(double*) data;
  site x(lattice);
  x.set_global(position);
#ifdef USE_DOUBLE_PRECISION
  float  *q=(float*) malloc(psize/2);
  position=(((x(0)*lattice.size(3)+x(3))*lattice.size(2)+x(2))*lattice.size(1)+x(1));
  if(fseek(fp, position*psize/2+header_size, SEEK_SET) ||
     fread(q, psize/2,1, fp)!=1) return false;
  for(uint i=0; i<psize/sizeof(double); i++) {
    switch_endianess_byte4(q[i]);
    p[i]=q[i];
  }
  free(q);
#else 
  float  *q=(float*) malloc(psize);
  position=(((x(0)*lattice.size(3)+x(3))*lattice.size(2)+x(2))*lattice.size(1)+x(1));

  if(fseek(fp, position*psize+header_size, SEEK_SET) ||
     fread(q, psize,1, fp)!=1) {
       return false;
  }
  for(uint i=0; i<psize/sizeof(float); i++) {
    switch_endianess_byte4(q[i]);
    p[i]=q[i];
  }
  free(q);
#endif
  return true;
}


bool load_milc(gauge_field &U, string filename, 
	       mdp_int max_buffer_size=128, int processIO=0) {

  struct {
    mdp_int magic_number;               /* Identifies file format */
    mdp_int dims[4];                    /* Full lattice dimensions */
    char   time_stamp[64];           /* Date and time stamp - used to
					check consistency between the
					ASCII header file and the
					lattice file */
    mdp_int header_bytes;               /* NOT WRITTEN TO THE FILE but
					helpful for finding the data */
    mdp_int what_is_this;
    mdp_int order;                      /* 0 means no coordinate list is
					attached and the values are in
					coordinate serial order.
					Nonzero means that a
					coordinate list is attached,
					specifying the order of values */
    // mdp_int boo;
  } milc_header;

  bool ew=false;
  int size=sizeof(milc_header);
  FILE *fp=fopen(filename.c_str(),"r");
  if(fp==0) return false;
  if(fread(&milc_header, 1, size, fp)!=size) {
    fclose(fp);
    return false;
  }

  if(milc_header.magic_number==0x874e0000) {
    switch_endianess_byte4(milc_header.dims[0]);
    switch_endianess_byte4(milc_header.dims[1]);
    switch_endianess_byte4(milc_header.dims[2]);
    switch_endianess_byte4(milc_header.dims[3]);
    ew=true;
  }

  if(U.lattice().ndim!=4 ||
     milc_header.dims[0]!=U.lattice().size(1) ||
     milc_header.dims[1]!=U.lattice().size(2) ||
     milc_header.dims[2]!=U.lattice().size(3) ||
     milc_header.dims[3]!=U.lattice().size(0)) {
    fclose(fp);
    return false;
  }
  fclose(fp);

  if(ew)
    return U.load(filename,processIO, max_buffer_size,false,size,
		  milc_read_as_float_switch, false);
  else
    return U.load(filename,processIO, max_buffer_size,false,size,
		  milc_read_as_float_noswitch, false);
}

