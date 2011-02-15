#include"stdio.h"
#include"math.h"
#include "complex.h"
#include "time.h"

template<class T>
void switch_endianess_byte4(T &a) {
  char *p=(char*) &a;
  static char q[4];
  if(sizeof(T)==4) {
    q[0]=p[0];
    q[1]=p[1];
    q[2]=p[2];
    q[3]=p[3];
    p[0]=q[3];
    p[1]=q[2];
    p[2]=q[1];
    p[3]=q[0];
  } else printf("error endianess\n");
}

#define Complex complex<float>

#define TRUE  1;
#define FALSE 0;

#ifndef _type32
#define _type32
#ifdef SHORT32
typedef short type32;
#else
typedef int type32;
#endif
#endif

/**********************************************************************/
/* Binary lattice formats                                             */
/**********************************************************************/
/* In version 5 we have two binary lattice file formats:

      serial files     Data written in coordinate natural order
      checkpoint files Data written in node dump order

      Further descriptive information is kept in a separate ASCII header
      file.  See below.

      */


/*--------------------------------------------------------------------*/
/* version 5 binary file format */

#define GAUGE_VERSION_NUMBER 20103
#define MAX_TIME_STAMP 64

/* 1. Header comes first    */

typedef struct {
  type32 magic_number;               /* Identifies file format */
  type32 dims[4];                    /* Full lattice dimensions */
  char   time_stamp[MAX_TIME_STAMP]; 
                                     /* Date and time stamp - used to
					check consistency between the
					ASCII header file and the
					lattice file */
  type32 header_bytes;               /* NOT WRITTEN TO THE FILE but
					 helpful for finding the data */
  type32 what_is_this;
  type32 order;                      /* 0 means no coordinate list is
				        attached and the values are in
				        coordinate serial order.
				        Nonzero means that a
				        coordinate list is attached,
				        specifying the order of values */
  type32 boo;
} gauge_header;


class _generic_field_file_header {
 public:
  char file_id[60];
  char program_version[60];
  char creation_date[60];
  float endianess;
  int  ndim;
  int  box_size[10];
  long bytes_per_site;
  long sites;
  _generic_field_file_header() {
    strcpy(file_id, "File Type: MDP FIELD\n");
  };
};

int number(char *x) {
  return 10*(((int) x[0])-48) + (((int) x[1]) -48); 
};

class short_field {
public:
  Complex *m;
  long size;
  int  dim[7];
  short_field() {
    m=0;
  };
  int initialize(int x1, int x2, int x3, int a=1, int b=1, int c=1, int d=1) {
    size=x1*x2*x3*a*b*c*d;
    dim[0]=x1;
    dim[1]=x2;
    dim[2]=x3;
    dim[3]=a;
    dim[4]=b;
    dim[5]=c;
    dim[6]=d;
    if(m!=0) delete[] m;
    m=new Complex[size];
  };
  Complex &operator() (int x1, int x2, int x3, int a=0, int b=0, int c=0, int d=0) {
    return m[(((((x1*dim[1]+x2)*dim[2]+x3)*dim[3]+a)*dim[4]+b)*dim[5]+c)*dim[6]+d];
  };
};

int nx[4];

int main(int argc, char **argv) {

    printf("======================================================\n");
    printf("Program for converting MILC gauge configurations and\n"); 
    printf("and propagators into MDP files\n");
    printf("Conversion of propagators not implemented yet\n");
    printf("======================================================\n");

  if(argc<4) {
    printf("usage:\n\n");
    printf("milc2mdp -gauge 16x08x08x08 input output\n\n");
    printf("milc2mdp -fermi 16x08x08x08 input output\n\n");
    exit(0);
  };
  sscanf(argv[2],"%ix%ix%ix%i", nx, nx+1, nx+2, nx+3);

  int x0, x1, x2, x3, mu, Ndim=4;
  long position;
  long time0= clock()/CLOCKS_PER_SEC;

  if(strcmp(argv[1],"-gauge")==0) {

    printf("Lattice: %i x %i x %i x %i\n", nx[0], nx[1], nx[2], nx[3]);
    printf("opening the MILC file: %s (read)\n", argv[3]);
    FILE *MILC_fp=fopen(argv[3], "r");
    printf("opening the MDP file: %s (write) \n", argv[4]);
    FILE *MDP_fp =fopen(argv[4], "w");

    gauge_header               MILC_header;
    _generic_field_file_header myheader;
    fread(&MILC_header, sizeof(gauge_header)/sizeof(char),1,MILC_fp);

    printf("%x\n", MILC_header.magic_number);

    if(MILC_header.magic_number==0x874e0000) {
      printf("switching endianess\n");
      switch_endianess_byte4(MILC_header.dims[0]);
      switch_endianess_byte4(MILC_header.dims[1]);
      switch_endianess_byte4(MILC_header.dims[2]);
      switch_endianess_byte4(MILC_header.dims[3]);
    }

    if((MILC_header.dims[0]!=nx[1]) &&
       (MILC_header.dims[1]!=nx[2]) &&
       (MILC_header.dims[2]!=nx[3]) &&
       (MILC_header.dims[3]!=nx[0])) {

      printf("%xx%xx%xx%x => Incorrect dimensions!\n",
	     MILC_header.dims[3],
	     MILC_header.dims[0],
	     MILC_header.dims[1],
	     MILC_header.dims[2]);

      exit(1);
    };

    myheader.ndim=4;
    int ii;
    for(ii=0; ii<4; ii++)  myheader.box_size[ii]=nx[ii];
    for(ii=4; ii<10; ii++) myheader.box_size[ii]=0;
    myheader.sites=nx[0]*nx[1]*nx[2]*nx[3];
    if(strcmp(argv[1],"-gauge")==0) myheader.bytes_per_site=288; 
    if(strcmp(argv[1],"-fermi")==0) myheader.bytes_per_site=1152; 
    myheader.endianess=0x87654321;
    strcpy(myheader.program_version, "Converted from MILC");
    strcpy(myheader.creation_date, MILC_header.time_stamp);
    int offset=sizeof(_generic_field_file_header)/sizeof(char);
    fwrite(&myheader, sizeof(char),offset, MDP_fp);

    short_field U;
    U.initialize(nx[1], nx[2], nx[3], 4, 3, 3);

    int matrix_size=72;
    char buffer[matrix_size];  // assumes single precision: 72 = 9 x 2 x 4
    for(x0=0; x0<nx[0]; x0++) {
      for(x3=0; x3<nx[3]; x3++) 
	for(x2=0; x2<nx[2]; x2++) 
	  for(x1=0; x1<nx[1]; x1++) 
	    for(mu=1; mu<=4; mu++) 
	      fread(&U(x1,x2,x3, mu % 4,0,0), matrix_size, 1, MILC_fp);
      
    if(MILC_header.magic_number==0x874e0000) 
      for(mu=0; mu<U.size; mu++) {
	switch_endianess_byte4(*((long*) U.m+2*mu));
	switch_endianess_byte4(*((long*) U.m+2*mu+1));
      }

      fwrite(U.m,U.size, sizeof(Complex), MDP_fp);
    };
    fclose(MILC_fp);
    fclose(MDP_fp);
    printf("\nAll sites are OK.\n");
    printf("Done in %i secs.\n", clock()/CLOCKS_PER_SEC-time0);
  } else {
    printf("I am sorry conversion of propagators not implemnetd yet!\n");
  };
};












