#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "string.h"
#include "time.h"
#include "complex.h"


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

#ifndef USE_DOUBLE_PRECISION
typedef complex<float> Complex;
#else
typedef complex<double> Complex;
#endif

class _generic_field_file_header {
 public:
  char file_id[60];
  char program_version[60];
  char creation_date[60];
  unsigned long endianess;
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
  int size;
  int  dim[7];
  short_field() {
    m=0;
  };
  void initialize(int x1, int x2, int x3, int a=1, int b=1, int c=1, int d=1) {
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

int main(int argc, char **argv) {

  printf("======================================================\n");
  printf("Program for converting LUIGI gauge configurations and\n"); 
  printf("and propagators into MDP files\n");
  printf("Conversion of propagators not implemented yet\n");
  printf("======================================================\n");
  
  int nx[4], nc;
    
    
  
  if(argc<4) {
    printf("usage:\n\n");
    printf("luigi2mdp -gauge 16x08x08x08,NC input output\n\n");
    exit(0);
  };
  sscanf(argv[2],"%ix%ix%ix%i,%i", nx, nx+1, nx+2, nx+3, &nc);
  
  int x0, x1, x2, x3, mu, nu;
  long time0= clock()/CLOCKS_PER_SEC;

  if(strcmp(argv[1],"-gauge")==0) {
    
    printf("Lattice: %i x %i x %i x %i\n", nx[0], nx[1], nx[2], nx[3]);
    printf("opening the LUIGI file: %s (read)\n", argv[3]);
    FILE *LUIGI_fp=fopen(argv[3], "r");
    printf("opening the MDP file: %s (write) \n", argv[4]);
    FILE *MDP_fp =fopen(argv[4], "w");
    
    _generic_field_file_header myheader;
    
    myheader.ndim=4;
    int ii;
    for(ii=0; ii<4; ii++)  myheader.box_size[ii]=nx[ii];
    for(ii=4; ii<10; ii++) myheader.box_size[ii]=0;
    myheader.sites=nx[0]*nx[1]*nx[2]*nx[3];
    if(strcmp(argv[1],"-gauge")==0) 
      myheader.bytes_per_site=4*nc*nc*sizeof(Complex); 
    myheader.endianess=0x87654321;
    strcpy(myheader.program_version, "Converted from LUIGI");
    time_t tt;
    time(&tt);
    strcpy(myheader.creation_date, ctime(&tt));
    int offset=sizeof(_generic_field_file_header)/sizeof(char);
    fwrite(&myheader, sizeof(char),offset, MDP_fp);

    short_field U;
    U.initialize(nx[1], nx[2], nx[3], 4, nc, nc);

    int matrix_size=nc*nc*sizeof(Complex);

    for(x0=0; x0<nx[0]; x0++) {
      for(x3=0; x3<nx[3]; x3++) 
	for(x2=0; x2<nx[2]; x2++) 
	  for(x1=0; x1<nx[1]; x1++) 
	    for(mu=0; mu<4; mu++) {
	      switch(mu) {
		case 0: nu=0; break;
		case 1: nu=3; break;
		case 2: nu=2; break;
		case 3: nu=1; break;
	      }
	      fread(&U(x1,x2,x3, nu,0,0), matrix_size, 1, LUIGI_fp);
	    }
      
      fwrite(U.m, U.size, sizeof(Complex), MDP_fp);
    };
    fclose(LUIGI_fp);
    fclose(MDP_fp);
    printf("\nAll sites are OK.\n");
    printf("Done in %li secs.\n", clock()/CLOCKS_PER_SEC-time0);
  } else {
    printf("I am sorry conversion of propagators not implemnetd yet!\n");
  };
};












