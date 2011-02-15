#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "string.h"
#include "time.h"

#define TRUE  1;
#define FALSE 0;

struct _generic_field_file_header {
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

int nx[4];

int main(int argc, char **argv) {

    printf("=======================================================\n");
    printf("Program for convertingMILC ascii gauge configurations\n"); 
    printf("and quark (fermi field) into MDP files\n");
    printf("=======================================================\n");

  if(argc<4) {
    printf("usage:\n\n");
    printf("asciimilc2mdp -gauge 16x08x08x08 input output\n\n");
    printf("asciimilc2mdp -quark 16x08x08x08 input output (to be tested) \n\n");
    exit(0);
  };

  sscanf(argv[2],"%ix%ix%ix%i", nx, nx+1, nx+2, nx+3);

  int x0, x1, x2, x3, mu, Ndim=4, a, i, j;
  long position;
  long time0= clock()/CLOCKS_PER_SEC;

  if(strcmp(argv[1],"-gauge")==0) {

    printf("Lattice: %i x %i x %i x %i\n", nx[0], nx[1], nx[2], nx[3]);
    printf("opening the MILC (ascii) file: %s (read)\n", argv[3]);
    FILE *TONY_fp=fopen(argv[3], "r");
    printf("opening the MDP file: %s (write) \n", argv[4]);
    FILE *MDP_fp =fopen(argv[4], "w");

    _generic_field_file_header myheader;

    myheader.ndim=4;
    int ii;
    for(ii=0; ii<4; ii++)  myheader.box_size[ii]=nx[ii];
    for(ii=4; ii<10; ii++) myheader.box_size[ii]=0;
    myheader.sites=nx[0]*nx[1]*nx[2]*nx[3];
    myheader.bytes_per_site=288; 
    myheader.endianess=0x87654321;
    strcpy(myheader.program_version, "Converted from MILC (ascii) file");
    strcpy(myheader.creation_date, "unknown");
    int offset=sizeof(_generic_field_file_header)/sizeof(char);
    fwrite(&myheader, offset, sizeof(char), MDP_fp);

    long a,b[4];
    char s[5][20];

    fscanf(TONY_fp, "%i%s%s%s%s%s%i%i%i%i",&a,s[0],s[1],s[2],s[3],s[4],b,b+1,b+2,b+3);

    printf("%i\n%s %s %s %s %s\n%i %i %i %i\n",a,s[0],s[1],s[2],s[3],s[4], b[0],b[1],b[2], b[3]);

    float buffer[18]; // this assumes data in single precision: 72 = 9 x 2 x 4
    for(x0=0; x0<nx[0]; x0++) 
      for(x3=0; x3<nx[3]; x3++) 
	for(x2=0; x2<nx[2]; x2++) 
	  for(x1=0; x1<nx[1]; x1++) 
	    for(mu=1; mu<=4; mu++) {
	      for(i=0; i<3; i++)
		for(j=0; j<3; j++) {
		  fscanf(TONY_fp, "%f%f", 
			 &(buffer[6*i+2*j]), 
			 &(buffer[6*i+2*j+1]));
		};

	      // this map to the MDP ordering
	      position=(((x0*nx[1]+x1)*nx[2]+x2)*nx[3]+x3)*4+(mu % Ndim);
	      fseek(MDP_fp, 72*position+offset, SEEK_SET);
	      fwrite(buffer, 72, 1, MDP_fp);
	    };
    fclose(TONY_fp);
    fclose(MDP_fp);
    printf("\nAll sites are OK.\n");
    printf("Done in %i secs.\n", clock()/CLOCKS_PER_SEC-time0);
  } else if(strcmp(argv[1],"-quark")==0) {

    printf("Lattice: %i x %i x %i x %i\n", nx[0], nx[1], nx[2], nx[3]);
    printf("opening the MILC (ascii) file: %s (read)\n", argv[3]);
    FILE *TONY_fp=fopen(argv[3], "r");
    printf("opening the MDP file: %s (write) \n", argv[4]);
    FILE *MDP_fp =fopen(argv[4], "w");

    _generic_field_file_header myheader;

    myheader.ndim=4;
    int ii;
    for(ii=0; ii<4; ii++)  myheader.box_size[ii]=nx[ii];
    for(ii=4; ii<10; ii++) myheader.box_size[ii]=0;
    myheader.sites=nx[0]*nx[1]*nx[2]*nx[3];
    myheader.bytes_per_site=96; 
    myheader.endianess=0x87654321;
    strcpy(myheader.program_version, "Converted from MILC (ascii) file");
    strcpy(myheader.creation_date, "unknown");
    int offset=sizeof(_generic_field_file_header)/sizeof(char);
    fwrite(&myheader, sizeof(char),offset, MDP_fp);

    float buffer[24];
    for(x0=0; x0<nx[0]; x0++) 
      for(x3=0; x3<nx[3]; x3++) 
	for(x2=0; x2<nx[2]; x2++) 
	  for(x1=0; x1<nx[1]; x1++) { 
	    for(a=0; a<4; a++) 
	      for(i=0; i<3; i++)
		fscanf(TONY_fp,"%f%f", &buffer[6*a+2*i], &buffer[6*a+2*i+1]);

	    // this map to the MDP ordering
	    position=(((x0*nx[1]+x1)*nx[2]+x2)*nx[3]+x3);
	    fseek(MDP_fp, 96*position+offset, SEEK_SET);
	    fwrite(buffer, 96, 1, MDP_fp);
	  };
    
    fclose(TONY_fp);
    fclose(MDP_fp);
    printf("\nAll sites are OK.\n");
    printf("Done in %i secs.\n", clock()/CLOCKS_PER_SEC-time0);
  } else if(strcmp(argv[1],"-gauge:d")==0) {
    
    printf("Lattice: %i x %i x %i x %i\n", nx[0], nx[1], nx[2], nx[3]);
    printf("opening the MILC (ascii) file: %s (read)\n", argv[3]);
    FILE *TONY_fp=fopen(argv[3], "r");
    printf("opening the MDP file: %s (write) \n", argv[4]);
    FILE *MDP_fp =fopen(argv[4], "w");
    
    _generic_field_file_header myheader;
    
    myheader.ndim=4;
    int ii;
    for(ii=0; ii<4; ii++)  myheader.box_size[ii]=nx[ii];
    for(ii=4; ii<10; ii++) myheader.box_size[ii]=0;
    myheader.sites=nx[0]*nx[1]*nx[2]*nx[3];
    myheader.bytes_per_site=576; 
    myheader.endianess=0x87654321;
    strcpy(myheader.program_version, "Converted from MILC (ascii) file");
    strcpy(myheader.creation_date, "unknown");
    int offset=sizeof(_generic_field_file_header)/sizeof(char);
    fwrite(&myheader, offset, sizeof(char), MDP_fp);
    
    double buffer[18]; //this assumes data in double precision: 144 = 9 x 2 x 8
    for(x0=0; x0<nx[0]; x0++) 
      for(x3=0; x3<nx[3]; x3++) 
	for(x2=0; x2<nx[2]; x2++) 
	  for(x1=0; x1<nx[1]; x1++) 
	    for(mu=0; mu<4; mu++) {
	      for(i=0; i<3; i++)
		for(j=0; j<3; j++) {
		  fscanf(TONY_fp, "%lf%lf", 
			 &(buffer[6*j+2*i]), 
			 &(buffer[6*j+2*i+1]));
		  buffer[6*j+2*i+1]*=-1;
		};
	      // this map to the MDP ordering
	      position=(((x0*nx[1]+x1)*nx[2]+x2)*nx[3]+x3)*4+((mu+1) % Ndim);
	      fseek(MDP_fp, 144*position+offset, SEEK_SET);
	      fwrite(buffer, 144, 1, MDP_fp);
	    };
    fclose(TONY_fp);
    fclose(MDP_fp);
    printf("\nAll sites are OK.\n");
    printf("Done in %i secs.\n", clock()/CLOCKS_PER_SEC-time0);
  } else if(strcmp(argv[1],"-quark:d")==0) {
    
    printf("Lattice: %i x %i x %i x %i\n", nx[0], nx[1], nx[2], nx[3]);
    printf("opening the MILC (ascii) file: %s (read)\n", argv[3]);
    FILE *TONY_fp=fopen(argv[3], "r");
    printf("opening the MDP file: %s (write) \n", argv[4]);
    FILE *MDP_fp =fopen(argv[4], "w");

    _generic_field_file_header myheader;

    myheader.ndim=4;
    int ii;
    for(ii=0; ii<4; ii++)  myheader.box_size[ii]=nx[ii];
    for(ii=4; ii<10; ii++) myheader.box_size[ii]=0;
    myheader.sites=nx[0]*nx[1]*nx[2]*nx[3];
    myheader.bytes_per_site=192; 
    myheader.endianess=0x87654321;
    strcpy(myheader.program_version, "Converted from MILC (ascii) file");
    strcpy(myheader.creation_date, "unknown");
    int offset=sizeof(_generic_field_file_header)/sizeof(char);
    fwrite(&myheader, sizeof(char),offset, MDP_fp);
    
    double buffer[24];
    for(x0=0; x0<nx[0]; x0++) 
      for(x3=0; x3<nx[3]; x3++) 
	for(x2=0; x2<nx[2]; x2++) 
	  for(x1=0; x1<nx[1]; x1++) { 
	    for(a=0; a<4; a++) 
	      for(i=0; i<3; i++)
		fscanf(TONY_fp,"%lf%lf", &buffer[6*a+2*i], &buffer[6*a+2*i+1]);
	    // this map to the MDP ordering
	    position=(((x0*nx[1]+x1)*nx[2]+x2)*nx[3]+x3);
	    fseek(MDP_fp, 192*position+offset, SEEK_SET);
	    fwrite(buffer, 192, 1, MDP_fp);
	  };
    
    fclose(TONY_fp);
    fclose(MDP_fp);
    printf("\nAll sites are OK.\n");
    printf("Done in %i secs.\n", clock()/CLOCKS_PER_SEC-time0);
  };
};












