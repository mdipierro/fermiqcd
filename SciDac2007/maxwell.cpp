#include "mdp.h"
#include "mdp_all.h"
#define X 0
#define Y 1
#define Z 2

void dump(mdp_field<int>& s, string filename="default.vtk") {
  char header[1024], number[1024];
  int LX=s.lattice().size(0), LY=s.lattice().size(1), LZ=s.lattice().size(2);
  sprintf(header,
	  "# vtk DataFile Version 2.0\n"
	  "Really cool data\n"
	  "ASCII\n"
	  "DATASET STRUCTURED_POINTS\n"
	  "DIMENSIONS %i %i %i\n"
	  "ORIGIN     0   0   0\n"
	  "SPACING    1   1   1\n"
	  "POINT_DATA %i\n"
	  "SCALARS scalars0 float 1\n"
	  "LOOKUP_TABLE default\n\n",LX,LY,LZ,LX*LY*LZ);
  
  int sfd=open("tmp.vtk",O_WRONLY);
  cout << "saving... " << filename << " as fd=" << sfd << endl;
  setFileLock(sfd);
  write(sfd,header,strlen(header));
  site p(s.lattice());
  for(int i=0; i<LX; i++)
    for(int j=0; j<LY; j++)
      for(int k=0; k<LZ; k++) {
	p.set(i,j,k);
	sprintf(number,"%i\n", (int) s(p));
	write(sfd, number, strlen(number));
      }
  setFileUnlock(sfd);
  system((string("cp tmp.vtk ")+filename).c_str());
  close(sfd);
}

int main(int argc, char** argv) {
  mdp.open_wormholes(argc,argv);
  int L[]={20,20,20};
  mdp_lattice cube(3,L);
  complex_field E(cube,3);
  complex_field B(cube,3);
  complex_field q(cube);
  complex_field j(cube,3);

  forallsites(x) {
    E(x)=;;
    B(x)=;
  }

  mdp.close_wormholes();
  return 1;
}
