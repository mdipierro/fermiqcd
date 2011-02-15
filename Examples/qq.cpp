#include "fermiqcd.h"       // include FermiQCD libraries

// Usage example:
// ./a.out gauge_filename 2 4
// 
// computes average 2x4 loop looping over all space directions

int main(int argc, char** argv) {
  mdp.open_wormholes(argc,argv);             // START
  string gauge_filename=argv[1];
  int size1=atoi(argv[2]);
  int size2=atoi(argv[3]);
  mdp_field_file_header header;
  double result=0.0;

  // read file metadata
  if(is_file(gauge_filename)) header=get_info(gauge_filename);
  else error("Unable to access gauge configuration\n");
  if(header.ndim!=4) error("sorry, only in 4D");
  int nc=(int) sqrt((double) header.bytes_per_site/(4*sizeof(mdp_complex)));
  int *L=header.box; // lattice size

  // create lattice and read it in
  mdp_lattice lattice(4,L); // make a 4D lattice
  gauge_field U(lattice,nc); // make a gauge field U
  U.load(gauge_filename);

  
  int length=2*size1+2*size2;
  int path[length][2];

  // make a generic path
  for(int i=0; i<size1; i++) {
    path[i][0]=+1;
    path[i+size1+size2][0]=-1;
  }
  for(int i=size1; i<size1+size2; i++) {
    path[i][0]=+1;
    path[i+size1+size2][0]=-1;
  }
  // loop over all possible paths
  for(int mu=1; mu<4; mu++)
    for(int nu=mu+1; nu<4; nu++) {
      // build each path
      for(int i=0;i<size1;i++) path[i][1]=path[i+size1+size2][1]=mu;
      for(int i=size1;i<size1+size2;i++) path[i][1]=path[i+size1+size2][1]=nu;

      result+=real(average_path(U,length,path))/6;
    }
  cout << "average loop " << size1 << "x" << size2 << " = " << result << endl;
  
  mdp.close_wormholes();                     // STOP
  return 0;
}
