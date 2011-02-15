#define X 0
#define Y 1
#define Z 2

void dump(mdp_field<float>& s,
	  int site_idx=0,
	  string filename="default.vtk",
	  bool bASCII=true){
  static const char *tempfile="tmp.vtk";
  char header[1024];
  char number[1024];
  char vectorHeader[512]; 
  int  LX=s.lattice().size(0), 
       LY=s.lattice().size(1), 
       LZ=s.lattice().size(2);
  FILE *file=NULL;
  
  sprintf(header,
	  "# vtk DataFile Version 2.0\n"
	  "Really cool data\n"
	  "%s\n"
	  "DATASET STRUCTURED_POINTS\n"
	  "DIMENSIONS %i %i %i\n"
	  "ORIGIN 0 0 0\n"
	  "SPACING 1 1 1\n"
	  "POINT_DATA %i\n"
	  "SCALARS scalar0 float 1\n"
  	  "LOOKUP_TABLE default\n",
	  ((bASCII)?"ASCII":"BINARY"),
	  LX+1,LY+1,LZ+1,
	  (LX+1)*(LY+1)*(LZ+1));

  remove(tempfile); //remove tmp.vtk file if exists
  file=fopen(tempfile, "wb");
  int sfd=fileno(file);
  cout << "saving... " << filename << " as fd=" << sfd << endl;
  fwrite(header, sizeof(char), strlen(header),file);
  site p(s.lattice());
  
  float fval=0.0F; 
  for(int i=0; i<LX+1; i++)
    for(int j=0; j<LY+1; j++)
      for(int k=0; k<LZ+1; k++){	
        p.set(i % LX,j % LY,k % LZ);
        fval=(float)s(p,site_idx);
	if(bASCII){
          memset(number, 0, sizeof(number));
          sprintf(number,"%e\n",fval);
	  fwrite(number, sizeof(char), strlen(number), file);
        }
        else{	
          fwrite(&fval, sizeof(float), 1, file);
        }
      }
 
  fclose(file);
  file=NULL;
  remove(filename.c_str());
  rename(tempfile, filename.c_str());
  return;
}
