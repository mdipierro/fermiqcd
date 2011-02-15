
void ildg_gauge_reader(gauge_field& U, 
		       filename, 
		       header_bytes=0, 
		       precision=16) {
  
  mdp_site x(lattice());
  mdp_site y(lattice());
  mdp_matrix v[4];
  if(precision==16) U.load_as_float(filename,0,512,false,header_bytes);
  else if(precision==32) U.load_as_double(filename,0,512,false,header_bytes);
  forallsites(x) {
    if(x(1)<x(3)) {
      y.set(x(0),x(3),x(2),x(1));
      v[0]=U(x,3);
      v[1]=U(x,0);
      v[2]=U(x,1);
      v[3]=U(x,2);
      U(x,0)=U(y,3);
      U(x,1)=U(y,0);
      U(x,2)=U(y,1);
      U(x,3)=U(y,2);
      U(y,0)=v[0];
      U(y,1)=v[1];
      U(y,2)=v[2];
      U(y,3)=v[3];      
    } 
  }
}
