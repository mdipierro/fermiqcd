// Program: application1.cpp
#include "mdp.h" 

int main(int argc, char **argv) { 
  mdp.open_wormholes(argc,argv);
  int mybox[]={20,20,20}; 
  mdp_lattice vacuum(3,mybox,
		     default_partitioning<0>,
		     box_topology);
  mdp_field<float> u(vacuum); 
  mdp_field<float> q(vacuum);
  mdp_site x(vacuum);
  mdp_site A(vacuum);
  mdp_site B(vacuum);
  float precision, old_u;
  forallsitesandcopies(x) { 
    u(x)=0;
    if(x(0)==3 && x(1)==3 && x(2)==3) q(x)=3;
    else if(x(0)==17 && x(1)==17 && x(2)==17) q(x)=3;
    else          q(x)=0;
  }
  do {
    precision=0;
    forallsites(x) if((x(0)>0) && (x(0)<mybox[0]-1) &&
		      (x(1)>0) && (x(1)<mybox[1]-1) &&
		      (x(2)>0) && (x(2)<mybox[2]-1)) {
      old_u=u(x);
      u(x)=(q(x)+u(x+0)+u(x-0)+u(x+1)+u(x-1)+u(x+2)+u(x-2))/6;
      precision+=pow(u(x)-old_u,2);
    }
    u.update();
    mdp.add(precision);
  } while (sqrt(precision)>0.0001);
  u.save("potential.dat");
  mdp.close_wormholes(); 
}
