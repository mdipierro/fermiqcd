#include "mdp.h" 

int ruleofgame(int a00, int a01, int a02,
	       int a10, int a11, int a12,
	       int a20, int a21, int a22) {
  int sum=a00+a01+a02+a10+a12+a20+a21+a22;
  if(a11==0 && sum==3) return 1;
  else if(a11==1 && (sum==2 || sum==3)) return 1;
  return 0;
}

int myownpartitioning(int *x, 
		      int ndim, 
		      int *nx) {
  if(x[0]==2 && x[1]==3) return NOWHERE;
  return (x[0]+x[1]) % Nproc;
}

void board_with_hole(int mu, 
		    int *x_dw, 
		    int *x, 
		    int *x_up, 
		    int ndim, 
		    int *nx) {
  for(int nu=0; nu<ndim; nu++) if(nu==mu) {
    x_dw[mu]=(x[mu]-1+nx[mu]) % nx[mu];
    x_up[mu]=(x[mu]+1) % nx[mu];
  } else x_up[nu]=x_dw[nu]=x[nu];
  if(x_up[0]==2 && x_up[1]==3) {
    x_up[mu]=(x_up[mu]+1) % nx[mu];
  }
  if(x_dw[0]==2 && x_dw[1]==3) {
    x_dw[mu]=(x_dw[mu]-1+nx[mu]) % nx[mu];
  }
}


int main(int argc, char **argv) {
  mdp.open_wormholes(argc,argv);
  int sides[]={10,10};
  mdp_lattice board(2,sides,myownpartitioning,board_with_hole);
  mdp_field<int> S(board);
  mdp_field<int> newS(board);
  mdp_site x(board);

  forallsites(x) S(x)=(board.random(x).plain()>0.5)?1:0;

  int sum;
  for(int k=0; k<1000; k++) {
    S.update();
    forallsites(x)
      newS(x)=ruleofgame(S((x-1)-0),S(x-1),S((x-1)+0),
			 S(x-0),    S(x),  S(x+0),
			 S((x+1)-0),S(x+1),S((x+1)+0));    
    sum=0;
    forallsites(x) {
      S(x)=newS(x);
      sum=sum+S(x);
    }
    mdp.add(sum);
    mdp << sum << endl;
    if(sum==0) break;
  }
  

  mdp.close_wormholes();
  return 0;
}
