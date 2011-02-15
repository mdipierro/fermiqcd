#include "mdp.h" 

int ruleofgame(int a00, int a01, int a02,
	       int a10, int a11, int a12,
	       int a20, int a21, int a22) {
  int sum=a00+a01+a02+a10+a12+a20+a21+a22;
  // if dead (0) and 3 neighbors alive => alive (1)
  if(a11==0 && sum==3) return 1;
  // if alive (1) and 2 or 3 neighbors alive => alive (1)
  else if(a11==1 && (sum==2 || sum==3)) return 1;
  // die!
  return 0;
}

// START OPTIONAL

// parallelize by diagonal stripes and punch hole in board
int myownpartitioning(int *x, 
		      int ndim, 
		      int *nx) {
  if(x[0]==2 && x[1]==3) return NOWHERE; 
  return (x[0]+x[1]) % Nproc;
}

// make topolgy with a hole
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

// END OPTION

int main(int argc, char **argv) {
  mdp.open_wormholes(argc,argv); // START

  int sides[]={10,10};
  // a regular board parallelized by columns
  mdp_lattice board(2,sides);
  // or a board with a hole parallelized by diagonal stripes
  // mdp_lattice board(2,sides,myownpartitioning,board_with_hole);
  mdp_field<int> S(board);    // create field of int on board (S)
  mdp_field<int> newS(board); // create field of int in board (newS)
  mdp_site x(board);          // create variable to loop on board

  // initialize board at random
  forallsites(x) S(x)=(board.random(x).plain()>0.5)?1:0;

  int sum,diff;
  // iterate 1000 times max
  for(int k=0; k<1000; k++) {    
    S.update(); // communicate
    forallsites(x) // loop over all local sites and apply rules
      newS(x)=ruleofgame(S((x-1)-0),S(x-1),S((x-1)+0),
			 S(x-0),    S(x),  S(x+0),
			 S((x+1)-0),S(x+1),S((x+1)+0));    

    sum=0,diff=0;
    forallsites(x) {
      diff=diff+abs(S(x)-newS(x));
      S(x)=newS(x); // copy new board into board
      sum=sum+S(x); // compute number of sites alive
    }
    mdp.add(sum);   // add the local sums

    mdp << sum << endl; // print total number of sites alive

    if(diff==0) break; // if board unchanged stop
  }
  

  mdp.close_wormholes(); // STOP
  return 0;
}
