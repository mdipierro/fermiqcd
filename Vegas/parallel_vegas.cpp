/*

  MyVegas.C  (requires mdp.h)

  parallel impementation of the Vegas algorithm 

  invented by Peter Lepage @ Cornell University

  from: vegas_int.c by Sinisa Veseli @ Fermilab, 12/1996
	
  rewritten in C++ by Massimo Di Pierro @ Fermilab, 4/2000

  A number of errors corrected, in particular loops now 
  run from 0 to n-1, while in the original version they
  run from 1 to n, because transalted badly from fortran.
  Moreover it does not need numerical recepies any more but relies
  on the Marsagla random number generator in mdp.h
  The output is much nicer now.
  Moreover all the variables and functions required are
  local or private members of VegasClass.

  Many of the original varibale names have been left 
  unchanged form the original version.
   
  Goals:
  Performs Monte Carlo integration of a user supplied
  ndim-dimensional function function over a rectangular volume
  specified by { IntegrationLimitsMin[i], IntegrationLimitsMax[i] }
  The integration consists of niterations iterations, each with 
  approximately ncalls calls to the function.
  After each iteration the grid is refined; more than 5 or 10
  iterations are rarely useful. The input flag init signals
  whether this call is a new start, or a subsequent call for
  additional iterations. 
  The member variable OutputFile (by default = stdout) is the 
  (FILE*) were the basic output is directed.
  The member variable OutputFileAdvanced (normally = stdout) is 
  the (FILE*) were the advanced output is directed.

  The main function for integration is

  double VegasClass::Integrate(double (*function)(double*, void*),
                               int init, int ncalls, int niter)

  where function is the function to integrate
        init is 0,1 or 2
        ncalls if the number of calls to function per iteration
	niter is the total number of iterations

  The output of the program is in the member variables:

  double Integral;
  double StandardDeviation;
  double ChiSquare;
  
  Integral is also return by Vegas.Integrate()

*/

#include "mdp.h"

class ClassVegasPrototype {
public:
  virtual double function(double*)=0;
private:
  // constants
  static const double ALPHA=1.5;
  static const int    NDMX=200;      
  static const int    MXDIM=10;     // max grid size
  static const double TINY=1.0e-30;

  // ////////////////////////////////////////////
  // Make all variables members of the class to
  // allow restart
  // ////////////////////////////////////////////

  // variables for ClassVegas:Integrate
  int it, it1, it2, mds, nd, ndo, ng, npg;
  int ia[MXDIM], kg[MXDIM];
  double calls, dv2g, dxg, ti, tsi, xjac, xnd;
  double schi, si, swgt;
  double d[NDMX][MXDIM],  di[NDMX][MXDIM];
  double dt[MXDIM],       dx[MXDIM];
  double r[NDMX],         x[MXDIM];
  double xi[MXDIM][NDMX], xin[NDMX];

  void PrintInputParameters(int init, int ncalls, int niterations) {
    int j;
    if(OutputFile!=0 && ME==0) {
      fprintf(OutputFile,
	      "======================================\n");
      fprintf(OutputFile,
	      "Vegas Montecarlo Numerical Integration\n");
      fprintf(OutputFile,
	      "======================================\n");
      fprintf(OutputFile,
	      "Input Parameters:\n");
      fprintf(OutputFile,
	      "Number of dimensions     = %i\n", NumberOfDimensions);
      fprintf(OutputFile,"Present iteration        = %i\n", it);
      fprintf(OutputFile,"Number of iterations     = %i\n", niterations);
      fprintf(OutputFile,"Number of function calls = %i\n", ncalls );
      fprintf(OutputFile,"ALPHA                    = %i\n", ALPHA);
      fprintf(OutputFile,"mds(?)                   = %i\n", mds);
      fprintf(OutputFile,"nd                       = %i\n\n",nd);
      fprintf(OutputFile,"Integration Limits: i\tmin(x[i])\tmax(x[i])\n");
      fprintf(OutputFile,
	      "=================================================\n");
      for(j=0; j<NumberOfDimensions; j++)
	fprintf(OutputFile,"                    %i\t%f\t%f\n",
		j, IntegrationLimitsMin[j],IntegrationLimitsMax[j]);
      fprintf(OutputFile,"\n");
      fflush(OutputFile);
    }
  }
  void PrintOutputParameters() {
    int i, j;
    if(OutputFile!=0 && ME==0) {
      if(it==it1) {
	fprintf(OutputFile,
		"Iteration    Integral    StandardDeviation ChiSquare\n");
	fprintf(OutputFile,
		"====================================================\n");
      }
      fprintf(OutputFile,"%3d\t%14.7f\t%14.7f\t%10.3f\n",
	      it, Integral, StandardDeviation,ChiSquare);
      if(it==it2-1) {
	fprintf(OutputFile,
		"====================================================\n");
      }
      fflush(OutputFile);
      PrintGrid();
    }
  }
  inline int maximum(int a, int b) {
    return (a>b)?a:b;
  }
  inline int minimum(int a, int b) {
    return (a<b)?a:b;
  }

  // /////////////////////////////////////////////////////
  // Utility routine used by integrate(), to rebin
  // a vector of densities xi into new bins defined by a
  // vector r.
  // /////////////////////////////////////////////////////
  void rebin(double rc, int nd, double r[], double xin[], double xi[]) {
    int i;
    int k = 0;
    double xo;
    double dr = 0.0;
    double xn = 0.0;
    for(i=0; i<nd-1; i++) {
      while(rc > dr) {
	dr += r[k];
	xo = xn;
	xn = xi[k];
	k++;
      };
      dr -= rc;
      xin[i] = xn - (xn-xo)*dr/r[k-1];
    };
    for(i=0; i<nd-1; i++) xi[i] = xin[i];
    xi[nd-1] = 1.0;
  }

  int is_local(int k, int nk) {
    return (k % Nproc)==ME;
  }

public:
  double IntegrationLimitsMin[MXDIM];
  double IntegrationLimitsMax[MXDIM];
  int NumberOfDimensions;
  double Integral;
  double StandardDeviation;
  enum {RelativePrecision, AbsolutePrecision} ConvergenceCriteria;
  double TargetPrecision;
  double ChiSquare;
  FILE*  OutputFile;
  FILE*  OutputFileAdvanced;
  ClassVegasPrototype() {
    OutputFile=stdout;
    OutputFileAdvanced=0;
    ConvergenceCriteria=RelativePrecision;
    TargetPrecision=1e-4;
  }
  void SetConvergenceCriteriaToRelativePrecision(double x=1e-4) {
    ConvergenceCriteria=RelativePrecision;
    TargetPrecision=x;
  }
  void SetConvergenceCriteriaToAbsolutePrecision(double x=1e-4) {
    ConvergenceCriteria=AbsolutePrecision;
    TargetPrecision=x;
  }
  void SaveGrid(char filename[]) {
    if(ME==0) {
      FILE *fp=fopen(filename, "w");
      fwrite(this, sizeof(ClassVegasPrototype)/sizeof(char),1, fp);
      fclose(fp);
    }
  }
  void LoadGrid(char filename[]) {
    FILE *fp=fopen(filename, "r");
    fread(this, sizeof(ClassVegasPrototype)/sizeof(char),1, fp);
    fclose(fp);
  }
  void PrintGrid() {
    int i,j;
    if(OutputFileAdvanced!=0 && ME==0) {
      for(j=0; j<NumberOfDimensions; j++) {
	fprintf(OutputFileAdvanced,"DATA FOR axis %2d\n", j);
	fprintf(OutputFileAdvanced,"Dim\ti\t x\t\tDelta_i\n");
	fprintf(OutputFileAdvanced,
		"=============================================\n");	    
	for(i=0; i<nd; i++)
	  fprintf(OutputFileAdvanced,"%d\t%d\t%8.5f\t%12.4f\n",
		  j,i,xi[j][i], di[i][j]);
      }
    }
    fflush(OutputFile);
  }
  void PrintOutput() {
    if(OutputFile!=0 && ME==0) {
      fprintf(OutputFile, "Integral           = %f\n", Integral);
      fprintf(OutputFile, "Standard Deviation = %f\n", StandardDeviation);
      fprintf(OutputFile, "Chi^2              = %f\n", ChiSquare);
      fflush(OutputFile);
    }
  }

  void parallel_loop(double &fb, double &f2b) {
    
    static int i,j,k;
    static double wgt, xn, xo, rc, f, f2;
    static double local_d[NDMX][MXDIM],  local_di[NDMX][MXDIM];

    // /////////////////
    // initialization...
    // /////////////////
    
    fb=f2b=0.0;
    for(i=0; i<nd; i++) 
      for(j=0; j<NumberOfDimensions; j++)
	local_d[i][j]=local_di[i][j] = 0.0;
    
    // /////////////////
    // parallel loop
    // /////////////////
    for(k=0; k<npg; k++) if(is_local(k,npg)) {
      wgt=xjac;
      for(j=0; j<NumberOfDimensions; j++) {
	xn = (kg[j]-mdp_random.plain())*dxg;
	ia[j] = maximum(minimum((int)(xn), NDMX-1), 0); // check NDMX
	if(ia[j]>0) {
	  xo = xi[j][ia[j]] - xi[j][ia[j]-1];
	  rc = xi[j][ia[j]-1] + (xn-ia[j])*xo;
	} else {
	  xo = xi[j][ia[j]];
	  rc = (xn-ia[j]) * xo;
	};
	x[j] = IntegrationLimitsMin[j] + rc*dx[j];
	wgt *= xo*xnd;
      }
      f=wgt*function(x);
      f2=f * f;
      fb+=f;
      f2b+=f2;
      for(j=0; j<NumberOfDimensions; j++) {
	local_di[ia[j]][j] += f;
	if(mds >= 0) 
	  local_d[ia[j]][j] += f2;
      }
    }
    // ////////////////////////
    // Parallel sum here!
    // ////////////////////////
    mpi.add(fb);
    mpi.add(f2b);
    mpi.add((double*) local_d, nd*MXDIM);
    mpi.add((double*) local_di, nd*MXDIM);
    for(j=0; j<NumberOfDimensions; j++)
      for(i=0; i<nd; i++) {
	d[i][j]+=local_d[i][j];
	di[i][j]+=local_di[i][j];
      }
  }
  
  double Integrate(int init=0, 
		   int ncalls=1000, 
		   int niterations=100) { 
    int i, j, k;
    double wgt, f2b, fb, xo, xn, rc;
  
    if(init<=0) {
      // ////////////////////////////////////////////
      // Normal entry. Enter here on a cold start.
      // To disable stratified sampling, i.e. to use
      // importance sampling only, change to mds = 0.
      // ////////////////////////////////////////////
      mds=ndo=1;
      for(j=0; j<NumberOfDimensions; j++) xi[j][0] = 1.0;
    }
    if(init<=1) {	
      // ////////////////////////////////////////////
      // Enter here to inherit the grid from a previous
      // call, but not its answers.
      // ////////////////////////////////////////////
      it = 0;
      si = swgt = schi = 0.0;
    }
    if(init<=2) {
      // ////////////////////////////////////////////
      // Enter here to inherit the previous grid
      // and its answers.
      // ////////////////////////////////////////////
      nd = NDMX;
      ng = 1;
      if(mds) {	// Set up for stratification.
	ng = (int)pow(ncalls/2.0+0.25, 1.0/NumberOfDimensions);
	mds = 1;
	if((2*ng-NDMX)>=0) {
	  mds = -1;
	  npg = ng/NDMX + 1;
	  nd = ng / npg;
	  ng = npg * nd;
	}
      }
      for(k=1, i=0; i<NumberOfDimensions; i++) k *= ng;
      npg = maximum(ncalls/k, 2);
      calls = npg * k;
      dxg = 1.0 / ng;
      for(dv2g=1, i=0; i<NumberOfDimensions; i++) dv2g *= dxg;
      dv2g = pow(calls*dv2g/npg,2.0)/(npg-1.0);
      xnd = nd;
      dxg *= xnd;
      xjac = 1.0 / calls;
      for(j=0; j<NumberOfDimensions; j++) {
	dx[j] = IntegrationLimitsMax[j]-IntegrationLimitsMin[j];
	xjac *= dx[j];
      }
      // ////////////////////////////////////////////
      // Do binning if necessary.
      // ////////////////////////////////////////////
      if(nd!=ndo) {
	for(i=0; i<nd; i++) r[i]=1.0;
	for(j=0; j<NumberOfDimensions; j++) rebin(ndo/xnd, nd, r, xin, xi[j]);
	ndo = nd;
      }
      PrintInputParameters(init, ncalls, niterations);
    }
    
    // ////////////////////////////////////////////
    // Main iteration loop. Can enter here (init >= 3)
    // to do an additional niterations iterations with
    // all other parameters unchanged.
    // ////////////////////////////////////////////
    it1=it;
    it2=it+niterations;
    for(; it<it2; it++) {
      ti=tsi=0.0;
      for(j=0; j<NumberOfDimensions; j++) {
	kg[j]=1;
	for(i=0; i<nd; i++) 
	  d[i][j]=di[i][j] = 0.0;
      }
      int kk=0;
      do {	
        // mdp << kk << "looping\n";
        kk++;
	parallel_loop(fb, f2b);
	f2b = sqrt(f2b*npg);
	f2b = (f2b-fb) * (f2b+fb);
	if(f2b<=0.0) f2b = TINY;
	ti += fb;
	tsi += f2b;

	if(mds<0) {	// Use stratified sampling. 
	  for(j=0; j<NumberOfDimensions; j++) 
	    d[ia[j]][j] += f2b;
	}
	for(k=NumberOfDimensions-1; k>=0; k--) {
	  kg[k] %= ng;
	  if(++kg[k] != 1) break;
	}
	if(k < 0) break;
      } while(TRUE);
      // ////////////////////////////////////////////
      // Compute final results for this iteration.
      // ////////////////////////////////////////////
      tsi *= dv2g;
      wgt = 1.0 / tsi;

      si += wgt*ti;
      schi += wgt*ti*ti;
      swgt += wgt;
      Integral = si / swgt;
      if(it == 0)
	ChiSquare = 0.0;
      else 
	ChiSquare = (schi-si*(Integral))/it;
      if(ChiSquare < 0.0)
	ChiSquare = 0.0;
      StandardDeviation = sqrt(1.0/swgt);
      tsi = sqrt(tsi);
      PrintOutputParameters();
      
      // ////////////////////////////////////////////
      // Refine the grid. The refinement is damped,
      // to avoid rapid, destabilizing changes, and also
      // compressed in range by the exponent ALPHA.
      // ////////////////////////////////////////////
      for(j=0; j<NumberOfDimensions; j++) { 
	xo = d[0][j];
	xn = d[1][j];
	d[0][j] = (xo+xn)/2.0;
	dt[j] = d[0][j];
	for(i=1; i<nd-1; i++) {
	  rc = xo + xn;
	  xo = xn;
	  xn = d[i+1][j];
	  d[i][j] = (rc+xn)/3.0;
	  dt[j] += d[i][j];
	};
	d[nd-1][j] = (xo+xn)/2.0;
	dt[j] += d[nd-1][j];
      }
      for(j=0; j<NumberOfDimensions; j++) {
	rc = 0.0;
	for(i=0; i<nd; i++) {
	  if(d[i][j] < TINY) d[i][j] = TINY;
	  r[i] = pow((1.0-d[i][j]/dt[j])/(log(dt[j])-log(d[i][j])), ALPHA);
	  rc += r[i];
	}
	rebin(rc/xnd, nd, r, xin, xi[j]);
      }
      if(((ConvergenceCriteria==RelativePrecision) && 
	  (fabs(StandardDeviation/Integral)<TargetPrecision)) ||
	 ((ConvergenceCriteria==AbsolutePrecision) && 
	  (fabs(StandardDeviation)<TargetPrecision)))
	return Integral;
    } 
    if(((ConvergenceCriteria==RelativePrecision) && 
	(fabs(StandardDeviation/Integral)>TargetPrecision)) ||
       ((ConvergenceCriteria==AbsolutePrecision) && 
	(fabs(StandardDeviation)>TargetPrecision)) && ME==0)
      fprintf(OutputFile,"Vegas failed to reach target precision.\n");
    return Integral;
  }
};


class MyFunction : public ClassVegasPrototype {
public:  
  double function(double* x) {
    return sin(x[0])*cos(x[1]);
  }
} myfunction;

int main(int argc, char** argv) {
  mdp.open_wormholes(argc, argv);
  myfunction.NumberOfDimensions=2;
  myfunction.IntegrationLimitsMin[0]=0;
  myfunction.IntegrationLimitsMin[1]=0;
  myfunction.IntegrationLimitsMax[0]=1;
  myfunction.IntegrationLimitsMax[1]=1;
  mdp << "Intergal=" << myfunction.Integrate();
  mdp.close_wormholes();
  return 0;
}
