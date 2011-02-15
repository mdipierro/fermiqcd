/////////////////////////////////////////////////////////////////
/// @file mdp_postscript.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Yes...MDP can print and draw in postscript
///
/// Licensed under GPL2 license
/// Read attached license in file mdp_license.pdf
/// This file cannot be distributed without file mdp_license.pdf
//////////////////////////////////////////////////////////////////

/// @brief to output and draw in postscript
/// 
/// Example:
/// @verbatim
///    mdp_postscript ps("test.ps");
///    ps.color(0.2,0.2,0.7);
///    ps.line(0,0,  5,5);
///    ps.print(5,5,"a line from (0,0) to here");
/// @endverbatim
class mdp_postscript {
public:
  enum {BOLD=10};
  FILE *fp;
  float X0, Y0, Z0;
  float X1, Y1, Z1;
  float c0, c1, c2;
  float scale;
  float alpha;
  mdp_postscript() {
    fp=0;
    scale=1;
  };
  mdp_postscript(char filename[]) {
    open(filename);
  };
  virtual ~mdp_postscript() {
    if (fp!=0) close();
  };
  FILE* open(char filename[]) {
    printf("Making frame: %s\n", filename);
    fp=fopen(filename, "w");
    if(fp!=0) fprintf(fp,"%%!PS-Adobe-3.0 EPSF-3.0\n");
    fflush(fp);
    return fp;
  };
  void close() {
    fprintf(fp,"showpage\n");
    fprintf(fp,"%%%%Trailer\n");
    fprintf(fp,"%%%%EOF\n");
    fflush(fp);
    fclose(fp);
    fp=0;
  };
  void size(float x0, float y0, float x1, float y1) {
    X0=x0; Y0=y0;
    X1=x1; Y1=y1;
    fprintf(fp,"%%%%BoundingBox: %.0f %.0f %.0f %.0f\n", scale*x0, scale*y0, scale*x1, scale*y1);
    fprintf(fp,"%%%%EndComments\n");
    fflush(fp);
  };
  void line(float x0,  float y0, float x1, float y1) {
    fprintf(fp,"%.2f %.2f %.2f setrgbcolor\n", c0,c1,c2);
    fprintf(fp,"%.2f %.2f moveto\n", scale*x0, scale*y0);
    fprintf(fp,"%.2f %.2f lineto\n", scale*x1, scale*y1);
    fprintf(fp,"stroke\n");
    fflush(fp);
  };
  void box(float x0,  float y0, float x1, float y1, int fill=0) {
    if (fill==1) fprintf(fp,"stroke\n");
    fprintf(fp,"%.2f %.2f moveto\n", scale*x0, scale*y0);
    fprintf(fp,"%.2f %.2f lineto\n", scale*x0, scale*y1);
    fprintf(fp,"%.2f %.2f lineto\n", scale*x1, scale*y1);
    fprintf(fp,"%.2f %.2f lineto\n", scale*x1, scale*y0);
    fprintf(fp,"%.2f %.2f lineto\n", scale*x0, scale*y0);
    if (fill==1) fprintf(fp,"fill\n");
    fprintf(fp,"stroke\n");
    fflush(fp);
  };
  void arc(float x0,  float y0, float r, float alpha, float beta) {
    fprintf(fp,"%.2f %.2f %.2f %.2f %.2f arc\n", scale*x0, scale*y0, scale*r, alpha, beta);
    fprintf(fp,"stroke\n");
    fflush(fp);
  };
  void circle(float x0,  float y0, float r, int fill=0) {
    if (fill!=BOLD) fprintf(fp,"%.2f %.2f %.2f 0 360 arc\n", scale*x0, scale*y0, r);
    if (fill==1) fprintf(fp,"fill\n");
    if (fill==BOLD) {
      fprintf(fp,"1 -0.1 0.1 {\n");
      fprintf(fp,"/r exch def\n");
      fprintf(fp,"%.2f 1 r 0.5 mul sub mul ", c0);
      fprintf(fp,"%.2f 1 r 0.5 mul sub mul ", c1);
      fprintf(fp,"%.2f 1 r 0.5 mul sub mul ", c2);
      fprintf(fp,"setrgbcolor\n");
      fprintf(fp,"stroke\n");
      fprintf(fp,"%.2f %.2f %.2f r mul 0 360 arc fill\n", scale*x0, scale*y0, scale*r);
      fprintf(fp,"} for\n");
    };
    fflush(fp);
  };
  void pen(float size) {
     fprintf(fp,"%.2f setlinewidth\n", size);
  };
  // colors are numbers in [0,1] black=(0,0,0) white=(1,1,1)
  void color(float r, float g, float b) {
     fprintf(fp,"%.2f %.2f %.2f setrgbcolor\n", r,g,b);
     c0=r; c1=g; c2=b;
     fflush(fp);
  };
  void font(const char *text, int size) {
    fprintf(fp,"/%s findfont\n%i scalefont\nsetfont\n", text, size);
    fflush(fp);
  };
  void print(float x0, float y0, char text[]) {
    fprintf(fp,"%.2f %.2f moveto\n", scale*x0, scale*y0);
    fprintf(fp,"(%s) show\n", text);
    fprintf(fp,"stroke\n");
    fflush(fp);
  };
};

