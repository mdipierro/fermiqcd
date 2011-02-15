/////////////////////////////////////////////////////////////////
/// @file mdp_prompt.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Functions to parse user input of parameters
/// in a way safe to parallel programs
///
/// Licensed under GPL2 license
/// Read attached license in file mdp_license.pdf
/// This file cannot be distributed without file mdp_license.pdf
//////////////////////////////////////////////////////////////////

const char STD_INPUT[] = "";
const char STD_INPUT_FILE[] = "<stdin>";

/// Converts string to float
double val(string s) {
  return atof(s.c_str());
}


/// Try prompt("<stdin>","VALUE","4.0")
/// It will prompt the user for variable VALUE and take 4.0 as default
string prompt(string filename, 
	      string variable, 
	      string def_val="0.0", 
	      int p=0) {
  FILE *fp=0;
  static char tmp[1024];
  int active, length;
  char c,d;
  char response[1024];

#ifdef PARALLEL
  mpi.barrier();
#endif

  if(ME==p) {
    if(filename==STD_INPUT) {
      printf("Input value of %s (default is `%s'): ", 
	     variable.c_str(), def_val.c_str());
      cin >> response;
    } else {
      if(filename==STD_INPUT_FILE) { 
	fp=stdin;
	printf("Input variable (default is `%s %s')...\n", 
	       variable.c_str(), def_val.c_str());
      } else {
	fp=fopen(filename.c_str(),"r");
      }
      fseek(fp, 0, 0);
      active=1;
      length=0;
      d=0;
      tmp[0]='\0';
      response[0]='\0';
      do {
	c=fgetc(fp);	
	switch(active) {
	case 0: if(c=='\n') {active=1; length=0;} break;
	case 1: if(c=='#') { 
	  active=0;
	} else if((length>0) && ((c==' ') || (c=='\t'))) {
	  active=2; length=0;
	} else if((length>0) && (c=='\n')) {
	  active=1; length=0;
	  if(variable==tmp) c=EOF;
	} else if((c!=' ') && (c!='\t') && (c!='\n')) {
	  tmp[length]=c; tmp[length+1]='\0'; length++; 
	} break;	
	case 2: if(c=='#') { 
	  active=0;
	} else if((c=='\n') || (c==EOF)) {
	  active=1;
	  length=0;
	  if(variable==tmp) d=1;
	} else if((length==0) && ((c==' ') || (c=='\t') || (c=='='))) {
	  // do nothing
	} else if((c==' ') || (c=='\t')) {
	  active=0;
	  if(variable==tmp) d=1;
	} else { 
	  response[length]=c; response[length+1]='\0'; length++; 
	} break;
	}
	// printf("%i %c %s %s\n", active, c, tmp, response);
      } while ((c!=EOF) && (d!=1));

      if((d==0) || (strcmp(response,"")==0)) 
	strcpy(response, def_val.c_str());
      
      if(fp!=stdin && fp!=0) 
	fclose(fp);
    }
    printf("... Adopting %s equal to \"%s\"\n\n", variable.c_str(), response);
  }

#ifdef PARALLEL
  mpi.broadcast(response, 1024, p);
#endif
  string s(response);
  return s;
}





