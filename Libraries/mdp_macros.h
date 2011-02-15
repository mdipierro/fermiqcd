/////////////////////////////////////////////////////////////////
/// @file mdp_macros.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains class mdp_macros
///
/// Licensed under GPL2 license
/// Read attached license in file mdp_license.pdf
/// This file cannot be distributed without file mdp_license.pdf
//////////////////////////////////////////////////////////////////

#define CHECK_ALL
#define MDP_MPI
#define INCLUDE_DEPRECATED_IO 

/// Loop on all local siltes of this process
#define forallsites(x)                                     \
        for(x.start(); x.is_in(); x.next())

/// Loop on all local sites of this process with given parity
/// If pofx is EVENODD=2 then loops on even and odd sites
#define forallsitesofparity(x,pofx)                        \
  for(x.start(), x.idx=x.lattice().start[ME][pofx % 2];               \
      x.idx<x.lattice().stop[ME][(pofx+(pofx % 2))/2];     \
      x.idx++)

/// Loop on all sites stored by this process
#define forallsitesandcopies(x)                            \
     for(x.start(), x.idx=0; x.idx<x.lattice().nvol; x.idx++)

/// Loop on all sites stored by this process with given parity

// if pofx is EVENODD=2 then loops on even and odd sites
#define forallsitesandcopiesofparity(x,pofx)                     \
     for(int __process=0; __process<Nproc; __process++)          \
     for(x.start(), x.idx=x.lattice().start[__process][pofx % 2];           \
         x.idx<x.lattice().stop[__process][(pofx+(pofx % 2))/2]; \
         x.idx++)

/// Returns the unique id of this process
#define ME     mpi.me()
/// Returns the total number of parallel processes for this job
#define Nproc  mpi.nproc()

/// Reports a runtime error and the line that caused it
#define error(a) _mpi_error_message(a,__FILE__, __LINE__);

#ifndef TRUE
#define TRUE true
#endif
#ifndef FALSE
#define FALSE false
#endif


