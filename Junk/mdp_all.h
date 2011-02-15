// BEGIN FILE: mdp_all.h
// C headers
#include "sys/types.h"
#include "sys/socket.h"
#include "sys/time.h"
#include "time.h"
#include "netinet/in.h"
#include "arpa/inet.h"
#include "errno.h"
#include "fcntl.h"
#include "netdb.h"
#include "signal.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "sys/stat.h"
#include "sys/uio.h"
#include "unistd.h"
#include "sys/wait.h"
#include "sys/un.h"
#include "sys/select.h"
#include "poll.h"
#include "strings.h"
#include "pthread.h"
 
// C++ headers and STL headers
#include "iostream"
#include "string"
#include "vector"
#include "deque"
#include "map"
using namespace std;
 
#ifndef HAVE_INET_NTOP
#define inet_ntop(a,b) inet_ntoa(b)
#define inet_pton(a,b,c) inet_aton(b,c)
#endif
 
void exit_message(int en, string message) {
  cerr << "FROM PROCESS PID: " << getpid() << endl;
  cerr << "CHILD OF PROCESS PID: " << getppid() << endl;
  cerr << "FATAL ERROR: " << message << endl;
  cerr << "EXITING WITH ERROR NUMBER: " << en << endl;
  exit(en);
}
