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
#include "sys/file.h"
#include "sys/un.h"
#include "poll.h"
#include "strings.h"
#include "pthread.h"

#ifdef HAWK
#define socklen_t int
#define FASYNC 020000
#else
#include "sys/select.h"
#endif
 
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

class InternetAddress {
public:  
  struct sockaddr_in address;
  string ipaddress;
  int port; 
  InternetAddress(string hostname="127.0.0.1", int port=0) {

    char tmp[16];
    struct hostent *h=gethostbyname(hostname.c_str());
    if(h==0) throw string("Invalid hostname");
    if(h->h_length!=4) throw string("Invalid hostname");
    strncpy(tmp,inet_ntop(AF_INET,*((struct in_addr*)  h->h_addr_list[0])),16);
    ipaddress=tmp;    

    this->ipaddress=ipaddress;
    this->port=port;
    memset(&address,0,sizeof(address));
    address.sin_family=AF_INET;
    address.sin_port=htons(port);
    if(!inet_pton(AF_INET,ipaddress.c_str(),&address.sin_addr))
      exit_message(1,"invalid IP address");    
  } 
  int getPort() {
    return ntohs(address.sin_port);    
  }
  string getIPAddress() {
    return inet_ntop(AF_INET,address.sin_addr);
  }
  int Connect(int sfd, int timeout=0) {
    if(! timeout) {
      return connect(sfd,(struct sockaddr*)&address,
		     (socklen_t) sizeof(address));
    } else {      
      int arg = fcntl(sfd, F_GETFL, NULL); 
      arg |= (O_NONBLOCK); 
      fcntl(sfd, F_SETFL, arg);
      pollfd fds;
      fds.fd=sfd;
      fds.events=POLLWRNORM;
      int ret=connect(sfd,(struct sockaddr*)&address,
		     (socklen_t) sizeof(address));      
      if(ret<0) {
	if(errno!=EINPROGRESS) return ECONNREFUSED;
	ret=poll(&fds,1,timeout);
	arg = fcntl(sfd, F_GETFL, NULL); 
	arg &= (~O_NONBLOCK); 
	fcntl(sfd, F_SETFL, arg);
	int valopt;
	socklen_t lon=sizeof(int);
	if (getsockopt(sfd,SOL_SOCKET,SO_ERROR,(void*)&valopt,&lon)<0 || 
	    valopt) return ECONNREFUSED;	  
	if(ret>0) return 0;
	if(ret==0) return ETIMEDOUT;
	return ECONNREFUSED;	  
      } else {
	arg = fcntl(sfd, F_GETFL, NULL); 
	arg &= (~O_NONBLOCK); 
	fcntl(sfd, F_SETFL, arg);
	return ret;
      }
    }
    return 0;
  }
  int Accept(int fd) {
    socklen_t ssize=sizeof(struct sockaddr_in);
    return accept(fd,(struct sockaddr*) &address, &ssize);
  }
  int sendTo(int sfd, void* data, int size, int options=0) {
    return sendto(sfd,data,size,options,
		  (struct sockaddr*)&address,
		  (socklen_t) sizeof(address));
  }
  int recvFrom(int sfd, void* data, int size, int options=0) {
    socklen_t address_size=sizeof(address);
    return recvfrom(sfd,data,size,options,
		    (struct sockaddr*)&address,
		    (socklen_t*) &address_size);
  }
};

int newTcpSocket(int flags=0) {
  return socket(AF_INET,SOCK_STREAM,flags);
}

int newUdpSocket(int flags=0) {
  return socket(AF_INET,SOCK_DGRAM,flags);
}

int newTcpClientSocket(string ipaddress, int port, int sleep_time=10) {
  InternetAddress peer=InternetAddress(ipaddress,port);
  int sfd=newTcpSocket();
  while(1) {
    int r=peer.Connect(sfd);
    if(r==0) return sfd;
    if(r!=ETIMEDOUT) {
      close(sfd);
      return -1;
    }
    sleep(sleep_time);
  }
}

int Bind(int sfd, int port) {
  struct sockaddr_in address;
  memset(&address,0,sizeof(address));
  address.sin_family=AF_INET;
  address.sin_port=htons(port);
  address.sin_addr.s_addr=htonl(INADDR_ANY);
  return bind(sfd,(struct sockaddr*) &address,(socklen_t) sizeof(address));  
}

int setSocketKeepAlive(int sfd, int on=1) {
  return setsockopt(sfd,SOL_SOCKET,SO_KEEPALIVE,&on,sizeof(on));
}

int setFileLock(int sfd) {
  return flock(sfd,LOCK_EX);
}

int setFileLockShared(int sfd) {
  return flock(sfd,LOCK_SH);
}

int setFileUnlock(int sfd) {
  return flock(sfd,LOCK_UN);
}

int setSocketNonblock(int sfd, int on=1) {
  int flags=fcntl(sfd,F_GETFL,0);
  if(on) flags=flags | O_NONBLOCK;
  else   flags=flags & (~O_NONBLOCK);
  return fcntl(sfd,F_SETFL,flags);
}

int setSocketAsync(int sfd, int owner=0, int on=1) {
  if(owner==0) owner=getpid();
  int flags=fcntl(sfd,F_GETFL,0);
  if(on) flags=flags | O_NONBLOCK | FASYNC;
  else   flags=flags & (~O_NONBLOCK) & (~FASYNC);
  fcntl(sfd,F_SETFL,flags);
  fcntl(sfd,F_SETOWN,owner);
  return 1;
}

int setSocketReusable(int sfd, int on=1) {  
  return setsockopt(sfd,SOL_SOCKET,SO_REUSEADDR,&on,sizeof(on));
}

int setSocketSendBroadcast(int sfd, int on=1) {  
  return setsockopt(sfd,SOL_SOCKET,SO_BROADCAST,&on,sizeof(on));

}

int setSocketRecvBroadcast(int sfd, int on=1) {
  return setSocketReusable(sfd,1);
}

int setSocketSendMulticast(int sfd, int ttl=1) {
  return setsockopt(sfd,IPPROTO_IP,IP_MULTICAST_TTL,&ttl,sizeof(ttl));
}

int setSocketRecvMulticast(int sfd, string from) {
  setSocketReusable(sfd,1);
  struct ip_mreq mreq;
  mreq.imr_multiaddr.s_addr=InternetAddress(from,0).address.sin_addr.s_addr;
  mreq.imr_interface.s_addr=htonl(INADDR_ANY);
  return setsockopt(sfd,IPPROTO_IP,IP_ADD_MEMBERSHIP,&mreq,sizeof(mreq));
}

int forkTwice() {
  pid_t pid;
  int status;
  pid=fork();
  if(!pid) {
    int pid2=fork();
    switch (pid2) {
    case 0:  return 0;
    case -1: _exit(errno);    // assumes all errnos are <256 
    default: _exit(0);
    }
  }
  if(pid>0) {
    if(waitpid(pid,&status,0)<0)
      return -1;   
    else
      return 1;
  }
  if(pid<0) {
    return -1;
  }
  return -1;
}

class Thread {
private:
  pthread_t _thread_number;
  pthread_attr_t _thread_attributes;
 public:
  Thread();
  pthread_t threadNumber();
  void threadStart();
  void threadStop();
  void threadSetJoinable();
  void threadJoin();
  virtual void run()=0;
};

void* thread_function(void *p) {
  Thread* pt=(Thread*) p;
  pt->run();
}

Thread::Thread() {
    pthread_attr_init(&_thread_attributes);
}

void Thread::threadStart() {
  pthread_create(&_thread_number,&_thread_attributes,thread_function,(void*) this);
}

void Thread::threadStop() {
  pthread_exit(NULL);
}

void Thread::threadSetJoinable() {
  pthread_attr_setdetachstate(&_thread_attributes,PTHREAD_CREATE_JOINABLE);
}

void Thread::threadJoin() {
  pthread_join(_thread_number,NULL);
}

void _handler(int);

class SignalHandler {
 private:
  struct sigaction action;
 public:
  void catch_signal(int signalnum=SIGALRM);
  virtual void handler(int s)=0;
};

map<int,SignalHandler*> _signal_handlers;
void SignalHandler::catch_signal(int signalnum) {
  if(sigemptyset(&action.sa_mask)!=0)
    throw string("sigemptyset");
  if(sigaddset(&action.sa_mask,signalnum)!=0)
    throw string("sigaddset");
  action.sa_handler=_handler;
  if(sigaction(signalnum,&action,NULL)!=0)
    throw string("sigaction");
  _signal_handlers[signalnum]=this;
}

void _handler(int s) {
  _signal_handlers[s]->handler(s);
}

//
// send the signum signal to the process pid
//
void signalSend(pid_t pid, int signum) {
  kill(pid,signum);
}

//
// block the signal signum
//
sigset_t signalBlock(int signum) {
  sigset_t set, oset;
  if(sigemptyset(&set)<0) throw string("sigemptyset");
  if(sigaddset(&set,signum)<0) throw string("sigaddset");
  sigprocmask(SIG_BLOCK,&set,&oset); 
  return oset;
}

//
// check if signal signum is pending
//
bool signalPending(int signum) {
  sigset_t set;
  sigpending(&set); 
  return sigismember(&set,signum); 
}

//
// unblock the signal signum
//
sigset_t signalUnblock(int signum) {
  sigset_t set, oset;
  if(sigemptyset(&set)<0) throw string("sigemptyset");
  if(sigaddset(&set,signum)<0) throw string("sigaddset");
  sigprocmask(SIG_UNBLOCK,&set,&oset); 
  return oset;
}

//
// restore signal blocks
//
sigset_t signalRestoreBlocks(sigset_t set) {
  sigset_t oset;  
  sigprocmask(SIG_BLOCK,&set,&oset);   
  return oset;
}
