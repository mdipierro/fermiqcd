/////////////////////////////////////////////////////////////////
/// @file mdp_psim.h
/// @version 2009-12-21
/// @author Massimo Di Pierro <mdipierro@cs.depaul.edu>
///
/// Contains class mdp_psim (the parallel simulator)
///
/// Licensed under GPL2 license
/// Read attached license in file mdp_license.pdf
/// This file cannot be distributed without file mdp_license.pdf
//////////////////////////////////////////////////////////////////

#include "cstdio"
#include "cstdlib"
#include "cmath"
#include "cstring"
#include "string"
#include "iostream"
#include "string"
#include "vector"
#include "map"
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/socket.h>
#include <sys/errno.h>
#include <fcntl.h>
using namespace std;

#ifndef FERMIQCD
typedef int mdp_int; // this in case this lib is used without fermiqcd
#endif

/// @brief Parallel SIMulator used by class mdp_communicator
///
/// Attention: under MDP and/or FermiQCD this is already
/// Instantiated inside class mdp_communicator.
///
/// Example:
/// @verbatim
/// int main(int argc, char** argv) {
///    mdp_psim node(argc,argv);
///    int a=3, b=0;
///    if(node.id()==0) node.send(1,a);
///    if(node.id()==1) { node.recv(0,b); cout << b << endl;
///    return 0;
/// }
/// @endverbatim
/// Compile with
/// @verbatim
///    g++ [filename] -o a.out
/// @endverbatim
/// and run with
/// @verbatim
///    ./a.out -PSIM_NPROCS=2
/// @endverbatim
/// Output should be 3.
class mdp_psim {
private:
 
  // Typed Constants
  const static int PROCESS_COUNT_MIN= 1;  // minimum number of processes
  const static int PROCESS_COUNT_MAX= 128;// maximum number of processes
  const static int CONN_LIST_IGNORE = -1; // connections that cannot occur
  const static int PROCESS_PARENT   = 0;  // The parent process ID number
  const static int COMM_RECV = 0;      // socket array indicator for reading
  const static int COMM_SEND= 1;  // socket array indicator for writing
  const static int COMM_TIMEOUT_DEFAULT = 86400;  // 1 day default
  
  // common enum values for logging routines
  enum enumBegEnd       { LOG_BEGIN, LOG_END };
  enum enumSendRecv     { LOG_SR_SEND, LOG_SR_RECV };
  enum enumSendRecvStep { LOG_SR_START, LOG_SR_FAIL, LOG_SR_SUCCESS };
  
  // Set this to true if you want verbose testing description
  
  // Class variables
  int     _verbatim;              // 0 for no output, 1 for some, 2 more
  int     _processCount;          // Holds the number of processes
  string  _logFileName;           // filename of the log file
  bool    _doLogging;             // do logging or not?
  FILE*   _logfileFD;             // file descriptor for the logging file
  int     _processID;             // process ID of "this" process
  
  int   (*_socketFD)[2];          // array to hold all of the sockets
  int _commTimeout;               // defaults to COMM_TIMEOUT_DEFAULT
  
  map< string, vector<char> >* _hash;  
  // Hash Map to hold out of sequence (send/receive) data

  // ******************************************************************* 
  // ***         Private Method: psim_begin                          *** 
  // ***                                                             *** 
  // ***  Used by the constructor ONLY                               *** 
  // ******************************************************************* 

  void psim_begin(int processCount, string logFileName, int verbatim) {
    _processCount=processCount;
    _logFileName=logFileName;
    _verbatim=verbatim;

    open_log();
    if ((processCount<PROCESS_COUNT_MIN) || (processCount<PROCESS_COUNT_MIN)) {
      log("PSIM ERROR: Invalid number of processes");
      throw string("PSIM ERROR: Invalid number of processes");
    }
    
    initialize(processCount);
    create_sockets();
    fork_processes();
    close_sockets();

    char buffer[256];
    sprintf(buffer, "process %i of %i created with pid=%i",
	    _processID, _processCount, getpid());
    log(buffer,1);
    
  }

  // ******************************************************************* 
  // ***         Private Method: psim_end                            *** 
  // ***                                                             *** 
  // ***  Used by the destructor ONLY                                ***
  // ***                                                             *** 
  // ******************************************************************* 

  void psim_end() {
    for (int source=0; source<_processCount; source++) {
      for (int dest=0; dest<_processCount; dest++) {
	if(dest==_processID) close(_socketFD[_processCount*source+dest][COMM_SEND]);
	if(source==_processID) close(_socketFD[_processCount*source+dest][COMM_RECV]);
      }
    }      
    if (_socketFD != NULL)
      delete [] _socketFD;
    _socketFD = NULL;
    
    // only delete the _hash if it still exists
    if (_hash != NULL) {
      delete [] _hash;
      _hash = NULL;
    }
    
    char buffer[256];
    sprintf(buffer, "process %i terminating", _processID);
    log(buffer,1);

    close_log();
  }

  // ******************************************************************* 
  // ***         Private Method: initialize                          *** 
  // ***                                                             *** 
  // ***  Used by the constructor, this method sets up values and    *** 
  // ***  some of the needed resources.                              *** 
  // ***                                                             *** 
  // ******************************************************************* 

  void initialize(int processCount) {
    _processCount = processCount;    
    _processID = PROCESS_PARENT;    
    _commTimeout = COMM_TIMEOUT_DEFAULT;
    
    _hash = new map< string, vector<char> >[_processCount];    
    if (_hash == NULL) {
      log("PSIM ERROR: failure to allocate hash");
      throw string("PSIM ERROR: failure to allocate hash");
    }
    _socketFD = new int[_processCount*_processCount][2];    
    if (_socketFD == NULL) {
      log("PSIM ERROR: failed to create socket array");    
      throw string("PSIM ERROR: failed to create socket array");    
    }
  }
  
  // ******************************************************************* 
  // ***         Private Method: create_sockets                      *** 
  // ***                                                             *** 
  // ***  Opens all of the sockets necessary for communication and   *** 
  // ***  inserts them into an array.                                *** 
  // ***                                                             *** 
  // ******************************************************************* 

  void create_sockets() { 
    char filename[512];
    for (int source=0; source<_processCount; source++) {
      for (int dest=0; dest<_processCount; dest++) {
	
	sprintf(filename,".fifo.%i.%i", source, dest);
	while(mkfifo(filename,0666)<0) {
	  if(errno==EEXIST) unlink(filename);
	  else throw string("PSIM ERROR: unable to mkfifo ")+filename;
	}
	_socketFD[_processCount*source+dest][COMM_RECV] = open(filename,O_RDONLY|O_NONBLOCK);
	_socketFD[_processCount*source+dest][COMM_SEND] = open(filename,O_WRONLY);
	int flags = fcntl(_socketFD[_processCount*source+dest][COMM_RECV], F_GETFL, 0);
	fcntl(_socketFD[_processCount*source+dest][COMM_RECV], F_SETFL, flags & ~O_NONBLOCK);
	if (_socketFD[_processCount*source+dest][COMM_SEND]<=0 ||
	    _socketFD[_processCount*source+dest][COMM_RECV]<=0) {
	  throw string("PSIM ERROR: unable to open fifo");
	}
      }
    }      
    char buffer[256];    
    for (int source=0; source<_processCount; source++) 
      for (int dest=0; dest<_processCount; dest++) {
	sprintf(buffer,"_socketFD[%i*%i+%i]={%i,%i}",
		source,_processCount,dest,
		_socketFD[_processCount*source+dest][COMM_SEND],
		_socketFD[_processCount*source+dest][COMM_RECV]);
	log(buffer);
      }
  }
  void close_sockets() { 
    for (int source=0; source<_processCount; source++) 
      for (int dest=0; dest<_processCount; dest++) {
	if(dest!=_processID) close(_socketFD[_processCount*source+dest][COMM_RECV]);
	if(source!=_processID) close(_socketFD[_processCount*source+dest][COMM_SEND]);
      } 
  }     
  
  // ******************************************************************* 
  // ***         Private Method: fork_processes                      *** 
  // ***                                                             *** 
  // ***                                                             *** 
  // ******************************************************************* 
  void fork_processes() {
    _processID=0;
    for (int i=1; i<_processCount; i++) {
      int pid = fork();
      
      if (pid == -1) {
	log("PSIM ERROR: fork");
	throw("PSIM ERROR: fork");
      } else if (pid == 0) { // Child Process
	_processID = i;
	break;
      }
    }
  }    

  // ******************************************************************* 
  // ***         Private Method: check_process_id                    *** 
  // ***                                                             *** 
  // ***  Varifies that the destination process ID is valid. This    *** 
  // ***  is done before data is sent or received.                   *** 
  // ***                                                             *** 
  // ******************************************************************* 
    
  void check_process_id(int processID) {
    
    if ((processID == _processID) || 
	(processID < PROCESS_PARENT) ||
	(processID >= _processCount)) {
      
      char buffer[256];
      sprintf(buffer, "PSIM ERROR: process %i referencing %i.",
	      _processID, processID);
      log(buffer);
      throw string( buffer );
    }
  }
  
  // ******************************************************************* 
  // ***         Private Method: open_log                            *** 
  // ***                                                             *** 
  // ***  This method initializes the process log and sets it up for *** 
  // ***  appending messages.                                        *** 
  // ***                                                             *** 
  // ******************************************************************* 
  
  void open_log() {
    _doLogging = false;
    if (_logFileName.length()==0) {      
      return;
    }
    
    // open and reset file
    if ((_logfileFD = fopen(_logFileName.c_str(), "w")) == NULL) {
      log("PSIM ERROR: unable to create logfile");
      throw string("PSIM ERROR: unable to create logfile");
    }
    // close the log file
    close_log();
    
    // reopen the log file in append mode
    if ((_logfileFD = fopen(_logFileName.c_str(), "a")) == NULL) {
      log("PSIM ERROR: unable to open logfile");
      throw string("PSIM ERROR: unable to open logfile");
    }
    _doLogging=true;
    
  }
  
  
  // ******************************************************************* 
  // ***         Private Method: close_log                           *** 
  // ***                                                             *** 
  // ***  Closes the log file.                                       *** 
  // ***                                                             *** 
  // ******************************************************************* 
  
  void close_log() {
    if (_doLogging)
      fclose(_logfileFD);
    
  }
  
  
  // ******************************************************************* 
  // ***         Private Method: logSendRecv                         *** 
  // ***                                                             *** 
  // ***  Centralizes the repetitive task of logging the steps       *** 
  // ***  during send and receive.                                   *** 
  // ***                                                             *** 
  // ******************************************************************* 
  
  void logSendRecv(int sourcedestProcessID,
		   string tag,
		   enumSendRecv method, 
		   enumSendRecvStep step){
    
    char buffer[256];
    const char cmdSendRecv[2][8] = { "send", "recv" };
    const char stepSendRecv[3][12] = { "starting...", "failed!", "success." };
    
    sprintf(buffer, "%i %s(%i,%s) %s", 
	    _processID, cmdSendRecv[method], 
	    sourcedestProcessID, tag.c_str(), stepSendRecv[step]);
    log(buffer);
  }
  
  // ******************************************************************* 
  // ***         Private Method: get_source_index                    *** 
  // ***                                                             *** 
  // ******************************************************************* 
  
  int get_source_index(int source) {
    check_process_id(source);
    return _processCount*source+_processID;
  }

  // ******************************************************************* 
  // ***         Private Method: detDestIndex                        *** 
  // ***                                                             *** 
  // ******************************************************************* 
  
  int get_dest_index(int dest) {
    check_process_id(dest);
    return _processCount*_processID+dest;
  }
  
  // ******************************************************************* 
  // ***         Private Method: send_buffer                         *** 
  // ***                                                             *** 
  // ***  Handles the sending of binary data.                        *** 
  // ***                                                             *** 
  // ******************************************************************* 

  void send_buffer(int destProcessID, 
		   const void* pdataToSend, mdp_int dataSize) {
    static int counter=0;
    char filename[512];
    counter++;
    int destIndex = get_dest_index(destProcessID);
    sprintf(filename,".fifo.%i.%i.%i",_processID,destProcessID,counter);
    int fd = open(filename,O_WRONLY|O_CREAT,0700);    
    if (write(fd, pdataToSend, dataSize) != dataSize) {
      log("PSIM ERROR: failure to write to socket");
      throw string("PSIM ERROR: failure to write to socket");
    }
    if(write(_socketFD[destIndex][COMM_SEND], &counter, sizeof(counter))!= sizeof(counter)) {
      log("PSIM ERROR: failure to write to socket");
      throw string("PSIM ERROR: failure to write to socket");
    }
    close(fd);
  }
  
  
  // ******************************************************************* 
  // ***         Private Method: send_binary                         *** 
  // ***                                                             *** 
  // ***  Sends a data tag and a vector of chars (as binary data).   *** 
  // ***                                                             *** 
  // ******************************************************************* 
  
  void send_binary(int destProcessID, 
		   const string& tag, 
		   const vector<char>& data) {
    
    int tagSize = tag.size();
    int dataSize = data.size();
    send_buffer(destProcessID, &tagSize, sizeof(tagSize));    
    send_buffer(destProcessID, tag.c_str(), tagSize);    
    send_buffer(destProcessID, &dataSize, sizeof(dataSize));
    send_buffer(destProcessID, &data[0], dataSize);
  }
  
  
  // ******************************************************************* 
  // ***         Private Method: recv_buffer                         *** 
  // ***                                                             *** 
  // ***  Handles the receiving of binary data through the sockets.  *** 
  // ***                                                             *** 
  // ******************************************************************* 
  
  void recv_buffer(int sourceProcessID, 
		   void* pdataToReceive, mdp_int dataSize) {
    char filename[512];
    int counter;
    int sourceIndex = get_source_index(sourceProcessID);
    if(read(_socketFD[sourceIndex][COMM_RECV],&counter,sizeof(counter))!=sizeof(counter)) {
      log("PSIM ERROR: timeout error in readin from socket");
      throw string("PSIM ERROR: timeout error in readin from socket");
    }
    sprintf(filename,".fifo.%i.%i.%i",sourceProcessID,_processID,counter);
    int fd = open(filename,O_RDONLY);
    if(read(fd,(char*) pdataToReceive,dataSize)!=dataSize) {
      log("PSIM ERROR: timeout error in readin from socket");
      throw string("PSIM ERROR: timeout error in readin from socket");
    } else {
      unlink(filename);
    }
    close(fd);
  }
  
  
  // ******************************************************************* 
  // ***         Private Method: recv_binary                         *** 
  // ***                                                             *** 
  // ***  Receives data utilizing a data tag to make sure that the   *** 
  // ***  data coming in is what was expected.                       *** 
  // ***                                                             *** 
  // ******************************************************************* 
  
  void recv_binary(int sourceProcessID, 
		   const string& tag, 
		   vector<char>& data ) {
    
    static vector<char> dataTmp;
    static string tagReceived;
    int size;
    
    map< string, vector<char> >::iterator itr;
    
    while (true) {
      itr = _hash[sourceProcessID].find(tag);
      
      if (itr != _hash[sourceProcessID].end()) { // Found?
	data = _hash[sourceProcessID][tag];
	_hash[sourceProcessID].erase(itr);
	break;
      } else {
	recv_buffer(sourceProcessID, &size, sizeof(size));
	char* buffer= new char[size+1];
	recv_buffer(sourceProcessID, buffer, size);
	buffer[size] = 0;
	tagReceived = buffer;
	delete buffer;

	if (tagReceived == tag) {
	  recv_buffer(sourceProcessID, &size, sizeof(size));
	  data.resize(size);	  
	  recv_buffer(sourceProcessID, &data[0], size);
	  break;
	} else {
	  recv_buffer(sourceProcessID, &size, sizeof(size));
	  dataTmp.resize(size);
	  recv_buffer(sourceProcessID, &dataTmp[0], size);
	  _hash[sourceProcessID][tagReceived] = dataTmp;
	}
      }
    }
  }
  
  
  // ******************************************************************* 
  // ***         Private Method: doBegEndLog                         *** 
  // ***                                                             *** 
  // ***  Centralized log method used by some of the public methods  *** 
  // ***  to send a standardized message to the log at the beginning *** 
  // ***  and the end of the routine.                                *** 
  // ***                                                             *** 
  // ******************************************************************* 
  
  void doBegEndLog(string method, enumBegEnd begEnd) {
    char buffer[256];
    char* be;
    
    if (begEnd == LOG_BEGIN) be = (char*)"BEGIN";
    else  be = (char*)"END";
    
    sprintf(buffer, "%i %s [%s]", _processID, be, method.c_str());
    log(buffer);
  }
    
public:
  
  // ******************************************************************* 
  // ******************************************************************* 
  // ***                                                             *** 
  // ***                P U B L I C   M E T H O D S                  *** 
  // ***                                                             *** 
  // ******************************************************************* 
  // ******************************************************************* 
  
  
  // ******************************************************************* 
  // ***               Constructor: mdp_psim                         *** 
  // ***                                                             *** 
  // ***  Provide the number of processes to create and the name of  *** 
  // ***  the logfile if desired and "" if no logfile is needed.     *** 
  // ***                                                             *** 
  // ******************************************************************* 
  
  mdp_psim(int processCount, string logFileName=".psim.log", int verbatim=0) {
    psim_begin(processCount, logFileName, verbatim);
  }
  
  mdp_psim(int argc, char** argv) {
    int processCount=parse_argv_nprocs(argc,argv);
    string logFileName=parse_argv_logfile(argc,argv);
    int verbatim=parse_argv_verbatim(argc,argv);
    psim_begin(processCount, logFileName, verbatim);
  }
  
  
  // ******************************************************************* 
  // ***               Destructor: ~mdp_psim                         *** 
  // ***                                                             *** 
  // ***  Deallocates space that was created within the process,     *** 
  // ***  releases sockets, closes the log, etc.                       *** 
  // ***                                                             *** 
  // ******************************************************************* 
  
  virtual ~mdp_psim() {
    psim_end();
  }
  
  
  // ******************************************************************* 
  // ***         Public Method: log                                  *** 
  // ***                                                             *** 
  // ***  Accepts a string and appends the message to the common     *** 
  // ***  log file.  Note: locking is not necessary because of the   *** 
  // ***  deffinition of append.  It does not matter how many        *** 
  // ***  processes share file pointers, writing will always occur   *** 
  // ***  at the end of the file.                                    *** 
  // ***                                                             *** 
  // ******************************************************************* 
  
  void log(string message,int level=2) {    
    if (_doLogging) {
      int fd=fileno(_logfileFD);
      flock(fd,LOCK_EX);      
      fwrite("PSIM LOG: ", 10, 1, _logfileFD);
      fwrite(message.c_str(), message.length(), 1, _logfileFD);
      fwrite("\n", 1, 1, _logfileFD);
      // Clear out the file buffer
      fflush(_logfileFD);
      flock(fd,LOCK_UN);
    }
    if(_verbatim>=level) {
      cout << "PSIM LOG: " << message << endl;
      cout.flush();
    }
  }
  
  
  // ******************************************************************* 
  // ***         Public Method: id                                   *** 
  // ***                                                             *** 
  // ***  Returns an integer identifying which process is currently  *** 
  // ***  executing.                                                 *** 
  // ***                                                             *** 
  // ******************************************************************* 
  
  int id() {
    return _processID;
  }
  
  // ******************************************************************* 
  // ***         Public Method: nprocs                               *** 
  // ***                                                             *** 
  // ***  Returns an integer identifying the current number of       *** 
  // ***  active processes.                                          *** 
  // ***                                                             *** 
  // ******************************************************************* 
  
  int nprocs() {
    return _processCount;
  }
  
  
  // ******************************************************************* 
  // ***         Public Method: setCommTimeout                       *** 
  // ***                                                             *** 
  // ***  Sets the number of seconds that a process will wait to     *** 
  // ***  receive data from another process before throwing an       *** 
  // ***  exception.                                                 *** 
  // ***                                                             *** 
  // ******************************************************************* 
  
  void setCommTimeout(unsigned int commTimeout) {
    _commTimeout = commTimeout;
  }
  
  // ******************************************************************* 
  // ***         Public Method: send                                 *** 
  // ***                                                             *** 
  // ***  This aynchronous method sends the data referenced bu       *** 
  // ***  "dataToSend" to "destProcessID".  The size of the data     *** 
  // ***  is obtained by looking at the type "T".                    *** 
  // ***                                                             *** 
  // ******************************************************************* 
  
  template<class T>
  void send(int destProcessID, string dataTag, T &dataToSend) {
    logSendRecv(destProcessID, dataTag,LOG_SR_SEND, LOG_SR_START);  
    vector<char> data(sizeof(T));
    for(unsigned int k=0; k<sizeof(T); k++) 
      data[k]=((char*)&dataToSend)[k];
    send_binary(destProcessID, dataTag, data);    
    // cout << _processID << "->" << destProcessID << " " << dataTag << " " << dataToSend << endl;
    logSendRecv(destProcessID, dataTag, LOG_SR_SEND, LOG_SR_SUCCESS);
  }
  
  
  // ******************************************************************* 
  // ***         Public Method: send                                 *** 
  // ***                                                             *** 
  // ***  This aynchronous method sends the data at location         *** 
  // ***  "pdataToSend" to "destProcessID".  The size of the data    *** 
  // ***  being sent is provided in the integer: "dataSize".         *** 
  // ***                                                             *** 
  // ******************************************************************* 
  
  template<class T>
  void send(int destProcessID, string dataTag, 
	    T *pdataToSend, mdp_int dataSize) {
    logSendRecv(destProcessID, dataTag, LOG_SR_SEND, LOG_SR_START);  
    vector<char> data(sizeof(T)*dataSize);
    for(mdp_int k=0; k<data.size(); k++) 
      data[k]=((char*)pdataToSend)[k];
    send_binary(destProcessID, dataTag, data);    
    logSendRecv(destProcessID, dataTag,LOG_SR_SEND, LOG_SR_SUCCESS);
  }

  // ******************************************************************* 
  // ***         Public Method: recv                                 *** 
  // ***                                                             *** 
  // ***  This synchronous "blocking" method retrieves the data sent *** 
  // ***  to "destProcessID".  The size of the data being sent is    *** 
  // ***  provided in the integer: "dataSize".                       *** 
  // ***                                                             *** 
  // ******************************************************************* 
  
  template<class T>
  void recv(int sourceProcessID, string dataTag, T &dataToReceive) {
    logSendRecv(sourceProcessID, dataTag,LOG_SR_RECV, LOG_SR_START);    
    vector<char> data;
    recv_binary(sourceProcessID, dataTag, data);    
    if(data.size()!=sizeof(T)) {
      log("PSIM ERROR: recv invalid data)");
      throw string("PSIM ERROR: recv invalid data)");
    };
    for(unsigned int k=0; k<sizeof(T); k++) 
      ((char*)&dataToReceive)[k]=data[k];
    // cout << _processID << "<-" << sourceProcessID << " " << dataTag << " " << dataToReceive << endl;
    logSendRecv(sourceProcessID, dataTag, LOG_SR_RECV, LOG_SR_SUCCESS);   
  }
  
  
  // ******************************************************************* 
  // ***         Public Method: recv                                 *** 
  // ***                                                             *** 
  // ***  This synchronous "blocking" method retrieves the data sent *** 
  // ***  to "destProcessID".  "dataSize" bytes are copied to        *** 
  // ***  location "pdataToReceive".                                 *** 
  // ***                                                             *** 
  // ******************************************************************* 

  template<class T>
  void recv(int sourceProcessID, string dataTag, 
	    T *pdataToReceive, mdp_int dataSize) {
    logSendRecv(sourceProcessID, dataTag, LOG_SR_RECV, LOG_SR_START);    
    vector<char> data;
    recv_binary(sourceProcessID, dataTag, data);    
    if(data.size()!=sizeof(T)*dataSize) {
      log("PSIM ERROR: recv invalid data size");
      throw string("PSIM ERROR: recv invalid data size");
    }
    for(mdp_int k=0; k<data.size(); k++) 
      ((char*)pdataToReceive)[k]=data[k];
    logSendRecv(sourceProcessID, dataTag, LOG_SR_RECV, LOG_SR_SUCCESS);   
  }
  
  // ******************************************************************* 
  // ***         Public Method: broadcast                            *** 
  // ***                                                             *** 
  // ***  Allows broadcasting data to all processes.                 *** 
  // ***  sourceProcessID sends data (data) to all of the other      *** 
  // ***  processes who receive the data through the same variable). *** 
  // ***                                                             *** 
  // ******************************************************************* 
  
  template<class T>
  void broadcast(int sourceProcessID, T &data) {
    static string tag = "BROADCAST:0";
    doBegEndLog(tag, LOG_BEGIN);    
    if (_processID == sourceProcessID) {
      for (int i=0; i<_processCount; i++) {
	if (i != sourceProcessID) {
	  send(i, tag, data);
	}
      }
    } else {
      recv(sourceProcessID, tag, data);
    }    
    if(tag=="BROADCAST:0") tag="BROADCAST:1"; else tag="BROADCAST:0";
    doBegEndLog(tag, LOG_END);
  }

  template<class T>
    void broadcast(int sourceProcessID, T *data, int dataSize) {
    static string tag = "BROADCASTV:0";
    doBegEndLog(tag, LOG_BEGIN);    
    if (_processID == sourceProcessID) {
      for (int i=0; i<_processCount; i++) {
	if (i != sourceProcessID)
	  send(i, tag, data, dataSize);
      }
    } else {
      recv(sourceProcessID, tag, data, dataSize);
    }    
    if(tag=="BROADCASTV:0") tag="BROADCASTV:1"; else tag="BROADCASTV:0";
    doBegEndLog(tag, LOG_END);
  }
  

  // ******************************************************************* 
  // ***         Public Method: collect                              *** 
  // ***                                                             *** 
  // *** All parallel processes construct a list of the data passed  *** 
  // *** by each process.  The list is broadcasted and returned by   *** 
  // *** each processor.  This method is used to implement global    *** 
  // *** sum and some other global operations.                       *** 
  // ***                                                             *** 
  // ******************************************************************* 
  
  template<class T>
  vector<T> collect(int dest, T &data) {
    static string tag="COLLECT";
    vector<T> dataList;
    T dataToReceive;
    dataList.resize(_processCount);
    doBegEndLog(tag, LOG_BEGIN);
    
    if (_processID != dest) {
      send(dest, tag, data);
    } else {
      dataList[dest]=data;
      
      for (int i=0; i<_processCount; i++) {
	if(i!=dest) {
	  recv(i, tag, dataToReceive);
	  dataList[i]=dataToReceive;
	}
      }
    }
    if(tag=="COLLECT:0") tag="COLLECT:1"; else tag="COLLECT:0";
    doBegEndLog(tag, LOG_END);    
    return dataList;
  }

  // ******************************************************************* 
  // ***         Public Method: combine                              *** 
  // ***                                                             *** 
  // *** All parallel processes construct a list of the data passed  *** 
  // *** by each process.  The list is broadcasted and returned by   *** 
  // *** each processor.  This method is used to implement global    *** 
  // *** sum and some other global operations.                       *** 
  // ***                                                             *** 
  // ******************************************************************* 

  template<class T>
  vector<T> combine(T &data) {
    vector<T> dataList=collect(PROCESS_PARENT,data); 
    cout << id() << " size=" << dataList.size() << endl;
    
    broadcast(PROCESS_PARENT, &dataList[0], dataList.size());    
    cout << id() << " list=" << dataList[0] << dataList[1]<< dataList[2]<< endl;
    return dataList;
  }
  
  // ******************************************************************* 
  // ***         Public Method: barrier                              *** 
  // ***                                                             *** 
  // ***  Initiates a blocking point so that the processes pause     *** 
  // ***  until ALL processes have reached the barrier.              *** 
  // ***                                                             *** 
  // ******************************************************************* 
  
  void barrier() {
    int dummy;
    broadcast(PROCESS_PARENT,dummy);
    collect(PROCESS_PARENT,dummy);
  }
    
  // ******************************************************************* 
  // ***         Public Method: add                                  *** 
  // ***                                                             *** 
  // *** All parallel processes sum their data in parallel.  The sum *** 
  // *** is returned.                                                *** 
  // ***                                                             *** 
  // ******************************************************************* 
  
  template<class T>
  T add(T &item) {
    vector<T> dataList;
    T total=0;
    dataList=collect(PROCESS_PARENT,item);
    if(_processID==PROCESS_PARENT) 
      for (int i=0; i<dataList.size(); i++) {
	total += dataList[i];   
      }
    broadcast(PROCESS_PARENT,total);
    return total;
  }

  // ******************************************************************* 
  // ***         Public Method: auxiliary functions                  *** 
  // *** Examples:                                                   *** 
  // *** a.out -PSIM_NPROCS=2           (parallel processes)         *** 
  // *** a.out -PSIM_LOGFILE=./test.log (log into test.log)          *** 
  // *** a.out -PSIM_VERBATIM=1         (show all communications)    *** 
  // ***                                                             *** 
  // ******************************************************************* 

  static int parse_argv_nprocs(int argc, char** argv) {
    int n=1;
    for(int i=1; i<argc; i++)
      if(strncmp(argv[i],"-PSIM_NPROCS=",13)==0) {
	sscanf(argv[i]+13,"%i",&n);
	break;
      }
    return n;
  }

  static string parse_argv_logfile(int argc, char** argv) {
    for(int i=1; i<argc; i++)
      if(strncmp(argv[i],"-PSIM_LOGFILE=",14)==0) {
	return string(argv[i]+14);
      }
    return string("");
  }

  static int parse_argv_verbatim(int argc, char** argv) {
    int n=1;
    for(int i=1; i<argc; i++)
      if(strncmp(argv[i],"-PSIM_VERBATIM=",15)==0) {
	sscanf(argv[i]+15,"%i",&n);
	break;
      }
    return n;
  }
};


