import sys, re, os

allowed_commands={'cold':[['TxXxYxZ','\d+x\d+x\d+x\d+','16x4x4x4'],['nc','\d+','3']],
	 'hot':[['TxXxYxZ','\d+x\d+x\d+x\d+','16x4x4x4'],['nc','\d+','3']],
	 'load':[['filename','[a-zA-Z0-9_/.\-]*',None],['precision','(?:float)|(?:double)','float']],	 
         'lload':[['filename','[a-zA-Z0-9_/.\-\*]*','gauge.*.mdp'],['precision','(?:float)|(?:double)','float']],	 
         'save_partitioning':[['filename','[a-zA-Z0-9_/.\-]*','partitining']],
	 'heatbath':[['beta','\d+(.\d+)+','5.0'],['steps','\d+','1']],
         'save':[['filename','[a-zA-Z0-9_/.\-\*]*','gauge.*.mdp']],         
         'plaquette':[],
         'ape_smear':[['alpha','\d(.\d+)+','0.7'],['steps','\d+','20'],['cooling_steps','\d+','10']],
         'coulomb_gauge_fix':[['boost','\d(.\d+)+','1.0'],['steps','\d+','20'],['precision','\d+(.\d+)+(e\-\d+)*','1e-6'], ['z3','\d',0]],
         'landau_gauge_fix':[['boost','\d(.\d+)+','1.0'],['steps','\d+','20'],['precision','\d+(.\d+)+(e\-\d+)*','1e-6']],
         'topological_charge':[['filename','[a-zA-Z0-9_/.\-]*','topological_charge_*.vtk'],['t','\d','-1']],         
        }

comments={'cold': 'creates a cold gauge configuration, i.e. all links U(x,mu)=1\nexmaple: "-cold nc=3 TxXxYxZ=10x4x4x4"',
          'hot': 'creates a hot gauge configuration, i.e.  all links U(x,mu)=random.SU(nc)\nexmaple: "-hot nc=3 TxXxYxZ=10x4x4x4"',
          'load': 'loads an existing gauge configuration in mdp fermiqcd format\ndetermines lattice size from file',
          'lload': 'lloads in a loop an existing gauge configuration in mdp fermiqcd format\ndetermines lattice size from file',          
	  'save_partitioning':'save a vtk file showing parallel partitioning information',          
          'heatbath': 'performs a number of heatbath steps using the WilsonGaugeAction',
          'save': 'saves current gauge configuration',
          'plaquette':'compute average plquette',
          'ape_smear': 'APE smearing as defiend in ref...',
          'coulomb_gauge_fix': 'Coulomb Gauge Fixing algorithm as defiend in ref...\n(can also do z3 fixing if z3=1 as in ref...)',
          'landau_gauge_fix': 'Landau Gauge Fixing algorithm as defiend in ref...',
          'topological_charge': 'saves a vtk showing the topological charge as defined in ref...\n(one should ape_smear first)',
}

class Command:
    def __init__(self,command):
        self.command=command
        self.args={}
    def __repr__(self):
        return '%s[%s]' % (self.command,repr(self.args))

NC_SWITCH="""   int nc=0;          
   switch(header.bytes_per_site) {
      case 4*4*1: precision=4; nc=1; break; 
      case 8*4*1:  precision=8; nc=1; break; 
      case 4*4*4:  precision=4; nc=2; break; 
      case 8*4*4:  precision=8; nc=2; break; 
      case 4*4*9:  precision=4; nc=3; break; 
      case 8*4*9:  precision=8; nc=3; break; 
      case 4*4*16:  precision=4; nc=4; break; 
      case 8*4*16:  precision=8; nc=4; break; 
      case 4*4*25:  precision=4; nc=5; break; 
      case 8*4*25:  precision=8; nc=5; break; 
      case 4*4*36:  precision=4; nc=6; break; 
      case 8*4*36:  precision=8; nc=6; break; 
      case 4*4*49:  precision=4; nc=7; break; 
      case 8*4*49:  precision=8; nc=7; break; 
      case 4*4*64:  precision=4; nc=8; break; 
      case 8*4*64:  precision=8; nc=8; break; 
      case 4*4*81:  precision=4; nc=9; break; 
      case 8*4*81:  precision=8; nc=9; break; 
      case 4*4*100:  precision=4; nc=10; break; 
      case 8*4*100:  precision=8; nc=10; break;
   }
"""

def parse(argv=sys.argv[2:]):
    instruction=' '.join(sys.argv[2:])
    commands=[]
    loops=[]

    k=0
    errors=[]
    warnings=[]
    for s in argv:
        if s[0]=='+':        
            commands.append(Command(s[1:]))
        elif s=='{':
            loops.append(k-1)
            commands[-1].args['begin']=k
        elif s=='}':
            k1=loops[-1]
            del loops[-1]	
            commands[k1].args['end']=k
        else: 
            try: 
                a,b=s.split('=')
   	        commands[-1].args[a]=b;
            except:
                errors.append('unable to parse %s' % s)
        k=len(commands)
    k=0
    if not commands[0].command in ['cold','hot','load']:
        errors.append('the first algorithm must be -cold, -hot, or -load')
    for s in commands:
        if s.command=='loop':
            if not s.args.has_key('end'):
                errors.append('+loop error: missing }')
	    if not s.args.has_key('n'):
	        errors.append('+loop error: missing n=... argument')
        elif not allowed_commands.has_key(s.command):
            errors.append('+%s error: unkown algorithm' % s.command)
        else:
            d={}
            arguments=allowed_commands[s.command]
            for a,b,c in arguments:                
                d[a]=1
                if not s.args.has_key(a):
                    if c==None:
                        errors.append('+%s error: missing %s=... argument' % (s.command,a))
                    else:
                        s.args[a]=c
                        warnings.append('+%s warning: assuming default argument %s=%s' % (s.command,a,str(c)))
                elif re.sub(b,'',s.args[a])!='':  #fix this line
                    errors.append('+%s error: invalid argument %s=%s' % (s.command,a,s.args[a]))            
            for a in s.args:
                if not d.has_key(a):
                    errors.append('+%s error: unkown argument %s=...' % (s.command,a))
        if s.command in ['cold', 'hot'] and k>0:
            errors.append('+%s error: must be the first algorithm!' % (s.command)) 
        k+=1            
    return instruction, commands, errors, warnings

def report(errors, warnings):
    if errors:
        print 'YOU HAVE THE FOLLLOWING ERRORS:'
        for e in errors: print e
    if warnings:
        print 'OK BUT SOME WARNINGS:'
        for w in warnings: print w    
         
def generate_code(instruction,warnings,commands):
    program=''
    program+='/*\n'
    program+='    python qcd.py code '+instruction+'\n'
    program+='\n'
    for w in warnings: program+='    '+w+'\n'
    program+='*/\n'
    program+='#include "fermiqcd.h"\n\n'
    program+='int main(int argc, char** argv) {\n'
    program+='   mdp.open_wormholes(argc,argv);\n'
    program+='   string filename;\n'
    program+='   coefficients coeff;\n'

    have_gauge=False   
    indent=1
    loops=[]
    k=0
    for s in commands:
      SPACE='\n'+'   '*indent
      if not have_gauge:
        if s.command=='cold':
            program+='   int L[]={%s};\n' % s.args['TxXxYxZ'].replace('x',',')
	    program+='   mdp_lattice spacetime(4,L);\n'
            program+='   int nc=%s;\n' % s.args['nc']
            program+='   gauge_field U(spacetime,nc);\n'
            program+='   set_cold(U);\n'
            have_gauge=True
        if s.command=='hot':
            program+='   int L[]={%s};\n' % s.args['TxXxYxZ'].replace('x',',')
	    program+='   mdp_lattice spacetime(4,L);\n'
            program+='   int nc=%s;\n' % s.args['nc']
            program+='   gauge_field U(spacetime,nc);\n'
            program+='   set_hot(U);\n'
            have_gauge=True
        if s.command=='load':
            program+='   mdp_field_file_header header=get_info(%s);\n' % s.args['filename']
            program+='   int L[]={header.box[0],header.box[1],header.box[2],header.box[3]};\n'
            program+=NC_SWITCH
	    program+='   mdp_lattice spacetime(4,L);\n'
            program+='   gauge_field U(spacetime,nc);\n'
            if s.args['precision']=='float': program+='   U.load_as_float("%s");\n' % s.args['filename']
            if s.args['precision']=='double': program+='   U.load_as_double("%s");\n' % s.args['filename']
            have_gauge=True
      else:
        if s.command=='loop':
            j=len(loops)
            loops.append(s.args['end'])
            program+=SPACE+'for(int i%i=0; i%i<%s; i%i++) {\n' % (j,j,s.args['n'],j)
            indent+=1
        if s.command=='load':                      
            if s.args['precision']=='float': program+=SPACE+'U.load_as_float("%s");\n' % s.args['filename']
            if s.args['precision']=='double': program+=SPACE+'U.load_as_double("%s");\n' % s.args['filename']
        if s.command=='lload':                      
            if s.args['precision']=='float': program+=SPACE+'U.load_as_float(glob("%s")[i%i]);' % (s.args['filename'],len(loops)-1)
            if s.args['precision']=='double': program+=SPACE+'U.load_as_double(glob("%s")[i%i]);' % (s.args['filename'],len(loops)-1)
        if s.command=='save_partitioning':                      
            program+=SPACE+'save_partitioning_vtk(spacetime,"%s");' % (s.args['filename']) #NOT IMPLEMENTED
        if s.command=='save':                      
            program+=SPACE+'U.save("%s");' % s.args['filename']
	if s.command=='plaquette':
            program+=SPACE+'mdp << "plaquette=" << average_plaquette(U) << endl;'
	if s.command=='heatbath':
            program+=SPACE+'coeff["beta"]=%s;' % s.args['beta']
            program+=SPACE+'WilsonGaugeAction::heatbath(U,coeff,%s);' % s.args['steps']
	if s.command=='ape_smear':
            program+=SPACE+'ApeSmearing::smear(U,%s,%s,%s);' % (s.args['alpha'],s.args['steps'],s.args['cooling_steps'])
	if s.command=='coulomb_gauge_fix':
            program+=SPACE+'GaugeFixing::fix(U,GaugeFixing::Coulomb,%s,%s,%s,%s);' % (s.args['steps'],s.args['precision'],s.args['boost'],s.args['z3'])
	if s.command=='landau_gauge_fix':
            program+=SPACE+'GaugeFixing::fix(U,GaugeFixing::Landau,%s,%s,%s);' % (s.args['steps'],s.args['precision'],s.args['boost'])
        if s.command=='topological_charge':                      
            program+=SPACE+'{float tc=topological_charge_vtk(U,"%s",%s);' % (s.args['filename'],s.args['t']) # NOT IMPLEMENTED
            program+=SPACE+'mdp << "topological_charge=" << tc << endl; }'
      k+=1
      if loops and k==loops[-1]:
        del loops[-1]
        indent-=1
        program+='\n'+'   '*indent+'}\n'
    program+='\n\n'
    program+='   mdp.close_wormholes();\n'
    program+='   return 0;\n'
    program+='}\n'
    return program

def menu():
    if len(sys.argv)==1:
	print 'usage:'
	print '"python qcd.py help" to list all availbale algorithms'
	print '"python qcd.py help -heatbath" for help about the heatbath algorithm'
	print '"python qcd.py code +cold TxXxYxZ=10x4x4x4 nc=3" to create a program that makes one cold gauge configuration'
	print '"python qcd.py compile filename.cpp" to compile filename.cpp (for those who hate make like me)'
        # print '"python qcd.py submit filename.exe" to submit filename.exe to your cluster
        # print '"python qcd.py monitor pid" to monitor the executable running as pid
    elif sys.argv[1]=='code':
	instruction,commands,errors,warnings=parse()  
	if errors or warnings: 
	   report(errors,warnings)
	if not errors:
	   program=generate_code(instruction,warnings,commands)
	   open('latest.cpp','w').write(program)
	   #os.system('g++ test.cpp -I../Libraries -DOSX')       
	   print 'saved as latest.cpp'
    elif sys.argv[1]=='compile':
	try: filename=sys.argv[2]
        except: filename='latest.cpp'
	osx=raw_input('do you have an an OS X 10 Mac (y/n)? ').lower()
	sse=raw_input('do you have a processor that supports SSE2 (y/n)? ').lower()
	pos=raw_input('do you have FULL posix support (y/n)? ').lower()
	par=raw_input('do you want to compile with MPI support (y/n)? ').lower()
	command='g++ -O2 %s -I../Libraries -o %s' % (filename,filename.replace('.cpp','.exe'))
	if osx=='y': command+=' -DOSX'
	elif sse=='y': command+=' -DSSE2'
	if par=='y': command+=' -DPARALLEL'
	if pos=='n': command+=' -DNO_POSIX'
	print command
	os.system(command)
    elif sys.argv[1]=='help' and len(sys.argv)<3:
	print 'available algorithms:'
	keys=allowed_commands.keys()
	keys.sort()
	for key in keys:
	    print '   -'+key
    elif sys.argv[1]=='help' and len(sys.argv)==3:
	key=sys.argv[2][1:]
	if not allowed_commands.has_key(key):
	    print 'command -%s not supported' % key
	else:
	    print 'INFO FOR ALGORITHM -%s' % key
	    print 'description:'
	    print comments[key]
	    print 'arguments:'
	    for a,b,c in allowed_commands[key]:
		print '    "%s"' % a,
		if c==None: print 'is required'
		else: print 'is %s by default' % c
    else:
	print 'I do not understand you commands'

if __name__=='__main__': menu()
         
