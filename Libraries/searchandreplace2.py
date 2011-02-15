import os
import sys
import string

for filename in os.listdir('./'):
    if filename=='mdp_compatibility_macros.h': continue
    if filename[-2:]!='.h': continue
    print 'file ', filename
    choice='y' #raw_input(' replace (y/n)? ')
    if choice.upper() in ['Y', 'YES']:
        os.system ('cp %s %s.bak' % (filename,filename))
        try:
            file=open(filename,'r')
        except:
            continue
        s=''
        for line in file.readlines():
            if line[:4]=="/// " and len(line)>5 and line[4:5] in string.lowercase:
                c=line[4:5]
                print '< ', line
                line2=line[:4]+c.upper()+line[5:]
                print '> ', line2
                c=raw_input("replace?").upper()
                if c=='Y':
                    line=line2
                elif c=='Q':
                    sys.exit(0)
            s=s+line
            pass
        file.close()
        file=open(filename,'w')
        file.write(s)
    
            
    
