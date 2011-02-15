import os

sin=raw_input('pattern to replace: ')
sout=raw_input('replace with: ')
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
            if line.find(sin)>=0:
                print '< ', line
                line=line.replace(sin,sout)
                print '> ', line
            s=s+line
            pass
        file.close()
        file=open(filename,'w')
        file.write(s)
    
            
    
