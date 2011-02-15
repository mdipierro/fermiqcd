import sys, re, os

regex=re.compile('(?P<key>\S+)\s*=\s*(?P<value>\S+)\s*')
start=False
data={}
counter=0
filename=sys.argv[1]
for line in open(filename,'rb'):
    counter+=len(line)
    if not start and line[:12]=='BEGIN_HEADER': start=True
    elif not start: raise SyntaxError, "Invalid file format"
    elif line[:10]=='END_HEADER': break
    else:
        match=regex.match(line)
        print line.strip()
        if not match: raise SyntaxError, "Invalid file format"
        key,value=match.group('key'), match.group('value')
        data[key]=value

print 'header size:', counter
print data

if data['DATATYPE'].find('3x3')>0:
    os.system('./nersc2mdp.exe -skip %s -gauge %sx%sx%sx%s:3DY %s %s.mdp' % (counter,
       data['DIMENSION_4'],data['DIMENSION_1'],data['DIMENSION_2'],data['DIMENSION_3'],
        filename,filename))
else:
    os.system('./nersc2mdp.exe -skip %s -gauge %sx%sx%sx:2%sDY %s %s.mdp' % (counter,
       data['DIMENSION_4'],data['DIMENSION_1'],data['DIMENSION_2'],data['DIMENSION_3'],
        filename,filename))

os.system('./plaquette.exe %sx%s %s.mdp' % (data['DIMENSION_4'],data['DIMENSION_1'],filename))
