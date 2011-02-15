import copy

PLUS=+1
MINUS=-1

def WilsonGaugeAction():
    gauge=[]
    for mu in range(4):
	for nu in range(mu):
	    path=[[PLUS,mu],[PLUS,nu],[MINUS,mu],[MINUS,nu]]
	    gauge.append(path)
    return gauge

def IzubuchiGaugeAction():
    gauge=WilsonGaugeAction()
    for mu in range(4):
	for nu in range(mu):
	    path=[[PLUS,mu],[PLUS,mu],[PLUS,nu],
		  [MINUS,mu],[MINUS,mu],[MINUS,nu]]
	    gauge.append(path)
	    path=[[PLUS,nu],[PLUS,nu],[PLUS,mu],
		  [MINUS,nu],[MINUS,nu],[MINUS,mu]]
	    gauge.append(path)
    return gauge

def reverse_path(path):
    new_path=[]
    for item in eval(path):
	new_path.insert(0,[-item[0],item[1]])
    return str(new_path)

def process_gauge(gauge_paths,ndim=4):
    derivate_paths={}    
    for mu in range(ndim):
	derivate_paths[mu]=new_paths=[]
	for path in gauge_paths:
	    k=0
	    for dir,nu in path:
		if nu==mu and dir==PLUS:
		    new_paths.append(str(path[k+1:]+path[:k])) 
		if nu==mu and dir==MINUS:
		    new_paths.append(reverse_path(str(path[k+1:]+path[:k]))) 
		k+=1
    return derivate_paths

def contains(txt1,txt2):
    if txt1.find(txt2[1:-1])>-1: return True
    if txt1.find(reverse_path(txt2)[1:-1])>-1: return True
    return False

def optimize_gauge_paths(derivate_paths):
    str_paths=[]
    for key, paths in derivate_paths.items(): 
	str_paths+=paths
    for j in range(1,len(str_paths)):
	for i in range(j):
	    if contains(str_paths[i],str_paths[j]):
		str_paths[j],str_paths[i]=str_paths[i],str_paths[j]
    base=[]
    names={}
    codes={}
    for i in range(0,len(str_paths)):
	for j in range(0,len(str_paths)):
	    if i!=j and contains(str_paths[i],str_paths[j]):
		#print i,'contains',j, str_paths[j]
		if i<j: raise Exception
		name=str_paths[j]
		if not names.has_key(name): 
		    code=names[name]=str(len(base)+1)
		    codes[code]=name
		    base.append(name)
    base2=copy.copy(base)
    base2.reverse()
    for name in base:
	code=names[name]
	reversed_name=reverse_path(name)
	reversed_code='-'+str(code)
	for key, paths in derivate_paths.items():
	    for i in range(len(paths)):
		paths[i]=paths[i].replace(name[1:-1],code)
		paths[i]=paths[i].replace(reversed_name[1:-1],reversed_code)
    return base,names,codes,derivate_paths

def write_mul1(name,path):
    shift=[]
    if path[0][0]==PLUS: 
	s='y+%i' % path[0][1]
	print 'm_set(tmp(x,%s),U(y,%i)); y=%s;' % (name,path[0][1],s)	
    if path[0][0]==MINUS: 
	s='y-%i' % path[0][1]
	print ' y=%s; m_seth(tmp(x,%s),U(y,%i));' % (s,name,path[0][1])
    for dir,mu in path[1:]:
	if dir==PLUS:
	    s='y+%i' % mu
	    print 'm_mul(tmp(x,%s),U(y,%i)); y=%s;' % (name,mu,s)
	if dir==MINUS:
	    s='y-%i' % mu
	    print 'y=%s; m_mulh(tmp(x,%s),U(y,-1,%i));' % (s,name,mu)

def shift(path):
    d={}
    for x in eval(path):
	if not d.has_key(x[1]): d[x[1]]=0
	if x[0]>0: 
	    d[x[1]]+=1
	else:
	    d[x[1]]-=1
    txt='y'
    for key,value in d.items():
	while value!=0:
	    if value>0:
		txt='(%s+%i)' % (txt,key)
		value-=1
	    if value<0:
		txt='(%s-%i)' % (txt,key)
		value+=1    
    return txt

def reverse_shift(path):
    path=[[-x[0],x[1]] for x in eval(path)]
    path.reverse()
    return shift(str(path))

def write_add_mul1(mu,path,codes):
    SPACE='   '
    if type(path[0])==type([]):
	if path[0][0]==PLUS: 
	    print SPACE+'m_store(A,U(y,%i)); y=y+%i;' % (path[0][1],path[0][1])
	if path[0][0]==MINUS: 
	    print SPACE+'y=y-%i; m_storeh(A,U(y,%i));' % (path[0][1],path[0][1])
    elif path[0]>0:
	print SPACE+'m_set(A,tmp(y,%s)); y=%s;' % (path[0], shift(codes[str(path[0])]))
    elif path[0]<0:
	print SPACE+'y=%s; m_seth(A,tmp(y,%s));' % (reverse_shift(codes[str(-path[0])]),path[0])
    for item in path[1:]:
	if type(item)==type([]):
	    dir,mu=item
	    if dir==PLUS:
		print SPACE+'m_mul(A,U(y,%i)); y=y+%i;' % (mu,mu)
	    if dir==MINUS:
		print SPACE+'y=y-%i; m_mulh(A,U(y,%i);' % (mu,mu)
	elif item>0:
		print SPACE+'m_mul(A,tmp(y,%s)); y=%s;' % (item,shift(codes[str(item)]))
	elif item<0:
		print SPACE+'y=%s; m_mulh(A,tmp(y,%s));' % (reverse_shift(codes[str(-item)]),item)
    print SPACE+'m_add(staple,A);'

def write_multiplications(base,names,codes,derivate_paths,ndim=4):
    print 'mdp_nmatrix_field tmp(U.lattice(),%i,U.nc,U.nc);' % len(base)
    for item in base:
	write_mul1(names[item],eval(item))
    for mu in range(ndim):
	print 'if(mu==%i) {' % mu
	for path in derivate_paths[mu]:
	    write_add_mul1(mu,eval(path),codes)
	print '}'

b,n,c,d=optimize_gauge_paths(process_gauge(IzubuchiGaugeAction()))
write_multiplications(b,n,c,d)

