from math import *
import re, random, copy, sys, csv
from optparse import *

usage = \
    "ifit.py [OPTIONS] 'expression@values'\n" \
    "  Example: ifit 'a*x+b@a=3,b=0'\n" \
    "  default filename is ibootstrap_min_mean_max.csv\n" \
    "  ...., 'x', 'min', 'mean', 'max'\n" \
    "  ...., 23, 10, 11, 12\n" \
    "  ...., etc etc etc\n"

version = \
    "ifit v1.0\n" \
    "  Copyright (c) 2007 Massimo Di Pierro\n" \
    "  All rights reserved\n" \
    "  License: GPL 2.0\n\n" \
    "  Written by Massimo Di Pierro <mdipierro@cs.depaul.edu>\n"

description = \
    "This program takes data produced by ibootstrap and fits it\n" \
    "it also does correlated fits by using the built-in function\n" \
    "(a==b)"

def clean(text):
    return re.sub('\s+','',text.replace('/','_div_'))


def invert_squared_matrix(A,checkpoint=None):
    """Computes the inverse of A using Gauss-Jordan emilimination"""
    A=copy.deepcopy(A)
    n=len(A)
    B=[[0.0 for i in range(n)] for j in range(n)]
    for i in range(n): B[i][i]=1
    for c in range(n):
        if checkpoint: checkpoint('pivoting (%i) ...' % c)
        for r in range(c+1,n):
            if abs(A[r][c])>abs(A[c][c]):
                A[r],A[c],B[c],B[r]=A[c],A[r],B[r],B[c]
                pass
            pass
        p=float(A[c][c])
        for k in range(n):
            A[c][k],B[c][k]=float(A[c][k])/p,float(B[c][k])/p
            pass
        for r in range(0,c)+range(c+1,n):
            p=float(A[r][c])
            for k in range(n):
                A[r][k]-=p*A[c][k]
                B[r][k]-=p*B[c][k]
                pass
            pass
        pass
    return B
    
class IFitException(Exception): 
    """
    Exception thrown by class IFit, usually due to wrong user input
    """
    pass

def restricted_eval(expression,loc):
    exec('__result__=%s' % expression) in loc
    return loc['__result__']

class IFit:
    """
    The class contains fitting algorithm
    """
    def __init__(self,expression,points,symbols=None,condition='True',import_module=None):
	"""
	Example:

	>>> ifit=IFit("a*x+b",[(0,1,0.1),(1,2,0.1)],['x'],'True')
        >>> ifit.fit(a=0.0,b=0.0)
        {'a':1.0,'b':1.0}
	"""
        if not symbols and len(points[0])==4: symbols=['x']
        if not symbols and len(points[0])==5: symbols=['x','y']
        if not symbols and len(points[0])==6: symbols=['x','y','z']
        if not symbols and len(points[0])==7: symbols=['x','y','z','t']
        if not symbols: raise Exception
        self.expression=expression  ### expression to use for fit
        self.points=points          ### raw data (x,y,dy)
        self.symbols=symbols        ### ['x']
        self.condition=condition    ### 'True' if all points
        self.locals={}              ### the restricted environment
        self.delta=0.001            ### for derivatives
        self.precision=1e-6         ### for convergence
        self.errors=[]              ### errors       
        self.variables_samples=[]   ### in iterative mode
	self.scatter_point=0        ### do not do the scatter plot
	self.callback=None
	self.last_fit=None
	self.last_least_squares=None
	self.last_variables=None
	self.last_hessian=None
	self.last_trail=None
        exec('from math import *') in self.locals
        if import_module: 
	    exec('from %s import *' % import_module) in self.locals
        ### parsing expression
        e=re.compile('[a-zA-Z_]+\w*').findall(expression)
        ### symbols like 'x' that are in the data
        self.required_symbols=[k for k in e if k in symbols]
        ### symbols like 'a' and 'b' that are intended to be variables
        self.undefined_symbols=[k for k in e if k and not self.locals.has_key(k) and not k in symbols]

    def apply(self,**objects):
	"""
	stores the arguments in self.locals
	"""
        for key,values in objects.items():
            self.locals[key]=value

    def f(self,**variables):
	"""
	evaluates the expression after adding the arguments to self.locals
	"""
        try:
            for key,value in variables.items():
                self.locals[key]=value
            return restricted_eval(self.expression,self.locals)
        except:
            self.errors.append('unable to evaluate the expression')
            raise IFitException

    def least_squares(self,**variables):
	"""
	the quantity we want to minimize as function of the variables
	"""
	self.last_fit=[]
        least_squares=0.0 
        for p in self.points:
            try:
                for i in range(len(self.symbols)):
                   self.locals[self.symbols[i]]=p[i]
                   if not restricted_eval(self.condition,self.locals):
		       continue
            except:
                self.errors.append('unable to evaluate the condition')
                raise IFitException
            ye=self.f(**variables)
            yo=p[-2]
	    if ye>yo:
		err=p[-1]-yo
	    else:
		err=yo-p[-3]
            least_squares+=((ye-yo)/err)**2
	    off=(ye-yo)*2/(p[-1]-p[-3])
            self.last_fit.append(p+[ye,off])	
	for key,value in variables.items():
	    if self.baesyan.has_key(key):
	        least_squares+=((value-self.priors[key])/self.baesyan[key])**2
        return least_squares
    def move(self,**variables):        
	"""
	makes a moves in the space of parameters, in the direction 
	tangent to the point
	"""
        n=len(variables)
        rn=range(n)
        delta=self.delta
	two_delta=2.0*delta
        df=[0.0]*n
        ddf={}
        hessian=[[0.0]*n for i in rn]
        keys=variables.keys()
        keys.sort()
        i=0
        for key1 in keys:
	    variables[key1]-=delta
	    f_minus=self.least_squares(**variables)
            for key2 in keys[i:]:
                variables[key2]-=delta
	        f_minus_minus=self.least_squares(**variables)
                variables[key2]+=two_delta
	        f_minus_plus=self.least_squares(**variables)
                ddf[key2]=f_minus_plus-f_minus_minus
                variables[key2]-=delta
	    variables[key1]+=two_delta
	    f_plus=self.least_squares(**variables)
            j=i
            for key2 in keys[i:]:
                variables[key2]-=delta
	        f_plus_minus=self.least_squares(**variables)
                variables[key2]+=two_delta
	        f_plus_plus=self.least_squares(**variables)
                h=((f_plus_plus-f_plus_minus)-ddf[key2])/two_delta**2
                hessian[i][j]=hessian[j][i]=h
                variables[key2]-=delta
                j+=1
	    variables[key1]-=delta
            df[i]=(f_plus-f_minus)/(two_delta)            
            i+=1
	H=invert_squared_matrix(hessian)                
	sum_shift=0.0
        i=0
        for key1 in keys:
	    for j in rn:
                shift=H[i][j]*df[j] 
                variables[key1]-=shift
                sum_shift+=shift**2
            i+=1
        shift=sqrt(sum_shift/n)
	least_squares=self.least_squares(**variables)
        return variables, shift, least_squares, hessian

    def fit(self,**variables):
	"""
	fit the points that meet the condition with the expression
	variables define the starting point int the parameter space
	"""
        self.variables_samples=None
        self.apriori,self.baesyan={},{}
        for key in self.undefined_symbols:
           if not variables.has_key(key):             
               self.errors.append('variable %s must be initialized' % key)
               raise IFitException
        for key,value in variables.items():
	   if key[0]=='_':
               del variables[key]
               self.baesyan[key[1:]]=value      
        self.priors=copy.deepcopy(variables)
	self.last_trail=[]
        for k in range(100):
	    v=copy.copy(variables)
	    variables,shift,least_squares,hessian=self.move(**variables)
	    v['[least_squares]']=least_squares
	    self.last_trail.append(v)
	    if self.callback: self.callback()
	    if shift<self.precision: break
	self.last_variables=variables
	self.last_least_squares=least_squares
	self.last_hessian=hessian
        return variables, least_squares, hessian

    def extrapolate(self,**coordinates):
	"""
	assuming a fit has been done... extrapolate to the point at coordinates
	"""
	### attentions! this knows variables because they are in self.locals
        try:
            for key in self.required_symbols:
                self.locals[key]=coordinates[key]
	except:
            self.errors.append('undefined symbol %s in extrapolation' % key)
            raise IFitException
        return self.f()
    
    def extrapolate_with_errors(self,**coordinates):
	"""
	assuming and iterative_fit has been done...
	extrapolate to the point and compute the error on the extrapolaiton
	"""
        p=0.158
        if not self.variables_samples: 
            self.errors.append('requires iterative_fit_first')
            raise IFitException
        try:
            for key in self.required_symbols:
                self.locals[key]=coordinates[key]
	except:
            self.errors.append('undefined symbol %s in extrapolation' % key)
            raise IFitException
        r=[self.f(**items[0]) for items in self.variables_samples]
        r.sort()
	mean=self.extrapolate(**coordinates)
        return r[int(p*len(r))],mean,r[int((1.0-p)*len(r))]

    def iterative_fit(self,filename='ifit_scatter.csv',**variables):
        """
	repeats the fit on multiple datasets generated by shifting the 
	original points according to their sigma.
	it also produces a table of results that can be used to
	make a scatter plot of results and study the error in the
	variables/parameters
	"""
        variables_samples=[]
        points=self.points
        for i in range(self.scatter_points):
            kb=copy.deepcopy(variables)
	    self.points=[]
            for p in points:
                p=list(p)
                p[-2]=random.gauss(p[-2],p[-1])
                self.points.append(p)	  
	    items=self.fit(**kb)
	    print i,items[0]
            variables_samples.append(items)
        self.variables_samples=variables_samples	
        writer=csv.writer(open(filename,'w'),delimiter=',',quoting=csv.QUOTE_NONNUMERIC)
	writer.writerows([item[0].values() for item in variables_samples])
        return variables_samples

    def plot2d(self,key1='a',key2='b',scale1=0.1,scale2=0.1):
        """
	under development
	plot the least_squared around the minimum
	assuming a quadratic approximation given by the hessian
	"""
	import pylab
	key1,key2=key2,key1
	keys=self.last_variables.keys()
	keys.sort()
        nv=len(self.last_variables)
	for k in range(len(keys)):
	    if keys[k]==key1: i,vi=k, self.last_variables[key1]
	    if keys[k]==key2: j,vj=k, self.last_variables[key2]
	ax=pylab.arange(vi*(1.0-10*scale1),vi*(1.0+12*scale1),vi*scale1)[:21]
	ay=pylab.arange(vj*(1.0-10*scale2),vj*(1.0+12*scale2),vj*scale2)[:21]

        dv=[0.0]*nv
        z=[[0 for ik in ax] for jk in ay]
        i0=0
        for x in ax:
            j0=0
            for y in ay:
               dv[i],dv[j]=x-vi,y-vj
               least_squares=0.0
               for ik in range(nv):
                   for jk in range(nv): 
                       least_squares+=0.5*dv[ik]*dv[jk]*self.last_hessian[ik][jk]
               z[i0][j0]=least_squares-self.last_least_squares
               j0+=1
            i0+=1
        pylab.contour(z,extent=(min(ay),max(ay),min(ax),max(ax)))
        pylab.show()        

    def save_fit_trail(self,filename):
	keys=[k for k in self.last_trail[0].keys() if not k=='[least_squares]']
	keys.sort()
	keys.append('[least_squares]')
	writer=csv.writer(open(filename,'w'),delimiter=',',quoting=csv.QUOTE_NONNUMERIC)
	writer.writerow(keys)
	for item in self.last_trail:
	    writer.writerow([item[key] for key in keys])

    def save_fit(self,filename):
	"""
	under development
	"""
	writer=csv.writer(open(filename,'w'),delimiter=',',quoting=csv.QUOTE_NONNUMERIC)
	expression=self.expression+'@'+','.join(["%s=%g" % (k,v) for k,v in self.last_variables.items()])
	others=["[min]", "[mean]", "[max]","[%s]" % expression,"[error]"]
	writer.writerow(self.symbols+others)
	for p in self.last_fit:
	    writer.writerow(p)

def read_min_mean_max_file(filename):
    """
    reads a standard ibootstrap_min_mean_max.csv file,
    extract all the points and symetrizes the error bars
    """
    reader=csv.reader(open(filename,'r'),delimiter=',',quoting=csv.QUOTE_NONNUMERIC)
    lines=list(reader)
    symbols=lines[0][:-3]
    for i in range(len(symbols)):
	if symbols[i]=='[min]': 
	    symbols=symbols[:i]
	    break
    i+=3
    points=lines[1:]    
    return symbols,points

def test_ifit():
    print 'generating points with z=x*sin(y)+4*y and dz=1'
    points=[[x,y,x*sin(y)+4*y-1,x*sin(y)+4*y,x*sin(y)+4*y+1] for x in range(3) for y in range(100)]
    print 'fitting with a*sin(y)+b*y'
    print IFit("a*sin(y)+b*y",points,symbols=['x','y']).fit(a=0.0,b=0.0)

def test_correlated_ifit():
    print 'generating points with z=x*sin(y)+4*y and dz=1'
    points=[[x,y,x*sin(y)+4*y-1,x*sin(y)+4*y,x*sin(y)+4*y+1] for x in range(3) for y in range(100)]
    print 'fitting with (a0*(x==0)+a1*(x==1)+a2*(x==2))*sin(y)+c*y'
    ifit=IFit("(a0*(x==0)+a1*(x==1)+a2*(x==2))*sin(y)+b*y",points,symbols=['x','y'])
    print ifit.fit(a0=0.0,a1=0.0,a2=0.0,b=0.0)

def main_ifit():
    loc={}
    parser = OptionParser(usage, None, Option, version)
    parser.add_option("-c", "--condition",
		      type="string", dest="condition",
		      default="True",
		      help="sets a filter on the points to be fitted")
    parser.add_option("-s", "--scater_points",
		      type="string", dest="scatter_points",
		      default='0',
		      help="repeats the fit multiple times")
    parser.add_option("-p", "--plot",
		      dest="plot",type='string',
		      default='',
		      help="plots the hessian")
    parser.add_option("-t", "--test",
		      dest="test",action='store_true',
		      default=False,
		      help="test a fit")
    parser.add_option("-e", "--extrapolate",
		      type='string',dest="extrapolations",
		      default=[],action='append',
		      help="extrpolation point")
    parser.add_option("-i", "--input_prefix",
		      type='string',dest="input_prefix",
		      default='ibootstrap',
		      help="prefix used to build input filename [prefix]_min_mean_max.csv")
    options,args=parser.parse_args()
    if options.test:
	test_ifit()
	return
    try:
	filename=options.input_prefix+'_min_mean_max.csv'
	expression,initial=args[0].split('@')
	symbols,points=read_min_mean_max_file(filename)
	ifit=IFit(expression,points,symbols,condition=options.condition)
	variables=restricted_eval('dict(%s)' % initial,loc)
	variables,least_squares,hessian=ifit.fit(**variables)    
	ifit.scatter_points=int(options.scatter_points)
	if ifit.scatter_points: 	
	    ifit.iterative_fit(options.input_prefix+'_scatter.csv',**variables)
	for key,value in variables.items():
	    print '%s = %g' % (key, value)
	print 'least_squares=',least_squares
	print 'hessian=',hessian        
	for item in options.extrapolations:
	    coordinates=restricted_eval('dict(%s)' % item,loc)	
	    try: e=ifit.extrapolate_with_errors(**coordinates) 
	    except: e=ifit.extrapolate(**coordinates) 
	    print 'extrapolation %s -> %s' % (item,str(e))
	if options.plot:
	    print 'attention! plotting is under development'
	    key1,key2=options.plot.split(',')
	    ifit.plot2d(key1,key2)
	ifit.save_fit(options.input_prefix+'_fit_%s.csv' % clean(args[0]))
	ifit.save_fit_trail(options.input_prefix+'_fit_%s_trail.csv' % clean(args[0]))
    except IFitException:
	for item in ifit.errors: print item
	sys.exit(-1)

if __name__=='__main__': main_ifit()
