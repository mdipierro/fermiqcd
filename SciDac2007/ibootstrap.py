import sys, re, copy, random, glob, csv
from optparse import *
from math import *

usage = "ibootstrap.py *.log 'x[<a>]/y[<b>]' 'abs(a-b)==1'\n" \
    "  scans all files *.log for expressions of the form\n" \
    "    x[<a>]=<value> and y[<b>]=<value>\n" \
    "  and computes the average and bootstrap errors of x[<a>]/y[<b>]\n" \
    "  where <a> and <b> satisfy the condition abs(a-b)==1.\n"

version = "iboostrap v1.0\n" \
    "  Copyright (c) 2007 Massimo Di Pierro\n" \
    "  All rights reserved\n" \
    "  License: GPL 2.0\n\n" \
    "  Written by Massimo Di Pierro <mdipierro@cs.depaul.edu>\n"

description = "This is program to scan the log files of a Markov Chain Monte Carlo\n" \
    "Algorithm, parse for expressions and compute the avarge and bootstrap errors\n" \
    "of any function of those expressions." \
    "It also compute the convergence trails of the averages.\n" 
    
            
class IBootstrapException(Exception): pass

class IBootstrap:
    def __init__(self,filepattern="*.log",expression="x[<a>]",condition="a%2==0",
                 min_index=0,max_index=0,nsamples=100,percent=0.158,
		 output_prefix='ibootstrap', raw=False,advanced=False,import_module=None,indices=None):
        # indices is ignored by ibootstrap but has to be here 
        # becase of the design of iplotwx.py
        self.expression=re.sub('(\"|\<|\>)','',expression)
        if condition!='True': self.expression+=' where '+condition
	self.report=[]
        self.import_module=import_module
        self.indices=indices
        try:
	    self.parse_expression(expression,advanced)
            if not raw:
   	        self.parse_input(filepattern)
                self.log_raw_data(output_prefix+'_raw_data.csv')
		self.compute_autocorrelations()
		self.log_autocorrelations(output_prefix+'_autocorrelations.csv')
            else:
                self.reload_input(filepattern)
            self.find_expressions(expression,condition)
            if not max_index: max_index=self.length
            self.trails(min_index,max_index)
            self.bootstrap(min_index,max_index,nsamples,percent)          
            self.log_trails(output_prefix+'_trails.csv')
	    self.log_samples(output_prefix+'_samples.csv')
            self.log_min_mean_max(output_prefix+'_min_mean_max.csv')
	    self.status='success'        
    	except IBootstrapException:
            self.report.append('FATAL ERROR')        
            self.status='failure'

    def restricted_eval(self,expression,loc={}):
        try:  
            exec('__result__=(%s)' % str(expression)) in loc
            return loc['__result__']
        except:
            self.report.append('expression "%s" contains undefined variables' % expression)
            raise IBootstrapException

    def super_round(self,a,b,c):
        d=log(c-a)/log(10)-2
        if d<0: d=10**(-int(-d))
        else: d=10**int(d)
        if a<0: a=-d*int(-a/d)
        else: a=d*int(a/d)
        if b<0: b=-d*int(-b/d)
        else: b=d*int(b/d)
        if c<0: c=-d*int(-c/d)
        else: c=d*int(c/d)
        return a,b,c

    def make_regex(self,text):    
	for s in "\.%*+?()[]{}|": text=text.replace(s,'\\'+s)
	text=re.sub('<(?P<id>\w+)>','(?P<\g<id>>[-|+]?\d+)',text)
	return text

    def parse_expression(self,expression,advanced):
	items=re.compile('\"[^\"]+\"').findall(expression)
        self.items=items
        self.items_regex=[]
        self.items_pattern=[]
	for item in items:	
                if advanced: tag_pattern=item[1:-1]
		else: tag_pattern=self.make_regex(item[1:-1])
		float_pattern='[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?'
		regex=re.compile('(?P<__name__>%s)\s*(=|:|\s)?\s*(?P<__value__>%s)' % (tag_pattern,float_pattern))
                self.items_regex.append(regex)
		self.items_pattern.append(tag_pattern)

    def parse_input(self,filename):
        items=self.items
        items_regex=self.items_regex
        if not items:
            self.report.append("your expressions does not containes items in quotes")
            raise IBootstrapException
	symbols={}
	filenames=glob.glob(filename)
	if len(filenames)>1: self.report.append('found %i files' % len(filenames))
	elif len(filenames)==1: self.report.append('reading file %s' % filenames[0])
        else: 
            self.report.append('no input file matches the file pattern')
            raise IBootstrapException
	for name in filenames:
	    data=open(name,'r').read()
	    for i in range(len(items)):	
                item=items[i]
                regex=items_regex[i]
		groupindex=regex.groupindex
		nameindex=0
		valueindex=groupindex['__value__']-1
		occurrences=regex.findall(data)
		for occurrence in occurrences:
		    try: value=float(occurrence[valueindex])
		    except:
		       self.report.append('error: non-float value at %s=%s' % occurrence)
		       raise IBootstrapException
		    try: symbols[occurrence[nameindex]].append(value)
		    except: symbols[occurrence[nameindex]]=[value]        
	self.symbols=symbols

	length=-1
        keys=symbols.keys()
        keys.sort()
	for key in keys:     
	    ell=len(symbols[key])
	    self.report.append('%s occurs %i times' % (key, ell))
	    if length==-1: length=ell   
	    if ell!=length:
		self.report.append('warining, not all field with same occurrences')
	    if ell<length: length=ell
	for key in symbols.keys():     
	    symbols[key]=symbols[key][:length]
        self.length=length

    def find_expressions(self,expression,condition):
	items=self.items
        symbols=self.symbols
	items_pattern=self.items_pattern
	columns=[]
	items_regex=[]
	for tag_pattern in items_pattern:
	    keys=[]
	    regex=re.compile(tag_pattern)
	    items_regex.append(regex)
	    for key in symbols.keys():
		m=regex.match(key)
		if m: keys.append(m)
	    if not keys:
	        self.report.append('unable to match entire expression')
	        raise IBootstrapException
	    columns.append(keys)
	
	expressions=[]
        variables={}
	counters=[0]*len(items)
	if 0 in counters:
	     # DO SOMETHING
	     pass        
	done=False
	while 1:
	    e=copy.copy(expression)
	    isvalid=True
	    loc={}
	    for i in range(len(items)):
		m=columns[i][counters[i]]            
		e=e.replace(items[i],m.group(0),1)            
		for key in items_regex[i].groupindex.keys():
		    value=items_regex[i].groupindex[key]
		    v=float(m.group(value))
		    if loc.has_key(key) and loc[key]!=v: isvalid=False
		    else: loc[key]=v
	    variables[e]=copy.copy(loc)
	    exec('from math import *') in loc
            if self.import_module: exec('from %s import *' % self.import_module) in loc	    
            if not self.restricted_eval(condition,loc): isvalid=False
	    if loc['__result__']==False: isvalid=False
	    if isvalid:
		expressions.append(e)
            else:
                del variables[e]

	    counters[0]+=1
	    for i in range(len(counters)):
		if counters[i]>=len(columns[i]):
		    counters[i]=0
		    if i<len(counters)-1:
			counters[i+1]+=1
		    else:
			done=True
	    if done: break
	self.expressions=expressions
        self.variables=variables

    def compute_autocorrelations(self):
	### perhps this should compute Hurst exponent too
	symbols=self.symbols
	self.autocorrelations={}
	for key,sq in symbols.items():
	    ac=self.autocorrelations[key]=[]
	    n=len(sq)
	    mu=sum(sq)/n
	    sd2=sum([(k-mu)**2 for k in sq])/n	    
	    for i in range(0,n/2):
		auto=0.0
		for j in range(i,n):
		    auto+=(sq[j-i]-mu)*(sq[j]-mu)
		ac.append(auto/(n-i)/sd2)
	    self.report.append("autocorrelation for %s and d=1 is %g" % (key,ac[1]))


    def trails(self,min_index,max_index):
	""" under development """ 
        expressions=self.expressions
	symbols=self.symbols
        loc={}
        exec('from math import *') in loc
        if self.import_module: exec('from %s import *' % self.import_module) in loc	    
	trails={}
	avgs={}
        n=max_index-min_index        
	for key in symbols.keys():     
            sq=symbols[key]
            aq=avgs[key]=[sq[min_index]]
            for i in range(1,n):
                aq.append((i*aq[i-1]+sq[min_index+i])/(i+1))
	for expression in expressions:
	    trails[expression]=[0]*n
            for i in range(n):
                result=expression	    
    	        for key in avgs.keys():
		    result=result.replace(key,str(avgs[key][i]))                                        
		trails[expression][i]=self.restricted_eval(result,loc)
	self.trails=trails        

    def average(self,expressions,symbols,min_index,max_index):
	results={}
	avgs={}
	loc={}
        exec('from math import *') in loc
        if self.import_module: exec('from %s import *' % self.import_module) in loc	    
	for key in symbols.keys():        
	    avgs[key]=sum(symbols[key][min_index:max_index])/(max_index-min_index)
	for expression in expressions:
	    result=expression
	    for key in avgs.keys():
		result=result.replace(key,str(avgs[key]))                        
	    results[expression]=self.restricted_eval(result,loc)
	return results

    def bootstrap(self,min_index,max_index,nsamples=100,percent=0.158):
        expressions=self.expressions
        symbols=self.symbols
        if max_index>self.length: max_index=self.length
        if min_index>max_index: min_index,max_index=max_index,min_index
        if min_index<0: min_index=0 
        length=max_index-min_index
	samples={}
	means=results=self.average(expressions,symbols,min_index,max_index)
	for key in results.keys(): samples[key]=[]
	for sample in range(nsamples):
	    v=[random.randint(min_index,max_index-1) for i in range(length)]       
	    symbols2={}        
	    for key in symbols.keys():
		s=symbols[key]
		s2=symbols2[key]=[]            
		for i in v:
		    s2.append(s[i])                
	    results=self.average(expressions,symbols2,min_index,max_index)
	    for key in results.keys():
		 samples[key].append(results[key])
        min_mean_max={}
	for key in results.keys():
            samples[key].sort()      
            s=samples[key]
	    min_mean_max[key]=self.super_round(s[int(percent*nsamples)],means[key],s[int((1.0-percent)*nsamples)])
	keys=min_mean_max.keys()
        keys.sort()
        for key in keys:
            m=min_mean_max[key]
            self.report.append('< %s > = min: %g, mean: %g, max: %g' % (key,m[0],m[1],m[2]))
        self.samples=samples
        self.min_mean_max=min_mean_max

    def log_raw_data(self,filename='ibootstrap_raw_data.csv'):
        writer=csv.writer(open(filename,'w'),delimiter=',',quoting=csv.QUOTE_NONNUMERIC)
        keys=self.symbols.keys()
        keys.sort()
        for key in keys: writer.writerow([key]+self.symbols[key])
        self.report.append('raw data saved in %s' % filename)

    def log_autocorrelations(self,filename='ibootstrap_autocorrelations.csv'):
        writer=csv.writer(open(filename,'w'),delimiter=',',quoting=csv.QUOTE_NONNUMERIC)
        keys=self.autocorrelations.keys()
        keys.sort()
        for key in keys: writer.writerow([key]+self.autocorrelations[key])
        self.report.append('autocorrelations saved in %s' % filename)

    def reload_input(self,filename):
        reader=csv.reader(open(filename,'r'),delimiter=',',quoting=csv.QUOTE_NONNUMERIC)
        symbols={}
        length=-1
        for row in reader:
            symbols[row[0]]=data=row[1:]
            ell=len(data)
            if length==-1: length=ell
            elif ell<length: length=ell
        self.symbols=symbols
        self.length=length

    def log_samples(self,filename='ibootstrap_samples.csv'):
        writer=csv.writer(open(filename,'w'),delimiter=',',quoting=csv.QUOTE_NONNUMERIC)
        keys=self.samples.keys()
        keys.sort()
        for key in keys: writer.writerow([key]+self.samples[key])
        self.report.append('bootstrap samples saved in %s' % filename)

    def log_trails(self,filename='ibootstrap_trails.csv'):
        writer=csv.writer(open(filename,'w'),delimiter=',',quoting=csv.QUOTE_NONNUMERIC)
        keys=self.trails.keys()
        keys.sort()
        for key in keys: writer.writerow([key]+self.trails[key])
        self.report.append('average trails saved in %s' % filename)
    
    def log_min_mean_max(self,filename='ibootstrap_min_mean_max.csv'): 
        writer=csv.writer(open(filename,'w'),delimiter=',',quoting=csv.QUOTE_NONNUMERIC)
        keys=self.min_mean_max.keys()
        keys.sort()             
	variables=self.variables[keys[0]].keys()
	variables.sort()
	print self.variables[keys[0]]
	header=[self.expression]+variables+['[min]','[mean]','[max]']
	rows=[header]
	for key in keys:
	    rows.append([key]+[self.variables[key][var] for var in variables]+list(self.min_mean_max[key]))
	writer.writerows(rows)
        self.report.append('results saved in %s' % filename)

    def log_report(self,filename=None):
        if filename: file=open(filename,'w')
        else: file=sys.stdout
        for msg in self.report:
            file.write('%s\n' % msg)


def test_ibootstrap():
    file=open('test_samples.log','w')
    for i in range(100):
        for t in range(16):
            file.write('2pt[%.2i]: %f\n' % (t,(2.0+random.gauss(0,1))*exp(-0.2*t)))    
        for t in range(16):
            for t1 in range(16):
               file.write('3pt[%.2i][%.2i]: %f\n' % (t,t1,(4.0+random.gauss(0,1))*exp(-0.2*(t+t1))))
    file.close()
    IBootstrap('test_samples.log','"3pt[<t1>][<t2>]"/"2pt[<t1>]"/"2pt[<t2>]"','t1==t2',0,0,100).log_report()    
    return 0

def shell_ibootstrap():
    parser=OptionParser(usage,None,Option,version)
    parser.description=description
    parser.add_option('-b','--minimum_index',default='0',dest='min',
		      help='the first occurrence of expression to be considered')
    parser.add_option('-e','--maxmium_index',default='0',dest='max',
		      help='the last occurrence +1 of expression to be considered')
    parser.add_option('-n','--number_of_samples',default='100',dest='nsamples',
		      help='number of required bootstrap samples')
    parser.add_option('-p','--percentage',default='0.158',dest='percent',
		      help='percentage in the lower and upper tails')
    parser.add_option('-t','--test',action='store_true',dest='test',
		      default=False,help='make a test!')
    parser.add_option('-r','--raw',action='store_true',dest='raw',
		      default=False,help='Load raw data instead of parsing input')
    parser.add_option('-a','--advanced',action='store_true',dest='advanced',
		      default=False,help='In advanced mode use regular expressions for variable patterns')
    parser.add_option('-i','--import_module',dest='import_module',
		      default=None,help='import a python module for expression evaluation')
    parser.add_option('-o','--output_prefix',dest='output_prefix',
		      default='ibootstrap',help='path+prefix used to build output files')
    (options, args) = parser.parse_args()
    if options.test: return test_ibootstrap()
    if len(args)<2: return 1
    if len(args)<3: args.append('True')
    ibootstrap=IBootstrap(args[0],args[1],args[2],
			  int(options.min),int(options.max),  
			  int(options.nsamples),float(options.percent),
			  options.output_prefix,
                          options.raw,options.advanced,options.import_module)
    ibootstrap.log_report()
    if ibootstrap.status=='failure': return 1
    return 0

if __name__=='__main__':
    shell_ibootstrap()
