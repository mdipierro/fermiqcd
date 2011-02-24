from optparse import *
from rpy import *
import re, urllib, csv

usage = "python iplot.py\n" \

version = "iplotv1.0\n" \
          "  Copyright (c) 2007 Massimo Di Pierro\n" \
	  "  All rights reserved\n" \
          "  License: GPL 2.0\n\n" \
	  "  Written by Massimo Di Pierro <mdipierro@cs.depaul.edu>\n"

description = "plot the output of ibootstrap.py"

r.library('Hmisc')

def clean(text):
    return re.sub('\s+','',text.replace('/','_div_'))

class IPlot:
    def __init__(self,filename,plot_type,items=[]):
        self.type=plot_type
        if self.type=='quartz': r.quartz()	
	self.plot_raw_data(filename+'_raw_data.csv')
        self.plot_autocorrelations(filename+'_autocorrelations.csv')
        self.plot_trails(filename+'_trails.csv')
	self.plot_samples(filename+'_samples.csv')
	self.plot_min_mean_max(filename+'_min_mean_max.csv',items)

    def begin(self,filename):
        if self.type=='ps': r.postscript(filename+'.ps')        
        if self.type=='png': r.png(filename+'.png')        

    def end(self):
        if self.type=='ps': r.dev_off()
        else: raw_input('press enter to continue')


    def plot_raw_data(self,filename):
	for items in csv.reader(open(filename,'r'),delimiter=',',quoting=csv.QUOTE_NONNUMERIC):
	   tag=items[0]
	   data=items[1:]
	   self.begin(filename[:-4]+'_%s' % clean(tag))
	   r.plot(x=range(len(data)),y=data,xlab='step',ylab=tag,main='')
	   self.end()
	   self.begin(filename[:-4]+'_%s_hist' % clean(tag))
	   r.hist(data,n=len(data)/20,xlab=tag,ylab='frequency',main='',prob='T')
	   r.rug(data) 
	   self.end()
	   #self.begin(filename[:-4]+'_%s_qq.ps' % clean(tag))
	   #r.qqnorm(data,xlab=tag+' quantiles',main='')
	   #r.qqline(data)
	   #self.end()
	   mu=r.mean(data)
	   sd=r.sd(data)
	   probs=[min(x,1-x) for x in [r.pnorm((x-mu)/sd) for x in data]]
	   self.begin(filename[:-4]+'_%s_probability' % clean(tag))
	   r.plot(range(len(probs)),probs,xlab="step",ylab="probability "+tag,main="",type="p")
	   self.end()

    def plot_autocorrelations(self,filename):
	for items in csv.reader(open(filename,'r'),delimiter=',',quoting=csv.QUOTE_NONNUMERIC):
	   tag=items[0]
	   data=items[1:]
	   self.begin(filename[:-4]+'_%s' % clean(tag))
	   r.plot(x=range(len(data)),y=data,xlab='step',ylab=tag,main='')
	   self.end()

    def plot_trails(self,filename):
	for items in csv.reader(open(filename,'r'),delimiter=',',quoting=csv.QUOTE_NONNUMERIC):
	   tag=items[0]
	   data=items[1:]
	   self.begin(filename[:-4]+'_%s' % clean(tag))
	   r.plot(x=range(len(data)),y=data,xlab='step',ylab=tag,main='',type='p')
	   self.end()

    def plot_samples(self,filename):
	for items in csv.reader(open(filename,'r'),delimiter=',',quoting=csv.QUOTE_NONNUMERIC):
	   tag=items[0]
	   data=items[1:]
	   self.begin(filename[:-4]+'_%s_hist' % clean(tag))
	   r.hist(data,n=len(data)/10,xlab=tag,ylab='frequency',main='',prob='T')
	   r.rug(data) 
	   self.end()
	   #self.begin(filename[:-4]+'_%s_qq.ps' % clean(tag))
	   #r.qqnorm(data,xlab=tag+' quantiles',main='')
	   #r.qqline(data)
	   #self.end()

    def plot_min_mean_max(self,filename,xlab=['t']):
	lines=list(csv.reader(open(filename,'r'),delimiter=',',quoting=csv.QUOTE_NONNUMERIC))
	tags=lines[0] 
	if not xlab or xlab[0]=='': xlab=[tags[1]]
	index=0
	for tag in tags:
	    if tag[0]=='[': break
	sets={}
	for items in lines[1:]:
	   tag=items[0]
	   data=items[1:]       
	   legend=""
	   for i in range(1,len(tags)-3):
	       if not tags[i] in xlab: legend+="%s=%g " % (tags[i],data[i-1])
	   if not sets.has_key(legend):
	       x,y,yminus,yplus=[],[],[],[]
	       sets[legend]=(x,y,yminus,yplus)
	   else:
	       x,y,yminus,yplus=sets[legend]
	   t=data[index]
	   x.append(t)
	   y.append(data[-2])
	   yminus.append(data[-3])
	   yplus.append(data[-1])
	#v=r.FALSE	
	for legend in sets.keys():
	   x,y,yminus,yplus=sets[legend]
	   self.begin(filename[:-4]+'_%s' % clean(legend))
	   r.errbar(x,y,yminus,yplus,xlab=tags[index+1],ylab=tags[0],main="")
	   #v=r.TRUE
	   self.end()
  
def shell_iplot():
    parser=OptionParser(usage,None,Option,version)
    parser.description=description
    parser.add_option('-i','--input_prefix',default='ibootstrap',dest='input_prefix',
		      help='the prefix used to build input filenames')
    parser.add_option('-p','--plot_type',default='ps',dest='plot_type',
		      help='ps or quartz')
    parser.add_option('-v','--plot_variables',default='',dest='plot_variables',
		      help='plotting variables')
    parser.add_option('-f','--fit',default=[],dest='fits',
		      action='append',
		      help='fits to be performs on results')
    (options, args) = parser.parse_args()
    if options.fits:
	print 'sorry -f not implemented yet!'
    plot=IPlot(options.input_prefix,options.plot_type,options.plot_variables.split(','))

if __name__=='__main__': shell_iplot()
