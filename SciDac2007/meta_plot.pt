from math import *
from random import *
import datetime
import sys
import thread
try: 
    from Tkinter import *
    import tkFileDialog
    import Pmw
except: pass


class Plot:
    def __init__(self,root, width=600, height=400, xlabel='X', ylabel='Y', bg='black', fg='white', name="myCanvas"):        
        self.width=width
        self.height=height
        self.minr=3
        self.maxr=5
        self.root=root
        self.initialized=False
        self.parameters=w0,dw,h0,dh=50,100,50,80
        self.legend_counter=0
        self.xlabel=xlabel
        self.ylabel=ylabel
        self.bg=bg
        self.fg=fg
        self.name=name
        self.width, self.height=width, height
        self.startdate=None #in case x axis is a date
        self.selection_rectangle=None # for zooming
        self.sets_plot=[]
        self.sets_hist=[]
        if root!='js':
            self.root.resizable(False, False)
            self.frame=Frame(root,width=width,height=height)        
            self.canvas=Canvas(self.frame,width=width,height=height, bg=bg)
            self.canvas.pack()
            self.canvas.create_text(self.width/2,self.height/2,text="Loading...",fill=fg)
            self.balloon = Pmw.Balloon(root)
            self.button=Button(self.root,text='save',command=self.save)
        else:
            self.canvas=None
            self.js=""

    def selectNW(self,e):
        a,b,c,d,minx,maxx,miny,maxy=self.abcd()        
        self.x0=max(minx,min(maxx,(e.x-a)/b))
        self.y0=max(miny,min(maxy,(e.y-c)/d))
        
    def selectSE(self,e):
        a,b,c,d,minx,maxx,miny,maxy=self.abcd()        
        self.x1=max(minx,min(maxx,(e.x-a)/b))
        self.y1=max(miny,min(maxy,(e.y-c)/d))
        self.show()        

    def moveNSEW(self,e):
        a,b,c,d,minx,maxx,miny,maxy=self.abcd()        
        self.x1=max(minx,min(maxx,(e.x-a)/b))
        self.y1=max(miny,min(maxy,(e.y-c)/d))
        if self.selection_rectangle:
            self.canvas.delete(self.selection_rectangle)
        self.selection_rectangle=self.canvas.create_rectangle(a+b*self.x0,c+d*self.y0,a+b*self.x1,c+d*self.y1,fill=None,outline=self.fg)
        self.canvas.update()

    def show(self):
        if self.root!='js':
            return self.show_tk()
        else:
            return self.show_js()
            
    def show_tk(self):
        self.canvas.delete(ALL)
        self.canvas.bind('<Button-1>', lambda e, s=self: s.selectNW(e))
        self.canvas.bind('<B1-Motion>', lambda e, s=self: s.moveNSEW(e))
        self.canvas.bind('<ButtonRelease-1>', lambda e, s=self: s.selectSE(e))
        self.canvas.create_window(self.width,self.height,
                                  anchor='se',window=self.button)
        self.canvas.create_rectangle(0,0,self.width,self.height,fill=self.bg)
        self.plot_axis()
        for legend,color,data,width in self.sets_hist:
            self.hist_(data,color,legend,width)
        for legend,color,data in self.sets_plot:
            self.plot_(data,color,legend)

    def show_js(self):
        self.js="""
        <div id="%s" style="position:relative;height:%ipx;width:%ipx;">
        </div> 
        <script type="text/javascript">
        <!--
        // http://www.walterzorn.com/jsgraphics/jsgraphics_e.htm#download
        var jg = new jsGraphics("%s");
        jg.setPrintable(false);
        jg.setFont("arial","12px");
        """ % (self.name,self.height,self.width,self.name)
        if self.bg:
            self.js+="""
            jg.setColor("%s");
            jg.fillRect(0,0,%i,%i);
            """ % (self.bg,self.width,self.height)
        self.plot_axis()
        for legend,color,data,width in self.sets_hist:
            self.hist_(data,color,legend,width)
        for legend,color,data in self.sets_plot:
            self.plot_(data,color,legend)
        self.js+="""
        jg.paint();      
        //-->
        </script>
        """
        return self.js
    
    def pack(self,**ka):
        self.frame.pack(ka)
        
    def clear():
        self.canvas.delete(ALL)    

    def round(self,minx,maxx,n):
        dx=float(maxx-minx)/n        
        if dx<=0: raise SyntaxError
        p=10.0**int(log(dx)/log(10)-0.99)
        dx=p*int(dx/p+1)
        return minx,dx,maxx

    def find_axis(self,data):                
        if not self.initialized:
            minx,miny=maxx,maxy=data[0][0],data[0][1]
            if len(data[0])>2: miny=miny+data[0][2]
            if len(data[0])>3: maxy=maxy+data[0][3]
            self.initialized=True
        else:
            minx,maxx,miny,maxy=self.minx,self.maxx,self.miny,self.maxy
        for i in range(0,len(data)):
            x,y=data[i][0],data[i][1]
            if x<minx: minx=float(x)
            if x>maxx: maxx=float(x)
            if len(data[i])>2: y=data[i][1]+data[i][2]
            if y<miny: miny=float(y)
            if len(data[i])>3: y=data[i][1]+data[i][3]
            if y>maxy: maxy=float(y)        
        minx,dx,maxx=self.round(minx,maxx,10)
        miny,dy,maxy=self.round(miny,maxy,10)
        self.minx=minx
        self.miny=miny
        self.maxx=maxx
        self.maxy=maxy
        self.dx=dx
        self.dy=dy

    def find_axis_hist(self,data):                
        w=(data[1][0]-data[0][0])/2
        if not self.initialized:
            minx,miny=data[0][0]-w,data[0][1]
            maxx,maxy=data[0][0]+w,data[0][1]
            if len(data[0])>2: miny=miny+data[0][2]
            if len(data[0])>3: maxy=maxy+data[0][3]
            self.initialized=True
        else:
            minx,maxx,miny,maxy=self.minx,self.maxx,self.miny,self.maxy
        for i in range(0,len(data)):
            x,y=data[i][0],data[i][1]
            if x-w<minx: minx=float(x)-w
            if x+w>maxx: maxx=float(x)+w
            if len(data[i])>2: y=data[i][1]+data[i][2]
            if y<miny: miny=float(y)
            if len(data[i])>3: y=data[i][1]+data[i][3]
            if y>maxy: maxy=float(y)        
        minx,dx,maxx=self.round(minx,maxx,10)
        miny,dy,maxy=self.round(miny,maxy,10)
        self.minx=minx
        self.miny=miny
        self.maxx=maxx
        self.maxy=maxy
        self.dx=dx
        self.dy=dy

    def abcd(self):
        if self.minx<0: minx=-self.dx*int(-self.minx/self.dx+0.99)
        else: minx=self.dx*int(self.minx/self.dx-0.99)        
        maxx=self.dx*int(self.maxx/self.dx+0.99)
        if self.miny<0: miny=-self.dy*int(-self.miny/self.dy+0.99)
        else: miny=self.dy*int(self.miny/self.dy-0.99)        
        maxy=self.dy*int(self.maxy/self.dy+0.99)

        w0,dw,h0,dh=self.parameters
        b=float(self.width-dw)/(maxx-minx)
        a=w0-b*minx
        d=-float(self.height-dh)/(maxy-miny)
        c=self.height-h0-d*miny
        return a,b,c,d,minx,maxx,miny,maxy
        
    def plot_axis(self):
        try:
            if self.x1!=self.x0:
                self.minx,self.maxx=min(self.x0,self.x1),max(self.x0,self.x1)
            if self.y1!=self.y0:
                self.miny,self.maxy=min(self.y0,self.y1),max(self.y0,self.y1)
            if self.x1==self.x0 and self.y1==self.y0: raise Exception
            self.minx,self.dx,self.maxx=self.round(self.minx,self.maxx,10)
            self.miny,self.dy,self.maxy=self.round(self.miny,self.maxy,10)
            if self.startdate and self.dx<1: self.dx=1.0
        except:
            for legend,color,data in self.sets_plot: self.find_axis(data)
            for legend,color,data,width in self.sets_hist: self.find_axis_hist(data)
            pass        
        self.legend_counter=0
        minx=self.minx
        miny=self.miny
        maxx=self.maxx
        maxy=self.maxy
        dx=self.dx
        dy=self.dy
        a,b,c,d,minx,maxx,miny,maxy=self.abcd()
        #plot axis
        self.create_line(a+b*minx,c+d*miny,a+b*minx,c+d*maxy,width=2, fill=self.fg)
        self.create_line(a+b*minx,c+d*miny,a+b*maxx,c+d*miny,width=2, fill=self.fg)
        # plot y tics
        x0,y0=minx,miny
        while x0<=maxx+dx/2:
            self.create_line(a+b*x0,c+d*y0,a+b*x0,c+d*y0-5,width=2, fill=self.fg)
            if not self.startdate:
                self.create_text(a+b*x0,c+d*y0+5,text='%g' % x0,anchor='n', fill=self.fg)
            else:
                date=(self.startdate+datetime.timedelta(days=int(x0))).strftime('%d%b\n%Y')
                self.create_text(a+b*x0,c+d*y0+5,text=date,anchor='n', fill=self.fg)                
            x0=x0+dx
            if abs(x0)<1e-12: x0=0.0
            
        #plot x tics
        x0,y0=minx,miny
        while y0<=maxy+dy/2:
            self.create_line(a+b*x0,c+d*y0,a+b*x0+5,c+d*y0,width=2, fill=self.fg)
            self.create_text(a+b*x0-5,c+d*y0,text='%g' % y0,anchor='e', fill=self.fg)
            y0=y0+dy
            if abs(y0)<1e-12: y0=0.0

        self.create_text(a+b*minx,c+d*maxy-20,text=self.ylabel,anchor='c', fill=self.fg)
        self.create_text(self.width/2,self.height-20,text=self.xlabel,anchor='c', fill=self.fg)

    def pick_color(self):
        if self.legend_counter==0: return '#6699FF'
        if self.legend_counter==1: return '#FF9966'
        if self.legend_counter==2: return '#66FF99'
        if self.legend_counter==3: return '#99FF66'
        if self.legend_counter==4: return '#FF6699'
        if self.legend_counter==5: return '#9966FF'
        return self.fg

    def add_legend(self,color,legend):
        if legend:
            w0,dw,h0,dh=self.parameters
            rx=ry=3
            x0,y0=w0+10,h0+20*self.legend_counter
            self.create_oval(x0-rx,y0-ry,x0+rx,y0+ry,fill=color)
            self.create_text(x0+rx+5,y0,anchor='w',text=legend,fill=color)
            self.legend_counter+=1

    def plot_(self,data,color=None,legend=''):
        if not data: return
        if not color: color=self.pick_color()
        w0,dw,h0,dh=self.parameters
        a,b,c,d,minx,maxx,miny,maxy=self.abcd()
        #plot lines if no error bars
        if len(data[0])==2:
            x0,y0=data[0][0],data[0][1]
            for i in range(len(data)):
                x1,y1=data[i][0],data[i][1]
                self.create_line(a+b*x0,c+d*y0,a+b*x1,c+d*y1,fill=color,width=2)
                x0,y0=x1,y1
        #plot points (with error bars)
        if len(data[0])>2:
            rx=ry=min(max(self.minr,int(self.width/20/len(data))),self.maxr)
            for i in range(len(data)):
                x0,y0=data[i][0],data[i][1]
                yl=yh=y0
                if len(data[i])>2: yl=y0+data[i][2]
                if len(data[i])>3: yh=y0+data[i][3]            
                id=self.create_oval(a+b*x0-rx,c+d*y0-ry,a+b*x0+rx,c+d*y0+ry,fill=color)
                if len(data[i])>4:
                    text=data[i][4]
                    self.balloon_tagbind(self.canvas,id,text)
                self.create_line(a+b*x0,c+d*yl,a+b*x0,c+d*yh,fill=color,width=1)
                pass
            pass
        self.add_legend(color,legend)

    def markx(self,x0,color='white',legend=''):
        if not color: color=self.pick_color()
        a,b,c,d,minx,maxx,miny,maxy=self.abcd()
        self.create_line(a+b*x0,c+d*self.miny,a+b*x0,c+d*self.maxy,fill=color,width=2)
        self.add_legend(color,legend)
        
    def rectangle(self,x0,y0,x1,y1,color='white',legend='',text=''):
        if not color: color=self.pick_color()
        a,b,c,d,minx,maxx,miny,maxy=self.abcd()
        id=self.create_rectangle(a+b*x0,c+d*y0,a+b*x1,c+d*y1,fill=color)
        self.balloon_tagbind(self.canvas,id,text)
        self.add_legend(color,legend)

    def hist_(self,data,color=None,legend='',width=1.0):
        if not data: return
        w=(data[1][0]-data[0][0])*width/2
        if not color: color=self.pick_color()
        for x0,y0 in data:
            x0=float(x0)
            text='source:%s\nbin:%f\nheight:%f' % (legend,x0,y0)
            self.rectangle(x0-w,0,x0+w,y0,color=color,text=text)
        self.add_legend(color,legend)

    def plot(self,data,color=None,legend=''):
        self.sets_plot.append((legend,color,data))

    def hist(self,data,color=None,legend='',width=1.0):
        self.sets_hist.append((legend,color,data,width))

    def line(self,x0,y0,x1,y1,color='white',legend=''):
        a,b,c,d,minx,maxx,miny,maxy=self.abcd()
        self.create_line(a+b*x0,c+d*y0,a+b*x1,c+d*y1,fill=color,width=2)
        self.add_legend(color,legend)

    def balloon_tagbind(self,canvas,id,text):
        if canvas:
            self.balloon.tagbind(canvas,id,text)
        else:
            x0=id[2]+5
            y1=id[3]
            self.create_text(x0,y1,text,'w',id[4])

    def create_line(self,x0,y0,x1,y1,width,fill):
        if self.root!='js':
            self.canvas.create_line(x0,y0,x1,y1,width=width,fill=fill)
        else:
            self.js+="""
            jg.setColor("%s");
            jg.setStroke(%i);
            jg.drawLine(%i,%i,%i,%i);
            """ % (fill,width,x0,y0,x1,y1)

    def create_rectangle(self,x0,y0,x1,y1,fill):
        if self.root!='js':
            return self.canvas.create_rectangle(x0,y0,x1,y1,fill=fill)
        else:
            dx,dy=x1-x0,y1-y0
            self.js+="""
            jg.setColor("%s");
            jg.fillRectangle(%i,%i,%i,%i);
            """ % (fill,x0-dx/2,y0-dy/2,dx,dy)
            return (x0,y0,x1,y1,fill)

    def create_oval(self,x0,y0,x1,y1,fill):
        if self.root!='js':
            return self.canvas.create_oval(x0,y0,x1,y1,fill=fill)
        else:
            dx,dy=x1-x0,y1-y0
            self.js+="""
            jg.setColor("%s");
            jg.fillEllipse(%i,%i,%i,%i);
            """ % (fill,x0,y0,dx,dy)
            return (x0,y0,x1,y1,fill)
        
    def create_text(self,x0,y0,text,anchor,fill):
        if self.root!='js':
            self.canvas.create_text(x0,y0,text=text,anchor=anchor, fill=fill)
        else:
            if anchor=='n' or anchor=='c':
                y0,al=y0,"center"
                self.js+="""
                jg.setColor("%s");
                jg.drawStringRect("%s",%i,%i,%i,"%s");
                """ % (fill,text,x0-80,y0,160,al)            
            elif anchor=='e':
                y0,al=y0-8,"left"
                self.js+="""
                jg.setColor("%s");
                jg.drawStringRect("%s",%i,%i,%i,"%s");
                """ % (fill,text,x0-40,y0,40,al)            
            else:
                y0,al=y0-8,"left"                
                self.js+="""
                jg.setColor("%s");
                jg.drawString("%s",%i,%i);
                """ % (fill,text,x0,y0)                        
                
    def postscript(self):
        return self.canvas.postscript()

    def save(self,file=None):
        if not file:
            file= tkFileDialog.asksaveasfile(parent=self.root,mode='w',
                                             title='Save as PostScript',
                                             filetypes=[("all files","*.ps")])
        if file == None: return
        try:
            file.write(self.postscript())
        except:
            pass

        """"
    WORK IN PROGRESS
    def grub(self):
        self.canvas.update()
        x0 = self.canvas.winfo_rootx()
        y0 = self.canvas.winfo_rooty()
        x1 = x0 + self.canvas.winfo_width()
        y1 = y0 + self.canvas.winfo_height()
        image = ImageGrab.grab((x0, y0, x1, y1))
        return image        
    """
    
def bin(set,nbins=30,minimum=None,maximum=None):    
    bins=[0]*nbins
    mini,maxi=min(set),max(set)
    if minimum!=None: mini=minimum
    if maximum!=None: maxi=maximum
    step=(maxi-mini)/nbins/0.999
    c=0
    for item in set:
        i=int((item-mini)/step)
        if 0<=i<nbins:        
            bins[i]+=1
            c+=1
    for i in range(nbins):
        bins[i]=(mini+step*(float(i)+0.5),float(bins[i])/c/step)
    return bins

def chi_squared(set_exp,set_th):
    n=len(set_exp)
    if n!=len(set_th): raise SyntaxError
    sum=0.0
    for i in range(n):
        if abs(set_exp[i][0]-set_th[i][0])>0.01: raise SyntaxError
        sum+=((set_exp[i][1]-set_th[i][1])**2)/set_th[i][1]
    return sum

def test_plot():
    data=[]
    for i in range(1000):
        data.append((i/10,gauss(1,1),0,0))
    root=Tk()
    p=Plot(root,xlabel='X axis')
    p.plot([(5,2,-1,1),(10,2,-1,1),(20,3,-1,1),(40,4,-1,1)],legend='set1')
    p.plot(data,legend='set2')
    p.show()
    p.pack()
    #p.markx(70,legend='break')
    #p.rectangle(80,0,90,1,legend='break')
    #print p.postscript()
    root.mainloop()

def test_gauss():
    mu=20
    sigma=5    
    set=[]
    for i in range(100000): set.append(gauss(mu,sigma))
    data1=bin(set,20,0,40)
    data2=[]
    f=lambda x: 1.0/sqrt(2.0*3.1415)/sigma*exp(-((x-mu)**2)/2.0/sigma**2)
    for x,y in data1: data2.append((x,f(x)))
    print chi_squared(data1,data2)
    root=Tk()
    p=Plot(root,xlabel='# BIN')
    p.hist(data1,legend='experiment')
    p.plot(data2,legend='theory')
    p.show()
    p.pack()
    root.mainloop()
        

def test_stock():
    data=[(0,10)]
    for i in range(30):
        data.append((i+1,data[i][1]*(1.0+gauss(0.20,0.1))))
    root=Tk()
    p=Plot(root,xlabel='X axis')
    p.hist(data,legend='AMD')
    p.plot(data,color='blue')
    p.show()
    p.pack()
    root.mainloop()

        
if __name__=='__main__': test_plot()
