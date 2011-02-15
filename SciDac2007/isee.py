# Author: Antonio X. Garcia <antoniox.garcia@gmail.com>
# Modified by: Massimo Di Pierro <mdipierro@cs.depaul.edu>
# Copyright (c) 2007, DePaul University
# All rights reserved.
# License: BSD Style

import os
import wx
import sys
import glob
#from optparse import *
from enthought.mayavi.app import Mayavi
from enthought.mayavi.sources.vtk_file_reader import VTKFileReader
from enthought.mayavi.modules.outline import Outline
from enthought.mayavi.modules.text import Text
from enthought.mayavi.modules.iso_surface import IsoSurface
from enthought.mayavi.modules.surface import Surface
from enthought.mayavi.modules.orientation_axes import OrientationAxes

usage = "show.py FILE(s) [OPTION]"

version = "Show.py v1.0" \
          "\nCopyright (c) 2007, DePaul University" \
          "\nAll rights reserved." \
          "\nLicense: BSD Style" \
          "\n\nWritten by Antonio X. Garcia <antonioxgarcia@gmail.com>" \
          "\n\nand Massimo Di Pierro <mdipierro@cs.depaul.edu>"

def parse():
    args=sys.argv[1:]
    filepatterns=args[0]       
    files=[glob.glob(p) for p in filepatterns.split(',')]
    files=[glob.glob(p) for p in [filepatterns]]
    return files

class Watcher(wx.Timer):
    def __init__(self, interval, callable, *args, **kw_args):
        wx.Timer.__init__(self)
        self.callable = callable
        self.args = args
        self.kw_args = kw_args
        self.Start(interval)

    def Notify(self):
        self.callable(*self.args, **self.kw_args)
        

class MayaViShow(Mayavi):
    def __init__(self, filelists, title):
        self.filelists  = filelists
	self.title=title
	self.data=[]
        self.watcher = None        
	self.stop=True
        self.imagecount = 0        
        self.sceneTextCount = None
        self.filestatus=[]
        
    def run(self):
        win = self.script.engine.application.gui.GetTopWindow()
        win.CenterOnScreen()        
        self.script.new_scene()
        for idx in range(len(self.filelists)):
            f=self.filelists[idx][0]
            d = VTKFileReader()
            d.initialize(f)
            self.data.append(d)
            self.filestatus.append(os.stat(f)[-2])
            mayavi.add_source(d)
            self.addmodules(idx, idx==0)
        self.modifymenu()

            
    def modifymenu(self):
        window = mayavi.engine.application.gui.GetTopWindow()
        menubar = window.GetMenuBar()
        mnuVid = wx.Menu("")
        pos = menubar.GetMenuCount() - 1
        menubar.Insert( pos, mnuVid, "Animations")
        mnuAnimLoop = wx.MenuItem( mnuVid, wx.ID_ANY, "Loop")
        mnuVid.AppendItem( mnuAnimLoop )
        mnuAnimLoop.Enable(1)
        window.Bind(wx.EVT_MENU, self.AnimLoop, mnuAnimLoop)
        mnuAnimStop = wx.MenuItem( mnuVid, wx.ID_ANY, "Stop")
        mnuVid.AppendItem( mnuAnimStop )
        mnuAnimStop.Enable(1)
        window.Bind(wx.EVT_MENU, self.AnimStop, mnuAnimStop)
        mnuAnimWatch = wx.MenuItem( mnuVid, wx.ID_ANY, "Watch")
        mnuVid.AppendItem( mnuAnimWatch )
        mnuAnimWatch.Enable(1)
        window.Bind(wx.EVT_MENU, self.AnimWatch, mnuAnimWatch)
        mnuAnimMPEG = wx.MenuItem( mnuVid, wx.ID_ANY, "Save as MPEG")
        mnuVid.AppendItem( mnuAnimMPEG )
        mnuAnimMPEG.Enable(1)
        #window.Bind(wx.EVT_MENU, self.AnimMPEG, mnuAnimMPEG)

    def AnimStop(self,evt):
        self.stop=True

    def AnimLoop(self,evt):
        print 'AnimLoop', self.stop
        self.stop=False
        self.timestep=0
        mayavi.watcher=Watcher(1000,self.AnimLoopRec)
    
    def AnimLoopRec(self):
        print 'AnimLoopRec'
        if self.stop: 
	    self.watcher=None
            return
        filelist=self.filelists[0]
        for idx in range(len(self.filelists)):
            if self.timestep<len(self.filelists[idx]):
                self.data[idx].initialize(filelist[self.timestep])
        strNextVal = str(self.timestep)
        self.sceneTextCount.text = "          "
        self.sceneTextCount.text = strNextVal
        self.timestep+=1
        if self.timestep>=len(filelist): self.watcher=None

    def AnimWatch(self,evt):
        print 'AnimLoop', self.stop
        self.stop=False
        self.frame=0
        mayavi.watcher=Watcher(1000,self.AnimWatchRec)

    def AnimWatchRec(self):
        print 'AnimLoopRec'
        if self.stop: 
	    self.watcher=None
            return
        filelist=self.filelists[0]
        changed=False
        for idx in range(len(self.filelists)):
            k=0
            for f in self.filelists[idx][:1]:                
                if self.filestatus[idx]!=os.stat(f)[-2]:
                    self.data[idx].reader.modified()
                    self.data[idx].update()
                    self.data[idx].data_changed=True
                    changed=True
                k+=1
        if(changed):
            strNextVal = str(self.frame)
            self.sceneTextCount.text = "          "
            self.sceneTextCount.text = strNextVal
            self.frame+=1
                       
    def addmodules(self, idx, addTextIndexer):
        #add Outline module        
        o = Outline()
        mayavi.add_module(o)
        #add OrientationAxes module
        oa = OrientationAxes()
        mayavi.add_module(oa)

	"""
        #add IsoSurface module if requested
        i = IsoSurface()
        mayavi.add_module(i)
        """

        #add Surface module if requested            
        s = Surface()
        s.enable_contours = True
        s.actor.property.opacity = 0.5
        mayavi.add_module(s)

        #add Text module if requested
        t1 = Text()
        t1.text = self.title
        t1.actor.scaled_text = False
        t1.actor.text_property.font_size = 18
        mayavi.add_module(t1)
        t1.width = 1.0*t1.actor.mapper.get_width(t1.scene.renderer)/t1.scene.renderer.size[0]
        height = 1.0*t1.actor.mapper.get_height(t1.scene.renderer)/t1.scene.renderer.size[1]
        t1.x_position = 0.5 - t1.width/2
        t1.y_position = 1-height

        if(addTextIndexer):
            #add default Text module for indicating scene index
            self.sceneTextCount = Text()
            self.sceneTextCount.text = "0"
            mayavi.add_module(self.sceneTextCount)
            self.sceneTextCount.actor.scaled_text = False
            self.sceneTextCount.actor.text_property.font_size = 24
            self.sceneTextCount.x_position = .95
            self.sceneTextCount.y_position = .05

    def snapshot(self):
        #take snapshot
        if len(self.outputdir) > 0: #if a path was entered
            e = mayavi.engine
            s = mayavi.engine.current_scene
            s.scene.save_png("%s/output%05d.png" %(self.outputdir, self.imagecount))
            self.imagecount = self.imagecount + 1
            

if __name__=='__main__': 
    files=parse()
    MayaViShow(files,"TEST").main()
