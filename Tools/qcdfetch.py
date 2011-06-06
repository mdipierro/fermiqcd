#!/usr/bin/python
# -*- coding: iso-8859-1 -*-

LOGO ="""
 #######   ######  ########     ######## ######## ########  ######  ##     ## 
##     ## ##    ## ##     ##    ##       ##          ##    ##    ## ##     ## 
##     ## ##       ##     ##    ##       ##          ##    ##       ##     ## 
##     ## ##       ##     ##    ######   ######      ##    ##       ######### 
##  ## ## ##       ##     ##    ##       ##          ##    ##       ##     ## 
##    ##  ##    ## ##     ##    ##       ##          ##    ##    ## ##     ## 
 ##### ##  ######  ########     ##       ########    ##     ######  ##     ## 
Created by Massimo Di Pierro - License GPL2 - all-to-all convertion utility
"""

USAGE ="""Usage:

  $ qcdfetch [options] source

Source can be a link to a NERSC ensable, a local file or a glob pattern 
(for example *.milc, folder/*.ildg)

Example for help:
  $ qcdfetch -h
 
Example downloaed an ensable:
  $ qcdfetch source

Example fetch an ensable in fermiqcd (.mdp) format.
  $ qcdfetch -c mdp source

Example fetch an ensable in fermiqcd (.mdp) format an convert to float
  $ qcdfetch -c mdp -4 source

Example fetch an ensable in fermiqcd (.mdp) format an convert to double
  $ qcdfetch -c mdp -8 source

Example fetch an ensable in ildg (.ildg)
  $ qcdfetch -c ildg source

Example fetch an ensable in fermiqcd format and break it into timeslices
  $ qcdfetch -c split.mdp source
  (will make one file per timeslice)

Example fetch a propagator and convert in fermiqcd(.prop.mdp) format
  $ qcdfetch -c prop.mdp source

Example fetch a propagator and convert in fermiqcd format and break into timeslices
  $ qcdfetch -c split.prop.mdp source
  (will make one file per timeslice)

Notice: qcdfetch runs in the folder here files are located unless, in case of 
donwloads, a --destination folder is specified.

"""

##### imports #############################################################

import urllib
import hashlib
import cPickle
import os
import re
import sys
import time
import datetime
import optparse
import struct
import mmap
import glob
import cStringIO
import termios
import signal
import array
import fcntl
import logging
import traceback
import xml.dom.minidom as dom
import xml.parsers.expat as expat

##### global variables #############################################################

NOW = datetime.datetime.now()
MAXBYTES = 1000  # max number of bytes for buffered reading
PRECISION = {'f':32,'d':64}
(X,Y,Z,T) = (1,2,3,0) # the MDP index convetion, used intenrnally

def notify(*a):
    """
    just an alies for the print function, in case we need to move to Python 3.x
    """
    print ' '.join([str(x) for x in a])

##### class Lime #############################################################

class Lime(object):
    """ based on this: http://usqcd.jlab.org/usqcd-docs/c-lime/lime_1p2.pdf"""
    @staticmethod
    def xml_parser(data):
        def f(name,dxml = dom.parseString(data)):
            return dxml.getElementsByTagName(name)[0].childNodes[0].nodeValue        
        return f
    def __init__(self,filename,mode,version = 1):
        """
        >>> lime = Lime('filename','r' or 'w')
        """
        self.magic = 1164413355
        self.version = version
        self.filename = filename
        self.mode = mode
        self.file = open(filename,mode)
        self.records = [] # [(name,position,size)]
        if mode == 'r' or mode == 'rb':
            while True:
                header = self.file.read(144)
                if not header: break
                magic, null,size, name = struct.unpack('!iiq128s',header)
                if magic != 1164413355:
                    raise IOError, "not in lime format"
                name = name[:name.find('\0')] # clenup junk from file
                position = self.file.tell()
                self.records.append((name,position,size)) # in bytes
                padding = (8 - (size % 8)) % 8
                self.file.seek(size+padding,1)
        self.dump_info()

    def dump_info(self,filename=None):
        f = open(filename or (self.filename+'.fromlime.info'),'w')
        f.write("LIME records:\n")
        for a,b,c in self:
            for a,b,c in self:
                f.write('- %s [%sbytes]\n' % (a,c))
                if c<1000:
                    f.write('\n'+b.read(c)+'\n\n')

    def read(self,record):
        """
        reads a Lime record
        >>> lime = Lime('filename','r')
        >>> name, reader, size = lime.read(records = 0)
        """
        if not self.mode in ('r','rb'):
            raise RuntimeError, "not suported"
        (name,position,size) = self.records[record]
        self.file.seek(position)
        return (name, self.file, size)
    def __iter__(self):
        """
        >>> lime = Lime('filename','r')
        >>> for name, reader, size in lime:
        >>> print name, size, reader.read(size)
        """
        for record in range(len(self)):
            yield self.read(record)
    def write(self,name,reader,size = None,chunk = MAXBYTES):
        """
        write a Lime record
        >>> lime = Lime('filename','w')
        >>> lime.write('record name','data',size = 4)
        data can be a string or a file object
        """
        if not self.mode in ('w','wb'):
            raise RuntimeError, "not supported"
        if isinstance(reader,str):
            if size == None:
                size = len(reader)
            reader = cStringIO.StringIO(reader)
        # write record header
        position = self.file.tell()
        header = struct.pack('!iiq128s',self.magic,self.version,size,name)
        self.file.write(header)
        # read data from reader and write to file
        if hasattr(reader,'read'):
            for i in xrange(size / chunk):
                data = reader.read(chunk)
                if len(data) != chunk:
                    raise IOError
                self.file.write(data)
            chunk = size % chunk
            data = reader.read(chunk)
            if len(data) != chunk:
                raise IOError
            self.file.write(data)
        else:
            for data in reader:
                self.file.write(data)
        # add padding bytes
        padding = (8 - (size % 8)) % 8
        self.file.write('\0'*padding)
        self.records.append((name,size,position))
    def close(self):
        self.file.close()
    def __len__(self):
        """
        returns the number of lime records
        """
        return len(self.records)
    def keys(self):
        """
        returns the name of lime records
        """
        return [name for (name,position,size) in self.records]

def test_lime():
    notify('making a dummy LIME file and writing junk in it...')
    lime = Lime('test.zzz.0.lime','w')
    lime.write('record1','01234567')
    lime.write('record2','012345678')
    file = cStringIO.StringIO('0123456789') # memory file
    lime.write('record3',file,10) # write file content as record
    lime.close()

    notify('reading the file back...')
    lime = Lime('test.zzz.0.lime','r')
    notify('file contans %s records' % len(lime))
    notify('they have names: %s' % lime.keys())

    for name,reader,size in lime:
        notify('record name: %s\nrecord size: %s\nrecord data: %s' % \
                   (name, size, reader.read(size)))
    lime.close()


##### reunitarize #############################################################

def reunitarize(s):
    (a1re, a1im, a2re, a2im, a3re, a3im, b1re, b1im, b2re, b2im, b3re, b3im) = s
    c1re = a2re * b3re - a2im * b3im - a3re * b2re + a3im * b2im
    c1im = -(a2re * b3im + a2im * b3re - a3re * b2im - a3im * b2re)
    c2re = a3re * b1re - a3im * b1im - a1re * b3re + a1im * b3im
    c2im = -(a3re * b1im + a3im * b1re - a1re * b3im - a1im * b3re)
    c3re = a1re * b2re - a1im * b2im - a2re * b1re + a2im * b1im
    c3im = -(a1re * b2im + a1im * b2re - a2re * b1im - a2im * b1re)
    return (a1re, a1im, a2re, a2im, a3re, a3im,
            b1re, b1im, b2re, b2im, b3re, b3im,
            c1re, c1im, c2re, c2im, c3re, c3im)


def check_unitarity(items,tolerance = 1.0):
    """
    this does not quite checks for unitarity, only that numbers are in [-1,+1] range
    """
    errors = [x for x in items if x<-tolerance or x>tolerance]
    if errors:
        raise RuntimeError, "matrix is not unitary"

def reorder(data,order1,order2): # data are complex numbers
    """
    reorders a list as in the example:
    >>> assert ''.join(reorder('AABBCCDD',[X,Y,Z,T],[Z,Y,X,T])) == 'CCBBAADD'
    """
    k = len(data)      # 4*9*2
    m = len(order1)    # 4
    n = k/m            # 9*2
    items = [None]*k
    for i in range(k):
        items[n*order1[i/n]+i%n] = data[i]
    items = [items[n*order2[i/n]+i%n] for i in range(k)]
    return items

assert ''.join(reorder('AABBCCDD',[X,Y,Z,T],[Z,Y,X,T])) == 'CCBBAADD'

##### Field readers #############################################################

class QCDFormat(object):
    site_order = [T,Z,Y,X] ### always unused but for reference
    link_order = [X,Y,Z,T] ### this is the order of links at the site level
    is_gauge = True or False
    def unpack(self,data):
        """
        unpacks a string of bytes from file into a list of float/double numbers
        """
        if self.precision.lower() == 'f':
            n = len(data)/4
        elif self.precision.lower() == 'd':
            n = len(data)/8
        else:
            raise IOError, "incorrect input precision"
        items = struct.unpack(self.endianess+str(n)+self.precision,data)
        if self.is_gauge:
            items = reorder(items,self.link_order,(T,X,Y,Z))
            check_unitarity(items)
        return items
    def pack(self,items):
        """
        packs a list of float/double numbers into a string of bytes
        """
        if self.is_gauge:
            items = reorder(items,(T,X,Y,Z),self.link_order)
        n = len(items)
        return struct.pack(self.endianess+str(n)+self.precision,*items)
    def __init__(self,filename): 
        """set defaults"""
        pass
    def read_header(self): 
        """read file header or fails"""
        return ('f',8,4,4,4) # open file
    def read_data(self,t,x,y,z):
        """random access read"""
        return (0,0,0,0,'data')
    def write_header(self,precision,nt,nx,ny,nz):
        """write header for new file"""
        pass    
    def write_data(self,data):
        """write next site variables, in order"""
        pass
    def close(self):
        """closes the file"""
        self.file.close()    

class GaugeCold(QCDFormat):
    def __init__(self,nt = 8,nx = 4,ny = 4,nz = 4):
        self.precision = 'f'
        self.size = (nt,nx,ny,nz)
    def read_header(self):
        (nt,nx,ny,nz) = self.size
        return (self.precision,nt,nx,ny,nz)
    def read_data(self,t,x,y,z):
        (nt,nx,ny,nz) = self.size
        return [1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 1.0, 0.0]


class GaugeMDP(QCDFormat):
    site_order = [T,X,Y,Z]
    link_order = [T,X,Y,Z]
    def __init__(self,filename,dummyfilename = 'none'):
        self.filename = filename
        self.dummyfilename = dummyfilename
        self.header_format = '<60s60s60sLi10iii'
        self.endianess = '<'
        self.header_size = 60+60+60+14*4
        self.offset = None
        self.site_size = None
        self.base_size = 4*9*2
    def read_header(self):
        self.file = open(self.filename,'rb')
        header = self.file.read(self.header_size)
        items = struct.unpack(self.header_format,header)
        if items[3] != 1325884739:
            notify('warning, this does not appear a MDP file, but could be wrong')
        nt,nx,ny,nz = items[5:9]
        self.site_size = items[15]
        if self.site_size == self.base_size*4:
            self.precision = 'f'
        elif self.site_size == self.base_size*8:
            self.precision = 'd'
        else:
            raise IOError, "unable to determine input precision"
        self.offset = self.file.tell()
        self.size = (nt,nx,ny,nz)
        return (self.precision,nt,nx,ny,nz)
    def write_header(self,precision,nt,nx,ny,nz):
        self.file = open(self.filename,'wb')
        self.site_size = self.base_size*(4 if precision == 'f' else 8)
        data = struct.pack(self.header_format,'File Type: MDP FIELD',
                           self.dummyfilename,NOW.isoformat(),
                           1325884739,4,nt,nx,ny,nz,0,0,0,0,0,0,
                           self.site_size,nt*nx*ny*nz)
        self.file.write(data)
        self.size = (nt,nx,ny,nz)
        self.precision = precision
        self.offset = self.file.tell()
    def read_data(self,t,x,y,z):
        (nt,nx,ny,nz) = self.size
        i = self.offset + (z+nz*(y+ny*(x+nx*t)))*self.site_size
        self.file.seek(i)
        data = self.file.read(self.site_size)
        return self.unpack(data)
    def write_data(self,data,target_precision = None):
        if len(data) != self.base_size:
            raise RuntimeError, "invalid data size"
        return self.file.write(self.pack(data))
    def convert_from(self,other,target_precision = None):        
        (precision,nt,nx,ny,nz) = other.read_header()
        print '  (precision: %s, size: %ix%ix%ix%i)' % (precision,nt,nx,ny,nz)
        self.write_header(target_precision or precision,nt,nx,ny,nz)
        pbar = ProgressBar(widgets = default_widgets , maxval = self.size[0]).start()
        for t in xrange(self.size[T]):
            for x in xrange(nx):
                for y in xrange(ny):
                    for z in xrange(nz):
                        data = other.read_data(t,x,y,z)
                        self.write_data(data)
            pbar.update(t)
        pbar.finish()


class GaugeMDPSplit(GaugeMDP):
    def convert_from(self,other,target_precision = None):        
        (precision,nt,nx,ny,nz) = other.read_header()
        print '  (precision: %s, size: %ix%ix%ix%i)' % (precision,nt,nx,ny,nz)
        pbar = ProgressBar(widgets = default_widgets , maxval = nt).start()
        for t in xrange(nt):
            slice = GaugeMDP(self.filename.replace('split.prop.mdp',
                                                   't%.4i.prop.mdp' % t))
            slice.write_header(target_precision or precision,1,nx,ny,nz)
            for x in xrange(nx):
                for y in xrange(ny):
                    for z in xrange(nz):
                        data = other.read_data(t,x,y,z)
                        slice.write_data(data)
            slice.close()
            pbar.update(t)
        pbar.finish()


class PropagatorMDP(QCDFormat):
    site_order = [T,X,Y,Z]
    is_gauge = False
    def __init__(self,filename):
        self.filename = filename
        self.header_format = '<60s60s60sLi10iii'
        self.endianess = '<'
        self.header_size = 60+60+60+14*4
        self.offset = None
        self.site_size = None
        self.base_size = 16*9*2
    def read_header(self):
        self.file = open(self.filename,'rb')
        header = self.file.read(self.header_size)
        items = struct.unpack(self.header_format,header)
        if items[3] != '1325884739':
            notify('warning, this does not appear a MDP file, but could be wrong')
        nt,nx,ny,nz = items[5:9]
        self.site_size = items[15]
        if self.site_size == self.base_size*4:
            self.precision = 'f'
        elif self.site_size == self.base_size*8:
            self.precision = 'd'
        else:
            raise IOError, "file not in GaugeMDP format"
        self.offset = self.file.tell()
        self.size = (nt,nx,ny,nz)
        return (self.precision,nt,nx,ny,nz)
    def write_header(self,precision,nt,nx,ny,nz):
        self.file = open(self.filename,'wb')
        self.site_size = self.base_size*(4 if precision == 'f' else 8)
        data = struct.pack(self.header_format,'File Type: MDP FIELD',
                           self.filename,NOW.isoformat(),
                           1325884739,4,nt,nx,ny,nz,0,0,0,0,0,0,
                           self.site_size,nt*nx*ny*nz)
        self.file.write(data)
        self.size = (nt,nx,ny,nz)
        self.precision = precision
        self.offset = self.file.tell()
    def read_data(self,t,x,y,z):
        i = self.offset + (z+nz*(y+ny*(x+nx*t)))*self.site_size
        self.file.seek(i)
        data = self.file.read(self.site_size)
        return self.unpack(data)
    def write_data(self,data,target_precision = None):
        if len(data) != self.base_size:
            raise RuntimeError, "invalid data size"
        return self.file.write(self.pack(data))
    def convert_from(self,other,target_precision = None):        
        (precision,nt,nx,ny,nz) = other.read_header()
        print '  (precision: %s, size: %ix%ix%ix%i)' % (precision,nt,nx,ny,nz)
        self.write_header(target_precision or precision,nt,nx,ny,nz)
        pbar = ProgressBar(widgets = default_widgets , maxval = self.size[0]).start()
        for t in xrange(nt):
            for x in xrange(nx):
                for y in xrange(ny):
                    for z in xrange(nz):
                        data = other.read_data(t,x,y,z)
                        self.write_data(data)
            pbar.update(t)
        pbar.finish()


class PropagatorMDPSplit(QCDFormat):
    site_order = [T,X,Y,Z]
    is_gauge = False
    def __init__(self,filename):
        self.filename = filename
        self.header_format = '<60s60s60sLi10iii'
        self.endianess = '<'
        self.header_size = 60+60+60+14*4
        self.offset = None
        self.site_size = None
        self.base_size = 16*9*2
    def write_header(self,precision,nt,nx,ny,nz):
        self.file = open(self.filename,'wb')
        self.site_size = self.base_size*(4 if precision == 'f' else 8)
        data = struct.pack(self.header_format,'File Type: MDP FIELD',
                           self.filename,NOW.isoformat(),
                           1325884739,4,nt,nx,ny,nz,0,0,0,0,0,0,
                           self.site_size,nt*nx*ny*nz)
        self.file.write(data)
        self.size = (nt,nx,ny,nz)
        self.precision = precision
        self.offset = self.file.tell()
    def write_data(self,data,target_precision = None):
        if len(data) != self.base_size:
            raise RuntimeError, "invalid data size"
        return self.file.write(self.pack(data))
    def convert_from(self,other,target_precision = None):        
        (precision,nt,nx,ny,nz) = other.read_header()
        print '  (precision: %s, size: %ix%ix%ix%i)' % (precision,nt,nx,ny,nz)
        pbar = ProgressBar(widgets = default_widgets , maxval = nt).start()
        for t in xrange(nt):
            slice = PropagatorMDP(self.filename.replace('.split.mdp','.t%.4i.mdp' % t))
            slice.write_header(target_precision or precision,1,nx,ny,nz)
            for x in xrange(nx):
                for y in xrange(ny):
                    for z in xrange(nz):
                        data = other.read_data(t,x,y,z)
                        slice.write_data(data)
            slice.close()
            pbar.update(t)
        pbar.finish()


class GaugeILDG(QCDFormat):
    def __init__(self,filename,lfn = 'unkown'):
        self.filename = filename
        self.endianess = '>'
        self.lfn = lfn
        self.field = 'su3gauge'
        self.base_size = 4*9*2
    def read_header(self):
        self.lime = Lime(self.filename,'r')
        self.file = self.lime.file
        for name,stream,size in self.lime:
            if name == 'ildg-binary-data':
                self.offset = stream.tell()
        for name,stream,size in self.lime:
            if name == 'ildg-format':
                data = stream.read(size)
                while data.endswith('\0'): data = data[:-1] # bug in generating data
                dxml = Lime.xml_parser(data)
                field = dxml("field")
                if field != self.field:
                    raise IOError, 'not a lime GaugeILDG'
                precision = int(dxml("precision"))
                nt = int(dxml("lt"))
                nx = int(dxml("lx"))
                ny = int(dxml("ly"))
                nz = int(dxml("lz"))
                if precision == 32:
                    self.precision = 'f'
                    self.site_size = self.base_size*4
                elif precision == 64:
                    self.precision = 'd'
                    self.site_size = self.base_size*8
                else:
                    raise IOError, "unable to determine input precision"
                self.size = (nt,nx,ny,nz)
                return (self.precision,nt,nx,ny,nz)
        raise IOError, "file is not in lime format"
    def write_header(self,precision,nt,nx,ny,nz):
        self.precision = precision
        self.site_size = 4*2*9*(4 if precision == 'f' else 8)
        self.size = (nt,nx,ny,nz)
        self.lime = Lime(self.filename,'w')
        self.file = self.lime.file
        precision = 32 if precision == 'f' else 64
        d = dict(field = 'su3gauge',version = '1.0',
                 precision = precision,lx = nx,ly = ny,lz = nz,lt = nt)
        data = """<?xml version = "1.0" encoding = "UTF-8"?>
<ildgFormat>
<version>%(version)s</version>
<field>su3gauge</field>
<precision>%(precision)s</precision>
<lx>%(lx)s</lx><ly>%(ly)s</ly><lz>%(lz)s</lz><lt>%(lt)s</lt>
</ildgFormat>""" % d
        self.lime.write('ildg-format',data)        
    def read_data(self,t,x,y,z):
        (nt,nx,ny,nz) = self.size
        i = self.offset + (x+nx*(y+ny*(z+nz*t)))*self.site_size
        self.file.seek(i)
        data = self.file.read(self.site_size)
        return self.unpack(data)
    def write_data(self,data,target_precision = None):
        if len(data) != self.base_size:
            raise RuntimeError, "invalid data size"
        return self.file.write(self.pack(data))
    def convert_from(self,other,target_precision = None):        
        (precision,nt,nx,ny,nz) = other.read_header()
        print '  (precision: %s, size: %ix%ix%ix%i)' % (precision,nt,nx,ny,nz)
        self.write_header(target_precision or precision,nt,nx,ny,nz)
        pbar = ProgressBar(widgets = default_widgets , maxval = self.size[0]).start()
        def reader():
            for t in xrange(nt):
                for z in xrange(nz):            
                    for y in xrange(ny):
                        for x in xrange(nx):
                            data = other.read_data(t,x,y,z)
                            yield self.pack(data)
                pbar.update(t)
        self.lime.write('ildg-binary-data',reader(),nt*nx*ny*nz*self.site_size)
        self.lime.write('ildg-data-LFN',self.lfn)
        self.lime.close()
        pbar.finish()


class GaugeSCIDAC(QCDFormat):
    def __init__(self,filename):
        self.filename = filename
        self.base_size = 4*9*2
        self.endianess = '>'
    def read_header(self):
        self.lime = Lime(self.filename,'r')
        self.file = self.lime.file
        for name,stream,size in self.lime:
            if name == 'scidac-binary-data':
                self.offset = stream.tell()
        self.size = self.precision = None
        for name,stream,size in self.lime:
            if name == 'scidac-private-file-xml':
                data = stream.read(size)
                while data.endswith('\0'): data = data[:-1] # bug in generating data
                dxml = Lime.xml_parser(data)
                dims = dxml("dims").strip().split()
                nt = int(dims[3])
                nx = int(dims[0])
                ny = int(dims[1])
                nz = int(dims[2])
                self.size = (nt,nx,ny,nz)
        for name,stream,size in self.lime:
            if name == 'scidac-private-record-xml':
                data = stream.read(size)
                while data.endswith('\0'): data = data[:-1] # bug in generating data
                dxml = Lime.xml_parser(data)
                precision = dxml("precision").lower()
                if precision == 'f':
                    self.precision = 'f'
                    self.site_size = self.base_size*4
                elif precision == 'd':
                    self.precision = 'd'
                    self.site_size = self.base_size*8
                else:
                    raise IOError, "unable to determine input precision"
        if self.size and self.precision:
            (nt,nx,ny,nz) = self.size
            return (self.precision,nt,nx,ny,nz)
        raise IOError, "file is not in lime format"
    def read_data(self,t,x,y,z):
        (nt,nx,ny,nz) = self.size
        i = self.offset + (x+nx*(y+ny*(z+nz*t)))*self.site_size
        self.file.seek(i)
        data = self.file.read(self.site_size)
        return self.unpack(data)


class PropagatorSCIDAC(QCDFormat):
    is_gauge = False
    def __init__(self,filename):
        self.filename = filename
        self.base_size = 16*9*2
        self.endianess = '>'


class GaugeMILC(QCDFormat):
    def __init__(self,filename):
        self.filename = filename
        self.header_format = '<i4i64siii' # may change
        self.endianess = '<' # may change
        self.header_size = 96 
        self.offset = None
        self.site_size = None
        self.base_size = 4*9*2
    def read_header(self):
        self.file = open(self.filename,'rb')
        header = self.file.read(self.header_size)
        for self.header_format in ('<i4i64siii','>i4i64siii'):
            self.endianess = self.header_format[0]
            items = struct.unpack(self.header_format,header)
            if items[0] == 20103:
                nt,nx,ny,nz = [items[4],items[1],items[2],items[3]]
                self.site_size = (os.path.getsize(self.filename)-96)/nt/nx/ny/nz
                self.size = (nt,nx,ny,nz)
                if self.site_size == self.base_size*4:
                    self.precision = 'f'
                elif self.site_size == self.base_size*8:
                    self.precision = 'd'
                else:
                    raise IOError, "file not in GaugeMILC fomat"
                self.offset = self.file.tell()
                self.size = (nt,nx,ny,nz)                
                return (self.precision,nt,nx,ny,nz)            
        raise IOError, "file not in GaugeMDP format"
    def read_data(self,t,x,y,z):
        (nt,nx,ny,nz) = self.size
        i = self.offset + (x+nx*(y+ny*(z+nz*t)))*self.site_size
        self.file.seek(i)
        data = self.file.read(self.site_size)
        return self.unpack(data)


class GaugeNERSC(QCDFormat):
    def __init__(self,filename):
        self.filename = filename
        self.offset = None
        self.site_size = None
        self.base_size = 4*9*2
        self.endianess = '>'
    def read_header(self):
        self.file = open(self.filename,'rb')
        header = self.file.read(100000)
        self.offset = header.find('END_HEADER\n')+11
        if self.offset<20:
            raise IOError, 'not in nersc format'
        lines = header[:self.offset-1].split('\n')[1:-2]
        info = dict([[x.strip() for x in item.split(' = ',1)] for item in lines])
        nx = int(info['DIMENSION_1'])
        ny = int(info['DIMENSION_2'])
        nz = int(info['DIMENSION_3'])
        nt = int(info['DIMENSION_4'])
        if info['FLOATING_POINT'].endswith('SMALL'):
            self.endianess = '<'
        else: # assume default big endinan
            self.endianess = '>'
        if info['DATATYPE'] == '4D_SU3_GAUGE_3x3':
            self.reunitarize = False
        elif info['DATATYPE'] == '4D_SU3_GAUGE':
            self.reunitarize = True
            self.base_size = 4*6*2
        else:
            raise IOError, "not in a known nersc format"
        if info['FLOATING_POINT'].startswith('IEEE32'):
            self.precision = 'f'
            self.site_size = self.base_size*4
        elif info['FLOATING_POINT'].startswith('IEEE64'):
            self.precision = 'd'
            self.site_size = self.base_size*8
        else:
            raise IOError, "unable to determine input precision"
        self.size = (nt,nx,ny,nz)        
        return (self.precision,nt,nx,ny,nz)        
    def read_data(self,t,x,y,z):
        (nt,nx,ny,nz) = self.size
        i = self.offset + (x+nx*(y+ny*(z+nz*t)))*self.site_size
        self.file.seek(i)
        data = self.file.read(self.site_size)
        items = self.unpack(data)
        if self.reunitarize:
            new_items = []
            for i in range(4):
                new_items += reunitarize(items[i*12:(i+1)*12])
            items = new_items        
        return items

OPTIONS = {
    'mdp':(GaugeMDP,GaugeMDP,GaugeMILC,GaugeNERSC,GaugeILDG,GaugeSCIDAC),
    'ildg':(GaugeILDG,GaugeILDG,GaugeMILC,GaugeNERSC,GaugeMDP,GaugeSCIDAC),
    'prop.mdp':(PropagatorMDP,PropagatorMDP,PropagatorSCIDAC),
    'prop.ildg':(PropagatorSCIDAC,PropagatorSCIDAC,PropagatorMDP),
    'split.mdp':(GaugeMDPSplit,GaugeMDP,GaugeMILC,GaugeNERSC,GaugeILDG,GaugeSCIDAC),
    'split.prop.mdp':(PropagatorMDPSplit,PropagatorMDP,PropagatorSCIDAC),
    }

def universal_converter(path,target,precision):
    filenames = [f for f in glob.glob(path) if not f[-4:] == '.'+target]
    if not filenames: raise RuntimeError, "no files to be converted"
    processed = set()
    messages = []
    option = OPTIONS[target]    
    for filename in filenames:        
        for formatter in option[1:]:
            messages.append('trying to convert %s (%s)' %(filename,formatter))
            try:
                ofilename = filename+'.'+target
                option[0](ofilename).convert_from(formatter(filename),precision)
                processed.add(filename)
                break
            except Exception, e:
                messages.append('unable to convert:\n' + traceback.format_exc())
        if not filename in processed:
            notify('\n'.join(messages))
            sys.exit(1)

##### BEGIN PROGRESSBAR ######
# progressbar  - Text progressbar library for python.
# Copyright (c) 2005 Nilton Volpato
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA


"""Text progressbar library for python.

This library provides a text mode progressbar. This is tipically used
to display the progress of a long running operation, providing a
visual clue that processing is underway.

The ProgressBar class manages the progress, and the format of the line
is given by a number of widgets. A widget is an object that may
display diferently depending on the state of the progress. There are
three types of widget:
- a string, which always shows itself;
- a ProgressBarWidget, which may return a diferent value every time
it's update method is called; and
- a ProgressBarWidgetHFill, which is like ProgressBarWidget, except it
expands to fill the remaining width of the line.

The progressbar module is very easy to use, yet very powerful. And
automatically supports features like auto-resizing when available.
"""

__author__ = "Nilton Volpato"
__author_email__ = "first-name dot last-name @ gmail.com"
__date__ = "2006-05-07"
__version__ = "2.2"

# Changelog
#
# 2006-05-07: v2.2 fixed bug in windows
# 2005-12-04: v2.1 autodetect terminal width, added start method
# 2005-12-04: v2.0 everything is now a widget (wow!)
# 2005-12-03: v1.0 rewrite using widgets
# 2005-06-02: v0.5 rewrite
# 2004-??-??: v0.1 first version


class ProgressBarWidget(object):
    """This is an element of ProgressBar formatting.

    The ProgressBar object will call it's update value when an update
    is needed. It's size may change between call, but the results will
    not be good if the size changes drastically and repeatedly.
    """
    def update(self, pbar):
        """Returns the string representing the widget.

        The parameter pbar is a reference to the calling ProgressBar,
        where one can access attributes of the class for knowing how
        the update must be made.

        At least this function must be overriden."""
        pass


class ProgressBarWidgetHFill(object):
    """This is a variable width element of ProgressBar formatting.

    The ProgressBar object will call it's update value, informing the
    width this object must the made. This is like TeX \\hfill, it will
    expand to fill the line. You can use more than one in the same
    line, and they will all have the same width, and together will
    fill the line.
    """
    def update(self, pbar, width):
        """Returns the string representing the widget.

        The parameter pbar is a reference to the calling ProgressBar,
        where one can access attributes of the class for knowing how
        the update must be made. The parameter width is the total
        horizontal width the widget must have.

        At least this function must be overriden."""
        pass


class ETA(ProgressBarWidget):
    "Widget for the Estimated Time of Arrival"
    def format_time(self, seconds):
        return time.strftime('%H:%M:%S', time.gmtime(seconds))
    def update(self, pbar):
        if pbar.currval == 0:
            return 'ETA:  --:--:--'
        elif pbar.finished:
            return 'Time: %s' % self.format_time(pbar.seconds_elapsed)
        else:
            elapsed = pbar.seconds_elapsed
            eta = elapsed * pbar.maxval / pbar.currval - elapsed
            return 'ETA:  %s' % self.format_time(eta)


class FileTransferSpeed(ProgressBarWidget):
    "Widget for showing the transfer speed (useful for file transfers)."
    def __init__(self):
        self.fmt = '%6.2f %s'
        self.units = ['B','K','M','G','T','P']
    def update(self, pbar):
        if pbar.seconds_elapsed < 2e-6:# == 0:
            bps = 0.0
        else:
            bps = float(pbar.currval) / pbar.seconds_elapsed
        spd = bps
        for u in self.units:
            if spd < 1000:
                break
            spd /= 1000
        return self.fmt % (spd, u+'/s')


class RotatingMarker(ProgressBarWidget):
    "A rotating marker for filling the bar of progress."
    def __init__(self, markers = '|/-\\'):
        self.markers = markers
        self.curmark = -1
    def update(self, pbar):
        if pbar.finished:
            return self.markers[0]
        self.curmark = (self.curmark + 1)%len(self.markers)
        return self.markers[self.curmark]

class Percentage(ProgressBarWidget):
    "Just the percentage done."
    def update(self, pbar):
        return '%3d%%' % pbar.percentage()


class Bar(ProgressBarWidgetHFill):
    "The bar of progress. It will strech to fill the line."
    def __init__(self, marker = '#', left = '|', right = '|'):
        self.marker = marker
        self.left = left
        self.right = right
    def _format_marker(self, pbar):
        if isinstance(self.marker, (str, unicode)):
            return self.marker
        else:
            return self.marker.update(pbar)
    def update(self, pbar, width):
        percent = pbar.percentage()
        cwidth = width - len(self.left) - len(self.right)
        marked_width = int(percent * cwidth / 100)
        m = self._format_marker(pbar)
        bar = (self.left + (m*marked_width).ljust(cwidth) + self.right)
        return bar


class ReverseBar(Bar):
    "The reverse bar of progress, or bar of regress. :)"
    def update(self, pbar, width):
        percent = pbar.percentage()
        cwidth = width - len(self.left) - len(self.right)
        marked_width = int(percent * cwidth / 100)
        m = self._format_marker(pbar)
        bar = (self.left + (m*marked_width).rjust(cwidth) + self.right)
        return bar

default_widgets = [Percentage(), ' ', Bar()]


class ProgressBar(object):
    """This is the ProgressBar class, it updates and prints the bar.

    The term_width parameter may be an integer. Or None, in which case
    it will try to guess it, if it fails it will default to 80 columns.

    The simple use is like this:
    >>> pbar = ProgressBar().start()
    >>> for i in xrange(100):
    ...    # do something
    ...    pbar.update(i+1)
    ...
    >>> pbar.finish()

    But anything you want to do is possible (well, almost anything).
    You can supply different widgets of any type in any order. And you
    can even write your own widgets! There are many widgets already
    shipped and you should experiment with them.

    When implementing a widget update method you may access any
    attribute or function of the ProgressBar object calling the
    widget's update method. The most important attributes you would
    like to access are:
    - currval: current value of the progress, 0 <= currval <= maxval
    - maxval: maximum (and final) value of the progress
    - finished: True if the bar is have finished (reached 100%), False o/w
    - start_time: first time update() method of ProgressBar was called
    - seconds_elapsed: seconds elapsed since start_time
    - percentage(): percentage of the progress (this is a method)
    """
    def __init__(self, maxval = 100, widgets = default_widgets, term_width = None,
                 fd = sys.stderr):
        assert maxval > 0
        self.maxval = maxval
        self.widgets = widgets
        self.fd = fd
        self.signal_set = False
        if term_width is None:
            try:
                self.handle_resize(None,None)
                signal.signal(signal.SIGWINCH, self.handle_resize)
                self.signal_set = True
            except:
                self.term_width = 79
        else:
            self.term_width = term_width

        self.currval = 0
        self.finished = False
        self.prev_percentage = -1
        self.start_time = None
        self.seconds_elapsed = 0

    def handle_resize(self, signum, frame):
        h,w = array.array('h', fcntl.ioctl(self.fd,termios.TIOCGWINSZ,'\0'*8))[:2]
        self.term_width = w

    def percentage(self):
        "Returns the percentage of the progress."
        return self.currval*100.0 / self.maxval

    def _format_widgets(self):
        r = []
        hfill_inds = []
        num_hfill = 0
        currwidth = 0
        for i, w in enumerate(self.widgets):
            if isinstance(w, ProgressBarWidgetHFill):
                r.append(w)
                hfill_inds.append(i)
                num_hfill += 1
            elif isinstance(w, (str, unicode)):
                r.append(w)
                currwidth += len(w)
            else:
                weval = w.update(self)
                currwidth += len(weval)
                r.append(weval)
        for iw in hfill_inds:
            r[iw] = r[iw].update(self, (self.term_width-currwidth)/num_hfill)
        return r

    def _format_line(self):
        return ''.join(self._format_widgets()).ljust(self.term_width)

    def _need_update(self):
        return int(self.percentage()) != int(self.prev_percentage)

    def update(self, value):
        "Updates the progress bar to a new value."
        assert 0 <= value <= self.maxval
        self.currval = value
        if not self._need_update() or self.finished:
            return
        if not self.start_time:
            self.start_time = time.time()
        self.seconds_elapsed = time.time() - self.start_time
        self.prev_percentage = self.percentage()
        if value != self.maxval:
            self.fd.write(self._format_line() + '\r')
        else:
            self.finished = True
            self.fd.write(self._format_line() + '\n')

    def start(self):
        """Start measuring time, and prints the bar at 0%.

        It returns self so you can use it like this:
        >>> pbar = ProgressBar().start()
        >>> for i in xrange(100):
        ...    # do something
        ...    pbar.update(i+1)
        ...
        >>> pbar.finish()
        """
        self.update(0)
        return self

    def finish(self):
        """Used to tell the progress is finished."""
        self.update(self.maxval)
        if self.signal_set:
            signal.signal(signal.SIGWINCH, signal.SIG_DFL)

###### END PROGRESSBAR #########

def md5_for_large_file(filename, block_size = 2**20):
    if not os.path.exists(filename):
        return None
    f = open(filename,'rb')
    md5 = hashlib.md5()
    while True:
        data = f.read(block_size)
        if not data:
            break
        md5.update(data)
    return md5.hexdigest()

def get_list(url):
    try:
        json = urllib.urlopen(url).read()
        files = eval(json)['files']
        return files
    except:
        return None

def download(files,target_folder,options):
    notify('total files to download: %s' % len(files))
    for k,f in enumerate(files):
        path = f['filename']
        basename = os.path.basename(path)
        target_name = os.path.join(target_folder,basename)
        if not os.path.exists(target_name) or \
                os.path.getsize(target_name) != f['size']:
            input = None
            widgets = [basename, Percentage(), ' ', Bar(marker = '>'),' ',
                       ETA(), ' ', FileTransferSpeed()]
            while not input:
                notify('downloading %s' % basename)
                try:
                    input = urllib.urlopen(f['link'])
                except IOError:
                    notify('failure to download %s' % f['link'])
                    sys.exit(1)
                length = int(input.info().get('Content-Length',f['size']))
                if not input:
                    notify('unable to retrieve %s retrying in 5 minutes' % basename)
                    time.sleep(5*60)
            if not options.quiet:
                pbar = ProgressBar(widgets = widgets, maxval = length).start()
            output = open(target_name,'wb')
            i = 0
            while True:
                data = input.read(MAXBYTES)
                if not data: break
                output.write(data)
                i+= len(data)
                if not options.quiet: pbar.update(i)
            input.close()
            output.close()
            if not options.quiet: pbar.finish()
            notify('completed downloads: %s/%s' % (k,len(files)))
        else:
            notify('skipping file %s (already present)' % basename)

def test_conversions():
    try:
        passed = False        
        test_lime()
        GaugeMDP('test.zzz.1.mdp').convert_from(GaugeCold(4,4,4,4))
        GaugeILDG('test.zzz.1.ildg').convert_from(GaugeMDP('test.zzz.1.mdp'))
        GaugeMDP('test.zzz.2.mdp').convert_from(GaugeILDG('test.zzz.1.ildg'))
        assert open('test.zzz.1.mdp','rb').read() == open('test.zzz.2.mdp','rb').read()
        GaugeMDP('test.zzz.3.mdp').convert_from(GaugeMDP('test.zzz.2.mdp'))
        assert open('test.zzz.1.mdp','rb').read() == open('test.zzz.3.mdp','rb').read()
        GaugeILDG('test.zzz.2.ildg').convert_from(GaugeILDG('test.zzz.1.ildg'))
        assert open('test.zzz.1.ildg','rb').read() == open('test.zzz.2.ildg','rb').read()
        #GaugeMDP('test.zzz.4.mdp').convert_from(GaugeNERSC('demo.nersc'))
        #GaugeILDG('test.zzz.4.ildg').convert_from(GaugeNERSC('demo.nersc'))
        #GaugeMDP('test.zzz.5.mdp').convert_from(GaugeILDG('test.zzz.4.ildg'))
        #assert open('test.zzz.4.mdp','rb').read() == open('test.zzz.5.mdp','rb').read()
        passed = True
    except Exception, e:
        notify('tests failed:\n%s' % traceback.format_exc())
    finally:
        notify('cleaning up tmp files')
        os.system('rm test.zzz.?.*')
    if not passed:
        sys.exit(1)


class DummyProgressBar(object):
    """used to avoid the ProgressBar on -n, --noprogressbar"""
    def __init__(self,*a,**b):
        self.maxval=b['maxval']
    def update(self,i):
        notify('...%s/%s' % (i,self.maxval))
    def start(self):
        return self
    def finish(self):
        notify('...done!')

def main():
    """this is the main program"""
    print LOGO
    parser = optparse.OptionParser(usage = USAGE)
    parser.add_option("-q", "--quiet",dest = 'quiet',action = 'store_true',
                      default = False,help = 'no progress bars')
    parser.add_option("-d", "--destination",dest = 'destination',default = None,
                      help = "destination folder")
    parser.add_option("-c", "--cpnvert",dest = 'convert',default = False,
                      help = "converts a field to format (mdp,ildg)")
    parser.add_option("-4", "--float",dest = 'float_precision',default = False,
                      action = 'store_true',
                      help = "converts to float precision")
    parser.add_option("-8", "--double",dest = 'double_precision',default = False,
                      action = 'store_true',
                      help = "converts to double precision")
    parser.add_option("-t", "--tests",dest = 'tests',default = False,
                      action = 'store_true',
                      help = "runs some tests")
    parser.add_option("-n", "--noprogressbar",dest = 'noprogressbar',default = False,
                      action = 'store_true',
                      help = "disable progress bar")
    (options, args) = parser.parse_args()

    ### disable progress bar if necessary
    if options.noprogressbar:
        global ProgressBar
        ProgressBar = DummyProgressBar

    ### run tests if asked
    if options.tests:
        test_conversions()
        sys.exit(0)

    ### if no argument source passed, print usage
    try:
        options.source = args[0]
    except IndexError:
        print USAGE
        sys.exit(1)

    ### download data (http, https, ftp, sftp) or not
    if options.source.startswith('http://') or options.source.startswith('https://'):
        files = get_list(options.source)
        if files == None:
            notify('unable to connect')
            sys.exit(0)
        else:
            regex = re.compile('pattern\ = (?P<pattern>[^\&]*)')
            pattern = regex.search(options.source).group('pattern')
            target_folder = options.destination or urllib.unquote(pattern)
            notify('target folder:',target_folder)
            if not os.path.exists(target_folder):
                os.mkdir(target_folder)
            download(files,target_folder,options)
        conversion_path = os.path.join(target_folder,pattern.replace('nnnnn','*'))
    elif options.source.startswith('ftp://') or options.source.startswith('sftp://'):
        username = raw_input('username:')
        password = raw_input('password:')
        if not os.path.exists(target_folder):
            os.mkdir(target_folder)
        ftp_download(options.source,target_folder,username,password)
    else:
        conversion_path = options.source

    ### if options.source = 'gauge.cold.TxXxYxZ' make it
    if options.source.startswith('gauge.cold') and not os.path.exists(options.source):
        size = [int(x) for x in options.source[11:].split('x')]
        GaugeMDP(options.source).convert_from(GaugeCold(*size))

    ### if conversion required use the universal converter
    if options.convert:
        notify('converting: %s -> %s.%s' % \
                   (conversion_path, conversion_path, options.convert))
        precision = 'f' if options.float_precision else \
            'd' if options.double_precision else None
        universal_converter(conversion_path,options.convert,precision)


if __name__ == '__main__': main()
