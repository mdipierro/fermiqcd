import sys, os, struct, mmap, datetime, re, glob, optparse

X,Y,Z,T=1,2,3,0

def reunitarize(s):
    (a1re, a1im, a2re, a2im, a3re, a3im, b1re, b1im, b2re, b2im, b3re, b3im)=s

    c1re = a2re * b3re - a2im * b3im - a3re * b2re + a3im * b2im
    c1im = -(a2re * b3im + a2im * b3re - a3re * b2im - a3im * b2re)

    c2re = a3re * b1re - a3im * b1im - a1re * b3re + a1im * b3im
    c2im = -(a3re * b1im + a3im * b1re - a1re * b3im - a1im * b3re)

    c3re = a1re * b2re - a1im * b2im - a2re * b1re + a2im * b1im
    c3im = -(a1re * b2im + a1im * b2re - a2re * b1im - a2im * b1re)
    return (a1re, a1im, a2re, a2im, a3re, a3im,
            b1re, b1im, b2re, b2im, b3re, b3im,
            c1re, c1im, c2re, c2im, c3re, c3im)

class MilcField:
    site_order=[T,Z,Y,X]
    link_order=[X,Y,Z,T]
    def __init__(self,file):
        self.ifile=file
        self.lattice_size=[0,0,0,0]
        self.precision='f'
        self.precision_target='f'
        self.endianess='>'
        self.offset=None
    def get_header(self):
        ifile=self.ifile
        header=ifile.read(96)
        if not len(header)==96: raise SyntaxError, "file too small"
        self.endianess='<'
        milc_head=struct.unpack('<i4i64siii',header)            
        if not milc_head[0]==20103:
            self.endianess='>'
            milc_head=struct.unpack('>i4i64siii',header)
        if not milc_head[0]==20103:
            raise SyntaxError
        self.lattice_size=[milc_head[4],milc_head[1],milc_head[2],milc_head[3]]
        self.offset=96
        ifile.seek(self.offset)
        return None
    def detect(self):
        try: self.get_header()
        except: return false
        else: return True
    def convert(self,filename):
        self.get_header()
        ifile=self.ifile
        ofile=open(filename,'wb')
        oheader_format='<60s60s60sLi10iii'
        if self.precision=='f': site_size=288
        else: site_size=576
        link_format=self.precision*(9*2) # 2 for complex
        link_format_out='<'+self.precision_target*(9*2) # 2 for complex
        if self.precision=='f': link_size=9*8
        else: link_size=9*16
        nx=self.lattice_size[X]
        ny=self.lattice_size[Y]
        nz=self.lattice_size[Z]
        nt=self.lattice_size[T]
        n0=self.lattice_size[self.site_order[0]] #T
        n1=self.lattice_size[self.site_order[1]] #Z
        n2=self.lattice_size[self.site_order[2]] #Y
        n3=self.lattice_size[self.site_order[3]] #X
        print nt,nx,ny,nz
        oheader=struct.pack(oheader_format,'File Type: MDP FIELD',
                            filename,datetime.datetime.now().isoformat(),
                            1325884739,4,nt,nx,ny,nz,0,0,0,0,0,0,
                            site_size,nt*nx*ny*nz)
        ooffset=struct.calcsize(oheader_format)
        ofile.write(oheader)        
        offset=self.offset
        site=[0,0,0,0]        
        ifile.seek(offset)
        for p0 in range(n0):
            print 'timeslice',p0,'...'
            for p1 in range(n1):
                for p2 in range(n2):
                    for p3 in range(n3):
                        p=[0,0,0,0]
                        p[self.site_order[0]]=p0 #T
                        p[self.site_order[1]]=p1 #Z
                        p[self.site_order[2]]=p2 #Y
                        p[self.site_order[3]]=p3 #X
                        ofile.seek(ooffset+site_size*(p[Z]+nz*(p[Y]+ny*(p[X]+nx*p[T]))))
                        for mu in self.link_order:
                             site[mu]=struct.unpack(self.endianess+link_format,
                                                    ifile.read(link_size))
                        for mu in [T,X,Y,Z]:
                            ofile.write(struct.pack(link_format_out,*site[mu]))
                            for i,r in enumerate(site[mu]):
                                if abs(r)>1: raise SyntaxError, "Invalid Format"
        print 'done'

class QioField:
    site_order=[T,Z,Y,X]
    link_order=[X,Y,Z,T]
    def __init__(self,file):
        self.ifile=file
        self.lattice_size=[0,0,0,0]
        self.precision='f'
        self.precision_target='f'
        self.endianess='>'
        self.offset=None
    def get_header(self):
        ifile=self.ifile
        offset=0
        while True:
            header=ifile.read(144)
            if not len(header)==144: break
            lime_head=struct.unpack('!2i1q128s',header)            
            if not lime_head[0]==1164413355: raise SyntaxError
            data=lime_head[2]
            print 'Lime block:',data,lime_head[3]
            padding=(8 - (data % 8)) % 8
            if lime_head[3][:16]=='ildg-binary-data': self.offset=144+offset
            if lime_head[3][:11]=='ildg-format':
                 regex=re.compile('\<lx\>(?P<lx>\d+)\</lx\>\<ly\>(?P<ly>\d+)\</ly\>\<lz\>(?P<lz>\d+)\</lz\>\<lt\>(?P<lt>\d+)\</lt\>')
                 ell=regex.search(ifile.read(data))
                 self.lattice_size=[int(ell.group('lt')),int(ell.group('lx')),int(ell.group('ly')),int(ell.group('lz'))]
                 print 'Lattice size from header:',self.lattice_size
            offset+=144+data+padding
            ifile.seek(offset)
        return None
    def detect(self):
        try: self.get_header()
        except: return false
        else: return True
    def convert(self,filename):
        self.get_header()
        ifile=self.ifile
        ofile=open(filename,'wb')
        oheader_format='<60s60s60sLi10iii'
        if self.precision=='f': site_size=288
        else: site_size=576
        link_format=self.precision*(9*2) # 2 for complex
        link_format_out='<'+self.precision_target*(9*2) # 2 for complex
        if self.precision=='f': link_size=9*8
        else: link_size=9*16
        nx=self.lattice_size[X]
        ny=self.lattice_size[Y]
        nz=self.lattice_size[Z]
        nt=self.lattice_size[T]
        n0=self.lattice_size[self.site_order[0]] #T
        n1=self.lattice_size[self.site_order[1]] #Z
        n2=self.lattice_size[self.site_order[2]] #Y
        n3=self.lattice_size[self.site_order[3]] #X
        oheader=struct.pack(oheader_format,'File Type: MDP FIELD',
                            filename,datetime.datetime.now().isoformat(),
                            1325884739,4,nt,nx,ny,nz,0,0,0,0,0,0,
                            site_size,nt*nx*ny*nz)
        ooffset=struct.calcsize(oheader_format)
        ofile.write(oheader)        
        offset=self.offset
        site=[0,0,0,0]        
        ifile.seek(offset)
        for p0 in range(n0):
            print 'timeslice',p0,'...'
            for p1 in range(n1):
                for p2 in range(n2):
                    for p3 in range(n3):
                        p=[0,0,0,0]
                        p[self.site_order[0]]=p0 #T
                        p[self.site_order[1]]=p1 #Z
                        p[self.site_order[2]]=p2 #Y
                        p[self.site_order[3]]=p3 #X
                        ofile.seek(ooffset+site_size*(p[Z]+nz*(p[Y]+ny*(p[X]+nx*p[T]))))
                        for mu in self.link_order:
                             site[mu]=struct.unpack(self.endianess+link_format,
                                                    ifile.read(link_size))
                        for mu in [T,X,Y,Z]:
                            ofile.write(struct.pack(link_format_out,*site[mu]))
                            for i,r in enumerate(site[mu]):
                                if abs(r)>1: raise SyntaxError, "Invalid Format"
        print 'done'

class ILDGField:
    site_order=[T,Z,Y,X]
    link_order=[X,Y,Z,T]
    def __init__(self,file):
        self.ifile=file
        self.lattice_size=[0,0,0,0]
        self.precision='f'
        self.precision_target='f'
        self.endianess='>'
        self.offset=None
    def get_header(self):
        ifile=self.ifile
        offset=0
        while True:
            header=ifile.read(144)
            if not len(header)==144: break
            lime_head=struct.unpack('!2i1q128s',header)
            if not lime_head[0]==1164413355: raise SyntaxError
            data=lime_head[2]
            print 'Lime block:',data,lime_head[3]
            padding=(8 - (data % 8)) % 8
            if lime_head[3][:16]=='ildg-binary-data': self.offset=144+offset
            if lime_head[3][:11]=='ildg-format':
                 regex=re.compile('\<lx\>(?P<lx>\d+)\</lx\>\<ly\>(?P<ly>\d+)\</ly\>\<lz\>(?P<lz>\d+)\</lz\>\<lt\>(?P<lt>\d+)\</lt\>')
                 ell=regex.search(ifile.read(data))
                 self.lattice_size=[int(ell.group('lt')),int(ell.group('lx')),int(ell.group('ly')),int(ell.group('lz'))]
                 print 'Lattice size from header:',self.lattice_size
            if lime_head[3][:18]=='scidac-binary-data': self.offset=144+offset
            if lime_head[3][:23]=='scidac-private-file-xml':
                 regex=re.compile('\<dims\>\D*(?P<dims>\d+.*\d+)\D*\</dims\>')
                 ell=regex.search(ifile.read(data))
                 tt=re.split('\D+',ell.group('dims'))
                 self.lattice_size=map(int,[tt[-1]]+tt[0:-1])
                 print 'Lattice size from header:',self.lattice_size
            offset+=144+data+padding
            ifile.seek(offset)
        return None
    def detect(self):
        try: self.get_header()
        except: return false
        else: return True
    def convert(self,filename):
        self.get_header()
        ifile=self.ifile
        ofile=open(filename,'wb')
        oheader_format='<60s60s60sLi10iii'
        if self.precision=='f': site_size=288
        else: site_size=576
        link_format=self.precision*(9*2) # 2 for complex
        link_format_out='<'+self.precision_target*(9*2) # 2 for complex
        if self.precision=='f': link_size=9*8
        else: link_size=9*16
        nx=self.lattice_size[X]
        ny=self.lattice_size[Y]
        nz=self.lattice_size[Z]
        nt=self.lattice_size[T]
        n0=self.lattice_size[self.site_order[0]] #T
        n1=self.lattice_size[self.site_order[1]] #Z
        n2=self.lattice_size[self.site_order[2]] #Y
        n3=self.lattice_size[self.site_order[3]] #X
        oheader=struct.pack(oheader_format,'File Type: MDP FIELD',
                            filename,datetime.datetime.now().isoformat(),
                            1325884739,4,nt,nx,ny,nz,0,0,0,0,0,0,
                            site_size,nt*nx*ny*nz)
        ooffset=struct.calcsize(oheader_format)
        ofile.write(oheader)        
        offset=self.offset
        site=[0,0,0,0]        
        ifile.seek(offset)
        for p0 in range(n0):
            print 'timeslice',p0,'...'
            for p1 in range(n1):
                for p2 in range(n2):
                    for p3 in range(n3):
                        p=[0,0,0,0]
                        p[self.site_order[0]]=p0 #T
                        p[self.site_order[1]]=p1 #Z
                        p[self.site_order[2]]=p2 #Y
                        p[self.site_order[3]]=p3 #X
                        ofile.seek(ooffset+site_size*(p[Z]+nz*(p[Y]+ny*(p[X]+nx*p[T]))))
                        for mu in self.link_order:
                             site[mu]=struct.unpack(self.endianess+link_format,
                                                    ifile.read(link_size))
                        for mu in [T,X,Y,Z]:
                            ofile.write(struct.pack(link_format_out,*site[mu]))
                            for i,r in enumerate(site[mu]):
                                if abs(r)>1: raise SyntaxError, "Invalid Format"
        print 'done'



class Nersc3x3Field:
    site_order=[T,Z,Y,X]
    link_order=[X,Y,Z,T]
    def __init__(self,file):
        self.ifile=file
        self.lattice_size=[0,0,0,0]
        self.precision='f'
        self.precision_target='f'
        self.endianess='>'
        self.offset=None
    def get_header(self):
        key="END_HEADER\n"
        self.ifile.seek(0)
        header=self.ifile.read(10000)
        self.offset=header.find(key)+len(key)
        if self.offset<len(key): raise SyntaxError, "not in NERSC format"
        header=header[:self.offset]
        print header
        self.lattice_size=[
           int(re.compile('DIMENSION_4\s*=\s*(?P<t>\d+)').search(header).group('t')),
           int(re.compile('DIMENSION_1\s*=\s*(?P<x>\d+)').search(header).group('x')),
           int(re.compile('DIMENSION_2\s*=\s*(?P<y>\d+)').search(header).group('y')),
           int(re.compile('DIMENSION_3\s*=\s*(?P<z>\d+)').search(header).group('z'))]
        if re.compile('FLOATING_POINT\s*=\s*IEEE64BIG').search(header):
            self.precision='d'
            self.endianess='>'
        elif re.compile('FLOATING_POINT\s*=\s*IEEE64SMALL').search(header):
            self.precision='d'
            self.endianess='<'
        elif re.compile('FLOATING_POINT\s*=\s*IEEE32SMALL').search(header):
            self.precision='f'
            self.endianess='<'
        elif re.compile('FLOATING_POINT\s*=\s*IEEE32(BIG)?').search(header):
            self.precision='f'
            self.endianess='>'
        else: raise SyntaxError, "invalid endianess or precision"
        if not re.compile('DATATYPE = 4D_SU3_GAUGE_3x3').search(header):
            raise SyntaxError, "DATATYPE = 4D_SU3_GAUGE_3x3"
        self.ifile.seek(self.offset)
        return None
    def detect(self):
        try: self.get_header()
        except: return false
        else: return True
    def convert(self,filename):
        self.get_header()
        ifile=self.ifile
        ofile=open(filename,'wb')
        oheader_format='<60s60s60sIi10iii'
        if self.precision_target=='f': site_size=288
        else: site_size=576
        link_format='<'+self.precision_target*(9*2) # 2 for complex
        if self.precision_target=='f': link_size=9*8
        else: link_size=9*16
        link_format_in=self.precision*(9*2) # 2 for complex
        if self.precision=='f': link_size_in=9*8
        else: link_size_in=9*16
        nx=self.lattice_size[X]
        ny=self.lattice_size[Y]
        nz=self.lattice_size[Z]
        nt=self.lattice_size[T]
        n0=self.lattice_size[self.site_order[0]] #T
        n1=self.lattice_size[self.site_order[1]] #Z
        n2=self.lattice_size[self.site_order[2]] #Y
        n3=self.lattice_size[self.site_order[3]] #X
        oheader=struct.pack(oheader_format,'File Type: MDP FIELD',
                            filename,datetime.datetime.now().isoformat(),
                            1325884739,4,nt,nx,ny,nz,0,0,0,0,0,0,
                            site_size,nt*nx*ny*nz)
        ooffset=struct.calcsize(oheader_format)
        ofile.write(oheader)
        offset=self.offset
        site=[0,0,0,0]
        ifile.seek(offset)
        for p0 in range(n0):
            print 'timeslice',p0,'...'
            for p1 in range(n1):
                for p2 in range(n2):
                    for p3 in range(n3):
                        p=[0,0,0,0]
                        p[self.site_order[0]]=p0 #T
                        p[self.site_order[1]]=p1 #Z
                        p[self.site_order[2]]=p2 #Y
                        p[self.site_order[3]]=p3 #X
                        for mu in self.link_order:
                             site[mu]=struct.unpack(self.endianess+link_format_in,
                                                    ifile.read(link_size_in))
                        ofile.seek(ooffset+site_size*(p[Z]+nz*(p[Y]+ny*(p[X]+nx*p[T]))))
                        for mu in [T,X,Y,Z]:
                            ofile.write(struct.pack(link_format,*site[mu]))
                            for i,r in enumerate(site[mu]):
                                if abs(r)>1: raise SyntaxError, "Invalid Format"
                            if mu==T and p0+p1+p2+p3==0: print site[mu]

class Nersc3x2Field:
    site_order=[T,Z,Y,X]
    link_order=[X,Y,Z,T]
    def __init__(self,file):
        self.ifile=file
        self.lattice_size=[0,0,0,0]
        self.precision='f'
        self.precision_target='f'
        self.endianess='>'
        self.offset=None
    def get_header(self):
        key="END_HEADER\n"
        self.ifile.seek(0)
        header=self.ifile.read(10000)
        self.offset=header.find(key)+len(key)
        if self.offset<len(key): raise SyntaxError, "not in NERSC format"
        header=header[:self.offset]
        print header
        self.lattice_size=[
           int(re.compile('DIMENSION_4\s*=\s*(?P<t>\d+)').search(header).group('t')),
           int(re.compile('DIMENSION_1\s*=\s*(?P<x>\d+)').search(header).group('x')),
           int(re.compile('DIMENSION_2\s*=\s*(?P<y>\d+)').search(header).group('y')),
           int(re.compile('DIMENSION_3\s*=\s*(?P<z>\d+)').search(header).group('z'))]
        if re.compile('FLOATING_POINT\s*=\s*IEEE64BIG').search(header):
            self.precision='d'
            self.endianess='>'
        elif re.compile('FLOATING_POINT\s*=\s*IEEE64SMALL').search(header):
            self.precision='d'
            self.endianess='<'
        elif re.compile('FLOATING_POINT\s*=\s*IEEE32SMALL').search(header):
            self.precision='f'
            self.endianess='<'
        elif re.compile('FLOATING_POINT\s*=\s*IEEE32(BIG)?').search(header):
            self.precision='f'
            self.endianess='>'
        else:
            self.precision='f'
            self.endianess='>'
            #raise SyntaxError, "invalid endianess or precision"
        if not re.compile('DATATYPE = 4D_SU3_GAUGE').search(header):
            raise SyntaxError, "DATATYPE = 4D_SU3_GAUGE"
        self.ifile.seek(self.offset)
        return None
    def detect(self):
        try: self.get_header()
        except: return false
        else: return True
    def convert(self,filename):
        self.get_header()
        ifile=self.ifile
        ofile=open(filename,'wb')
        oheader_format='<60s60s60sIi10iii'
        if self.precision_target=='f': site_size=288
        else: site_size=576
        link_format='<'+self.precision_target*(9*2) # 2 for complex
        if self.precision_target=='f': link_size=9*8
        else: link_size=9*16
        link_format_in=self.precision*(6*2) # 2 for complex
        if self.precision=='f': link_size_in=6*8
        else: link_size_in=6*16
        nx=self.lattice_size[X]
        ny=self.lattice_size[Y]
        nz=self.lattice_size[Z]
        nt=self.lattice_size[T]
        n0=self.lattice_size[self.site_order[0]] #T
        n1=self.lattice_size[self.site_order[1]] #Z
        n2=self.lattice_size[self.site_order[2]] #Y
        n3=self.lattice_size[self.site_order[3]] #X
        oheader=struct.pack(oheader_format,'File Type: MDP FIELD',
                            filename,datetime.datetime.now().isoformat(),
                            1325884739,4,nt,nx,ny,nz,0,0,0,0,0,0,
                            site_size,nt*nx*ny*nz)
        ooffset=struct.calcsize(oheader_format)
        ofile.write(oheader)
        offset=self.offset
        site=[0,0,0,0]
        ifile.seek(offset)
        for p0 in range(n0):
            print 'timeslice',p0,'...'
            for p1 in range(n1):
                for p2 in range(n2):
                    for p3 in range(n3):
                        p=[0,0,0,0]
                        p[self.site_order[0]]=p0 #T
                        p[self.site_order[1]]=p1 #Z
                        p[self.site_order[2]]=p2 #Y
                        p[self.site_order[3]]=p3 #X
                        for mu in self.link_order:
                            site[mu]=struct.unpack(self.endianess+link_format_in,
                                                   ifile.read(link_size_in))
                            for i,r in enumerate(site[mu]):
                                if abs(r)>1: raise SyntaxError, "Invalid Format %f" % r
                        ofile.seek(ooffset+site_size*(p[Z]+nz*(p[Y]+ny*(p[X]+nx*p[T]))))
                        for mu in [T,X,Y,Z]:
                            site[mu]=reunitarize(site[mu])
                            ofile.write(struct.pack(link_format,*site[mu]))
                            if mu==T and p0+p1+p2+p3==0: print site[mu]

def universal_converter(path, formats):
    files=[f for f in glob.glob(path) if not f[-4:]=='.mdp']    
    if not files: raise RuntimeError, "no files to be converted"
    done=False
    for file in files:
        print 'trying to convert '+file
        for format in formats:
            #try:
                ofile=file+'.mdp'
                format(open(file,'rb')).convert(ofile)
                done=True
                break
            #except Exception, e:
            #    print e
        if not done:
            print 'ERROR... skipping!'
    if not done: raise RuntimeError, "failure to convert "+file

class ILDGPropField:
    site_order=[T,Z,Y,X]
    link_order=[X,Y,Z,T]
    def __init__(self,file):
        # should check for total file size
        self.ifile=file
        self.lattice_size=[0,0,0,0]
        self.precision='f'
        self.precision_target='f'
        self.endianess='>'
        self.offset=None
    def get_header(self):
        ifile=self.ifile
        offset=0
        while True:
            header=ifile.read(144)
            if not len(header)==144: break
            lime_head=struct.unpack('!2i1q128s',header)
            if not lime_head[0]==1164413355: raise SyntaxError
            data=lime_head[2]
            print 'Lime block:',data,lime_head[3]
            padding=(8 - (data % 8)) % 8
            if lime_head[3][:16]=='ildg-binary-data': self.offset=144+offset
            if lime_head[3][:11]=='ildg-format':
                 regex=re.compile('\<lx\>(?P<lx>\d+)\</lx\>\<ly\>(?P<ly>\d+)\</ly\>\<lz\>(?P<lz>\d+)\</lz\>\<lt\>(?P<lt>\d+)\</lt\>')
                 ell=regex.search(ifile.read(data))
                 self.lattice_size=[int(ell.group('lt')),int(ell.group('lx')),int(ell.group('ly')),int(ell.group('lz'))]
                 print 'Lattice size from header:',self.lattice_size
            if lime_head[3][:18]=='scidac-binary-data': self.offset=144+offset
            if lime_head[3][:23]=='scidac-private-file-xml':
                 regex=re.compile('\<dims\>\D*(?P<dims>\d+.*\d+)\D*\</dims\>')
                 ell=regex.search(ifile.read(data))
                 tt=re.split('\D+',ell.group('dims'))
                 self.lattice_size=map(int,[tt[-1]]+tt[0:-1])
                 print 'Lattice size from header:',self.lattice_size
            offset+=144+data+padding
            ifile.seek(offset)
        return None
    def detect(self):
        try: self.get_header()
        except: return false
        else: return True
    def convert(self,filename):
        self.get_header()
        ifile=self.ifile
        oheader_format='<60s60s60sLi10iii'
        if self.precision=='f': site_size=1152
        else: site_size=2304
        spinor_format=self.precision*(16*9*2) # 2 for complex
        spinor_format_out='<'+self.precision_target*(16*9*2) # 2 for complex
        if self.precision=='f': spinor_size=16*9*8
        else: spinor_size=12*16
        nx=self.lattice_size[X]
        ny=self.lattice_size[Y]
        nz=self.lattice_size[Z]
        nt=self.lattice_size[T]
        n0=self.lattice_size[self.site_order[0]] #T
        n1=self.lattice_size[self.site_order[1]] #Z
        n2=self.lattice_size[self.site_order[2]] #Y
        n3=self.lattice_size[self.site_order[3]] #X
        oheader=struct.pack(oheader_format,'File Type: MDP FIELD',
                            filename,datetime.datetime.now().isoformat(),
                            1325884739,4,nt,nx,ny,nz,0,0,0,0,0,0,
                            spinor_size,nt*nx*ny*nz)
        ooffset=struct.calcsize(oheader_format)
        offset=self.offset
        site=[0,0,0,0]        
        ifile.seek(offset)
        ofile=open(filename,'wb')
        ofile.write(oheader)        
        for p0 in range(n0):
            print 'prop timeslice',p0,'...'
            for p1 in range(n1):
                for p2 in range(n2):
                    for p3 in range(n3):
                        p=[0,0,0,0]
                        p[self.site_order[0]]=p0 #T
                        p[self.site_order[1]]=p1 #Z
                        p[self.site_order[2]]=p2 #Y
                        p[self.site_order[3]]=p3 #X
                        data=struct.unpack(self.endianess+spinor_format,
                                           ifile.read(spinor_size))
                        ofile.seek(ooffset+spinor_size*(p[0]+nt*(p[Z]+nz*(p[Y]+ny*(p[X])))))
                        ofile.write(struct.pack(spinor_format_out,*data))
        print 'done'


class ILDGPropFieldSplit:
    site_order=[T,Z,Y,X]
    link_order=[X,Y,Z,T]
    def __init__(self,file):
        # should check for total file size
        self.ifile=file
        self.lattice_size=[0,0,0,0]
        self.precision='f'
        self.precision_target='f'
        self.endianess='>'
        self.offset=None
    def get_header(self):
        ifile=self.ifile
        offset=0
        while True:
            header=ifile.read(144)
            if not len(header)==144: break
            lime_head=struct.unpack('!2i1q128s',header)
            if not lime_head[0]==1164413355: raise SyntaxError
            data=lime_head[2]
            print 'Lime block:',data,lime_head[3]
            padding=(8 - (data % 8)) % 8
            if lime_head[3][:16]=='ildg-binary-data': self.offset=144+offset
            if lime_head[3][:11]=='ildg-format':
                 regex=re.compile('\<lx\>(?P<lx>\d+)\</lx\>\<ly\>(?P<ly>\d+)\</ly\>\<lz\>(?P<lz>\d+)\</lz\>\<lt\>(?P<lt>\d+)\</lt\>')
                 ell=regex.search(ifile.read(data))
                 self.lattice_size=[int(ell.group('lt')),int(ell.group('lx')),int(ell.group('ly')),int(ell.group('lz'))]
                 print 'Lattice size from header:',self.lattice_size
            if lime_head[3][:18]=='scidac-binary-data': self.offset=144+offset
            if lime_head[3][:23]=='scidac-private-file-xml':
                 regex=re.compile('\<dims\>\D*(?P<dims>\d+.*\d+)\D*\</dims\>')
                 ell=regex.search(ifile.read(data))
                 tt=re.split('\D+',ell.group('dims'))
                 self.lattice_size=map(int,[tt[-1]]+tt[0:-1])
                 print 'Lattice size from header:',self.lattice_size
            offset+=144+data+padding
            ifile.seek(offset)
        return None
    def detect(self):
        try: self.get_header()
        except: return false
        else: return True
    def convert(self,filename):
        self.get_header()
        ifile=self.ifile
        oheader_format='<60s60s60sLi10iii'
        if self.precision=='f': site_size=1152
        else: site_size=2304
        spinor_format=self.precision*(16*9*2) # 2 for complex
        spinor_format_out='<'+self.precision_target*(16*9*2) # 2 for complex
        if self.precision=='f': spinor_size=16*9*8
        else: spinor_size=12*16
        nx=self.lattice_size[X]
        ny=self.lattice_size[Y]
        nz=self.lattice_size[Z]
        nt=self.lattice_size[T]
        n0=self.lattice_size[self.site_order[0]] #T
        n1=self.lattice_size[self.site_order[1]] #Z
        n2=self.lattice_size[self.site_order[2]] #Y
        n3=self.lattice_size[self.site_order[3]] #X
        oheader=struct.pack(oheader_format,'File Type: MDP FIELD',
                            filename,datetime.datetime.now().isoformat(),
                            1325884739,4,1,nx,ny,nz,0,0,0,0,0,0,
                            spinor_size,1*nx*ny*nz)
        ooffset=struct.calcsize(oheader_format)
        offset=self.offset
        site=[0,0,0,0]        
        ifile.seek(offset)
        for p0 in range(n0):
            ofile=open(filename[:-4]+'t%.4i.mdp' % p0,'wb')
            ofile.write(oheader)        
            print 'prop timeslice',p0,'...'
            for p1 in range(n1):
                for p2 in range(n2):
                    for p3 in range(n3):
                        p=[0,0,0,0]
                        p[self.site_order[0]]=p0 #T
                        p[self.site_order[1]]=p1 #Z
                        p[self.site_order[2]]=p2 #Y
                        p[self.site_order[3]]=p3 #X
                        data=struct.unpack(self.endianess+spinor_format,
                                           ifile.read(spinor_size))
                        ofile.seek(ooffset+spinor_size*(p[Z]+nz*(p[Y]+ny*(p[X]))))
                        ofile.write(struct.pack(spinor_format_out,*data))
        print 'done'



formats = (
    ('gauge-milc',MilcField),
    ('gauge-qio',MilcField),
    ('gauge-lime',ILDGField),
    ('gauge-ildg',ILDGField),
    ('gauge-nersc3x3',Nersc3x3Field),
    ('gauge-nersc3x2',Nersc3x2Field),
    ('prop-ildg',ILDGPropField),
    ('prop-ildg-split',ILDGPropFieldSplit),
    )

def main():
    usage = "usage: %prog [options] path"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-f", "--format", dest="format", default="ALL",
                      help="any of the supported formats (%s)" % ','.join(f[0] for f in formats))
    (options, args) = parser.parse_args()
    myformats = [f[1] for f in formats if options.format in ('ALL',f[0])]
    universal_converter(args[0],myformats)

if __name__=='__main__': main()
