import sys
import re

def cleanup_cpp(data):
    data = re.compile('(?<=[\w\]\)])(?P<op>[\=\+\-\*])').sub(' \g<op>',data)
    data = re.compile('(?P<op>[\=\+\-\*])(?P<txt>[\w\-\+\.\[\{\"])').sub('\g<op> \g<txt>',data)
    data = re.compile(' \+ \+').sub('++',data)
    data = re.compile(' \- \-').sub('++',data)
    data = re.compile('(?P<op>[,;])(?P<txt>[^\s])').sub('\g<op> \g<txt>',data)
    data = re.compile('this ->').sub('this->',data)
    data = re.compile('\s*\{',re.DOTALL).sub(' {',data)
    data = re.compile('[\n]{2,10}',re.DOTALL).sub('\n\n',data)
    return data

data=open(sys.argv[1],'rb').read()
print cleanup_cpp(data)
