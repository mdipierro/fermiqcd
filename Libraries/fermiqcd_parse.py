import re
data=open("fermiqcd.cpp","r").read()
regex = re.compile('\("(.*?)","(.*?)","?(.*?)"?\)')
d = {}
for item in regex.findall(data):
    d[item[0]]=d.get(item[0],[])+[(item[1],item[2])]
regex = re.compile('have\("(.*?)"\)')
for item in regex.findall(data):
    d[item]=d.get(item,[])
for key in sorted(d):
    print '    "%s\\n"' % key
    for a,b in d[key]:
        print '    "   %s=%s\\n"' % (a,b)
