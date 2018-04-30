import os
import sys
import fnmatch
import shutil

begin = sys.argv[1]
end = sys.argv[2]

fd = open('duration.txt','r')
fw = open('sembmaxvalue.txt','r')

fdc=0
for i in fd:
    if fdc == 0:
        line1 = i
    fdc+=1    
fd.close()

c=0
for i in fw:
    line = str.split(i,' ')
    if fnmatch.fnmatch(line[1], begin):
        beginline = c
        begincontext = line
    if fnmatch.fnmatch(line[1], end):
        endline = c
        endcontext = line
    
    c+=1    
    
fw.seek(0)   

fnd = open('tmp.txt','w')
fnd.write(line1)
line2 ='# Shifted event onset at i = '+ begincontext[0] +' time ='+begincontext[1]+' semblance = '+begincontext[2]+'\n'
fnd.write(line2)

counter=0
for i in fw:
    if counter >= beginline and counter <=endline:
        fnd.write(i)
        line = str.split(i)
    counter+=1
    
line3='# Shiftet event stop at i = '+endcontext[0]+' time = '+endcontext[1]+' semblance = '+endcontext[2]+'\n'
fnd.write(line3)
duration = int(end)-int(begin)
ll ='# Event duration = '+str(duration)    
fnd.write(ll)
fnd.close()

shutil.copyfile('tmp.txt', 'duration.txt')
cmd='./summary-plot.sh test'
os.system(cmd)