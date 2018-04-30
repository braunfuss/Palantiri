import os
import fnmatch
import subprocess
import sys

D=[]
fobj = open('duration.txt','r')
for line in fobj:
	if line[0] != '#':
		line = line.split()
		D.append(float(line[2]))
fobj.close()

minVal= min(D)
maxVal= max(D)

fobj = open('minmax.txt','w')
fobj.write(('%f,%f')%(minVal,maxVal))
fobj.close()

L=[]
for i in os.listdir(os.getcwd()):
	if fnmatch.fnmatch(i,'*_*.ASC'):
		i = i.split('_')
		L.append(i[0])
K = list(set(L))
print '{ndepth} Depths found'.format(ndepth=len(K))

if len(sys.argv) < 2:
    print 'process all depths\n'
    K= sorted(K)
else:
  print 'process only give depths\n'
  K=sys.argv[1].split(',')
print sorted(K),type(K)

for depth in sorted(K):
	print depth
	print "start"
        subprocess.call("./semblance-plot_multiple_depths.sh test "+depth, shell=True)
	print "end"

