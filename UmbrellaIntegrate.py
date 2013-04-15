#!/usr/bin/env python

usage="./UmbrellaIntegrate.py metadata.txt path.dat mean-force.dat pmf.dat"


### parse the args ###

from sys import argv
try:
  import numpy as np
except:
  print "You must have numpy to use this script! Try \">$ module load python\" on clusters"
  exit(0)

if len(argv)==5:
  metadata=argv[1]
  pathfile=argv[2]
  mfrcfile=argv[3]
  freefile=argv[4]
else:
  print usage
  exit(0)

#####################


### functions ###

def diff(a,b):
  return [b[i]-a[i] for i in range(len(a))]

def dist(a,b):
  v=diff(a,b)
  return np.sqrt(np.inner(v,v))

def zip(a,b):
  return [[a[i],b[i]] for i in range(len(a))]

def proj(a,b):
  k=np.dot(a,b)/np.inner(a,a)
  return [k*x for x in a]

#################

### get the data ###

files=[l.split()[0] for l in open(metadata,'r').readlines()]
ks=[float(l.split()[2]) for l in open(metadata,'r').readlines()]
path=[[float(x) for x in l.split()] for l in open('path.dat','r').readlines()]

####################

### compute the displacement from the window center ###
###             (projected onto the path)           ###

avgs=[]
i=0
for f in files:
  dat=[[float(x) for x in l.split()] for l in open(f,'r').readlines()[3:]]
  if i>0:
    bvec=diff(path[i],path[i-1])
  if i<len(files)-1:
    fvec=diff(path[i+1],path[i])
  # for i=0, i=len(path) default to the forward / backward projections
  disps=[]
  if i>0 and i<len(files):
    for pt in dat:
      disp=diff(pt,path[i])
      bprj=-dist(proj(bvec,disp),bvec)
      fprj=dist(proj(fvec,disp),[0]*len(fvec))
      if abs(bprj)>fprj:
        disp_prj=fprj
      else:
        disp_prj=bprj
      disps.append(disp_prj)
  elif i==0:
    for pt in dat:
      disp=diff(pt,path[i])
      disps.append(dist(proj(fvec,disp),[0]*len(fvec)))
  elif i==len(files)-1: 
    for pt in dat:
      disp=diff(pt,path[i])
      disps.append(-dist(proj(bvec,disp),bvec))
  avgs.append(-np.mean(disps)*ks[i])
  print "Average Force on pathpoint %d is %8.3f"%(i,avgs[i])
  i+=1

######################################################


### compute the pmf ###

# project the path onto a line
proj=[0]
pathd=0
for i in range(len(path)-1):
  pathd+=dist(path[i],path[i+1])
  proj.append(pathd)

# write a data file with the interpolated average force
xpts=[x*pathd/10000 for x in range(10000)]
mfrc=np.interp(xpts,proj,avgs)
mfrcplot=zip(xpts,mfrc)
forcefile=open(mfrcfile,'w')
for pt in mfrcplot:
  forcefile.write("%8.3f %8.3f\n"%tuple(pt))

# integrate to get the PMF
pmffile=open(freefile,'w')
for i in range(10000):
  pmffile.write("%8.3f %8.3f\n"%(xpts[i],np.trapz(mfrc[:i],xpts[:i])))

#######################

#avgf=open('string-force.dat','w')
#avgf.write('x\ty\t\z\tf\n')
#for i in range(len(avgs)):
#  avgf.write('%8.3f %8.3f %8.3f'%tuple(path[i]))
#  avgf.write(' %8.3f\n'%avgs[i])
