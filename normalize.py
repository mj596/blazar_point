#! /usr/bin/env python

import numpy as np
import sys

DIR=sys.argv[1]
f=open(DIR+'/Lflare','r')
arr=[]
t=[]
for line in f:
   t.append(line.split()[0])
   arr.append(map(float,line.split()))


narr=np.array(arr)

print "%s\t%1.2e" % (DIR, narr.max())

f.close()
data=[]
f=open(DIR+'/Lflare','r')
for line in f:
    data.append(line.split())

kk=[]
for i in data:
    kk.append(narr.max())

ndata=np.array(data)
nkk=np.array(kk)

f.close()

f=open(DIR+'/Lflaremod','w')

c1=np.transpose(ndata)[0]
c2=np.transpose(ndata)[1]
c3=nkk

fc1=[]
fc2=[]
fc3=[]
fc4=[]

for ii in range(len(c1)):
    fc1.append(float(c1[ii]))
    fc2.append(float(c2[ii]))
    fc3.append(float(c3[ii]))

nfc1=np.array(fc1)
nfc2=np.array(fc2)
nfc3=np.array(fc3)

print nfc1[np.max(nfc2.argmax())]

for i in range(len(c1)):
    fc4.append(nfc1[np.max(nfc2.argmax())])

nfc4=np.array(fc4)

a=np.vstack((nfc1,nfc2,nfc3,nfc4))
ta=np.transpose(a)

np.savetxt(f,ta)

#DIR=sys.argv[1]
#f=open(DIR+'s/Lflare','r')
#sarr=[]
#for line in f:
#    sarr.append(map(float,line.split()))
#
#snarr=np.array(sarr)
#
#g=open(DIR+'g/Lflare','r')
#garr=[]
#for line in g:
#    garr.append(map(float,line.split()))
#
#gnarr=np.array(garr)
#
#
#print "id\tsyn\t\tgamma"
#print "%s\t%1.2e\t%1.2e" % (DIR, snarr.max(), gnarr.max())

