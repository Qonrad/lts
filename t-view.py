import sys
from numpy import *
from matplotlib import pyplot as plt

if len(sys.argv) < 2 :
	print "USAGE {} datafile".format(sys.argv[0])
	exit(1)
	
vdata = genfromtxt(sys.argv[1])
ncell = vdata.shape[1]

ax = plt.subplot(212)
spks = []
for i in xrange(1,ncell):
	spk = where(vdata[:,i] > 0.0)[0]
	spk = [ (vdata[spk[0],0],i) ] + [ (vdata[k,0],i) for p,k in zip(spk[:-1],spk[1:]) if p != k-1  ]
	plt.plot(vdata[:,0], vdata[:,i])
	spks += spk
plt.subplot(211, sharex=ax)
spks = array(spks)
plt.plot(spks[:,0], spks[:,1], "k.")
plt.show()
