from numpy import *
from matplotlib import pyplot as plt
from matplotlib.cm import get_cmap

import sys,os
if len(sys.argv) < 3:
	print "USAGE: python {} number_of_neurons datafile datafile ....".format(sys.argv[0])
	exit(1)

try:
	n = int(sys.argv[1])
except:
	print "couldn't get {} to an integer number".format(sys.argv[2])
	exit(1)
datas = []
for names in sys.argv[2:]:
	sys.stderr.write(" reading {} ...".format(names ))
	datas.append(genfromtxt(names))
	sys.stderr.write(" DONE\n")
	
	if (datas[-1].shape[1]-1)%n != 0:
		print "couldn't find n={} neurons in the file".format(n,sys.argv[1])
		exit(1)

nvars = (datas[-1].shape[1]-1)/n
Nplots = (nvars + 1)*100+10


cmap = get_cmap('jet')
f=plt.figure(1)
ax=plt.subplot(Nplots+1)
for ixd, data in enumerate(datas):
	c = cmap(float(ixd)/float(len(datas)))
	spks = []
	for i in xrange(n):
		print data.shape, i, i*nvars+1
		spk = where(data[:,i*nvars+1] > 0.0)[0]
		if len(spk) > 0:
			spk = [ (data[spk[0],0],i) ] + [ (data[k,0],i) for p,k in zip(spk[:-1],spk[1:]) if p != k-1  ]
			spks += spk
	spks = array(spks)
	plt.plot(spks[:,0], spks[:,1], ".",mfc=c,mec=c)

neuron = 0
curves = []
for ixd, data in enumerate(datas):
	c = cmap(float(ixd)/float(len(datas)))
	curves.append([])
	for nv in xrange(nvars):
		plt.subplot(Nplots+2+nv,sharex=ax)
		curves[ixd].append( plt.plot(data[:,0],data[:,neuron*nvars+nv+1], "-", c=c)[0] )

def kpress(event):
	global neuron,curves
	if event.key == "up":
		neuron += 1
		if neuron >= n : neuron =0
	elif event.key == "down":
		neuron -= 1
		if neuron < 0 : neuron = n-1
	for ixd, data in enumerate(datas):
		for nv,c in enumerate(curves[ixd]):
			c.set_ydata( data[:,neuron*nvars+nv+1] )
	f.canvas.draw()

f.canvas.mpl_connect('key_press_event',kpress)


plt.show()
