from numpy import *
from matplotlib import pyplot as plt
import sys,os
if len(sys.argv) < 3:
	print "USANGE: python {} datafile number_of_neurons".format(sys.argv[0])
	exit(1)

try:
	n = int(sys.argv[2])
except:
	print "couldn't get {} to an integer number".format(sys.argv[2])
	exit(1)

sys.stderr.write(" reading {} ...".format(sys.argv[1]) )
data = genfromtxt(sys.argv[1])
sys.stderr.write(" DONE\n")


	
if (data.shape[1]-1)%n != 0:
	print "couldn't find n={} neurons in the file".format(n,sys.argv[1])
	exit(1)

nvars = (data.shape[1]-1)/n
Nplots = (nvars + 1)*100+10


f=plt.figure(1)
ax=plt.subplot(Nplots+1)
spks = []
for i in xrange(n):
	print data.shape, i, i*nvars+1
	spk = where(data[:,i*nvars+1] > 0.0)[0]
	if len(spk) > 0:
		spk = [ (data[spk[0],0],i) ] + [ (data[k,0],i) for p,k in zip(spk[:-1],spk[1:]) if p != k-1  ]
		spks += spk
spks = array(spks)
plt.plot(spks[:,0], spks[:,1], "k.")

neuron = 0
curves = []
for nv in xrange(nvars):
	plt.subplot(Nplots+2+nv,sharex=ax)
	curves.append( plt.plot(data[:,0],data[:,neuron*nvars+nv+1])[0] )

def kpress(event):
	global neuron,curves
	if event.key == "up":
		neuron += 1
		if neuron >= n : neuron =0
	elif event.key == "down":
		neuron -= 1
		if neuron < 0 : neuron = n-1
	for nv,c in enumerate(curves):
		c.set_ydata( data[:,neuron*nvars+nv+1] )
	f.canvas.draw()

f.canvas.mpl_connect('key_press_event',kpress)

	
plt.show()
