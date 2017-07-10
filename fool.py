from numpy import *
from matplotlib.pyplot import *

d00 = genfromtxt("18full.data")
d01 = genfromtxt("1full.data")
#~ d02 = genfromtxt("no.data")
#~ d03 = genfromtxt("3.data")
#~ dXY = genfromtxt("unpertfull.data")

for nt,idx in enumerate(xrange(d00.shape[1]-1)):

	if nt == 0: ax=subplot(711)
	else      : subplot(710+nt+1, sharex=ax)

	#~ plot(dXX[:,0], dXX[:,nt+1],"k-",label="X.X",lw=3)
	#~ plot(d02[:,0], d02[:,nt+1],"y-",label="no delay", lw=3)
	plot(d00[:,0], d00[:,nt+1],"r-",label="neuron 18 in 20")
	plot(d01[:,0], d01[:,nt+1],"b-",label="neuron 1 in 2")
	#~ plot(d03[:,0], d03[:,nt+1],"g-",label="3ms")

legend(loc=0)
show()


