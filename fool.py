from numpy import *
from matplotlib.pyplot import *

d00 = genfromtxt("0full.data")
d01 = genfromtxt("1full.data")
#~ d05 = genfromtxt("0.5full.data")
#~ d06 = genfromtxt("0.6full.data")
#~ dXY = genfromtxt("unpertfull.data")

for nt,idx in enumerate(xrange(d00.shape[1]-1)):

	if nt == 0: ax=subplot(711)
	else      : subplot(710+nt+1, sharex=ax)

	#~ plot(dXX[:,0], dXX[:,nt+1],"k-",label="X.X",lw=3)
	#~ plot(d05[:,0], d05[:,nt+1],"y-",label="0.0")
	plot(d00[:,0], d00[:,nt+1],"r-",label="0")
	plot(d01[:,0], d01[:,nt+1],"b-",label="1")
	#~ plot(dXY[:,0], dXY[:,nt+1],"g-",label="X.Y")

legend(loc=0)
show()


