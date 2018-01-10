import numpy as np
from matplotlib import pyplot as plt
import sys
print "format is python testline.py prc.data delay period"	

f, ax1= plt.subplots()
v=np.genfromtxt(sys.argv[1]) #v is a 2 * (INTERVAL + 1) array of the prc data, [[0, 0.2], [0.01, 0.21]...[1, 0.9]]
delay = float(sys.argv[2])
period = float(sys.argv[3])

"""
2do
maybe change y to use linspace for input instead of phi? better intersection point?
I definitely don't think I'm using the best linear interpolation for the instersection
But everything else should be ok. The basin of attraction is still strange...

make it work with delay 0

also ought to implement system arguments into main code

fix antiphase basin of attraction
"""
def prc(phi):
	if np.any(phi) > 1 or np.any(phi) < 0:
		raise ValueError("points were given to prc < 0 or > 1")
	idx0 = np.where(v[:,0] < phi)[0][-1]
	idx1 = idx0 + 1
	x0 = v[idx0][0]
	x1 = v[idx1][0]
	y0 = v[idx0][1]
	y1 = v[idx1][1]
	return (y0 * (x1 - phi) + y1 * (phi - x0)) / (x1 - x0)
def linterp(x, x0, y0, x1, y1):
	return (y0 * (x1 - x) + y1 * (x - x0)) / (x1 - x0)
def slopeonprc(preidx):
	x0 = v[:,0][preidx]
	x1 = v[:,0][preidx + 1]
	y0 = v[:,1][preidx]
	y1 = v[:,1][preidx + 1]
	return (y1 - y0) / (x1 - x0)
	"""
def prcintersx(preidx, otherline):			#assumes otherline has same length as prc file (which it should)
	k1 = (otherline[preidx + 1] - otherline[preidx])
	"""
np.vectorize(prc)
np.vectorize(linterp)
np.vectorize(slopeonprc)

#calculations
y = (2. * v[:,0]) - 1. - (2. * (delay / period))									#line for 2-cluster stability
bu = prc(1 - v[:,0] + v[:,1] + 2. * delay / period)									#basin of attraction for unequal time lags, formula f(1 - phi + f(phi) + 2 * delay / period 
sync = delay/period																	#x coordinate of the point for synchrony stability
syncnum = prc(sync)																	#y coordinate of the point for synchrony stability
intersidx = np.where(np.diff(np.sign(y - v[:,1])) != 0)[0]							#index of array holding line 2-cluster stability immediately BEFORE intersection w/ prc
inters = linterp((intersidx + intersidx + 1) / 2., intersidx, y[intersidx], intersidx + 1, y[intersidx + 1])
print "line for 2-cluster stability is detected to intersect with prc at these points:\n", np.dstack((v[:,0][intersidx], inters))
ycut = np.where(y >= np.amin(v[:,1]))[0][0]											#cutting off y so that it doesn't go below the min value of f(phi)
slope = slopeonprc(intersidx)														#approximate slope of the prc on the point that it intersects with the line for 2-cluster stability
print "slope of prc where it intersects with line for 2-cluster stability\n", slope

"""
print np.arange(0, len(v[:,0]) - 1, 1)
all = slopeonprc(np.arange(0, len(v[:,0]) - 1, 1))
print "all slopes\n", all
max = np.amax(all)
print "max", max
"""

#plotting
ax1.plot(v[:,0], v[:,1],"k-") 														#main prc
ax1.plot(((intersidx + intersidx + 1) / 2.) / (len(y) - 1), inters, 'ro')			#intersection of line for 2-cluster stability with prc
ax1.plot(sync, syncnum, 'ro')														#point of synchrony stability on the prc
ax1.plot(v[:,0], bu, "g--") 														#drawing line for basin of attraction for unequal time lags
ax1.plot(v[ycut:,0], y[ycut:], "r--")												#line for 2-cluster stability
plt.show()
