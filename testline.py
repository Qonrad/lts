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
implement checking for slope of 2-cluster stability line
maybe change y to use linspace for input instead of phi? better intersection point?
"""
def prc(phi):
	if np.any(phi) > 1 or np.any(phi) < 0:
		print "prc was given values < 0 or > 1"
		return -1
	idx0 = np.where(v[:,0] < phi)[0][-1]
	idx1 = idx0 + 1
	x0 = v[idx0][0]
	x1 = v[idx1][0]
	y0 = v[idx0][1]
	y1 = v[idx1][1]
	return (y0 * (x1 - phi) + y1 * (phi - x0)) / (x1 - x0)
def linterp(x, x0, y0, x1, y1):
	return (y0 * (x1 - x) + y1 * (x - x0)) / (x1 - x0)

np.vectorize(prc)
np.vectorize(linterp)

#calculations
y = (2. * v[:,0]) - 1. - (2. * (delay / period))									#line for 2-cluster stability
b = prc(1 - v[:,0] + v[:,1] + 2. * delay / period)									#basin of attraction for antiphase stability
sync = delay/period																	#x coordinate of the point for synchrony stability
syncnum = prc(sync)																	#y coordinate of the point for synchrony stability
intersidx = np.where(np.diff(np.sign(y - v[:,1])) != 0)[0]							#index of array holding line 2-cluster stability immediately BEFORE intersection w/ prc
inters = linterp((intersidx + intersidx + 1) / 2., intersidx, y[intersidx], intersidx + 1, y[intersidx + 1])
print "line for 2-cluster stability is detected to intersect with prc at these points:\n", inters
ycut = np.where(y >= np.amin(v[:,1]))[0][0]											#cutting off y so that it doesn't go below the min value of f(phi)

#plotting
ax1.plot(v[:,0], v[:,1],"k-") 														#main prc
ax1.plot(((intersidx + intersidx + 1) / 2.) / (len(y) - 1), inters, 'ro')			#intersection of line for 2-cluster stability with prc
ax1.plot(sync, syncnum, 'ro')														#point of synchrony stability on the prc
ax1.plot(v[:,0], b, "g--") 															#drawing line for basin of attraction
ax1.plot(v[ycut:,0], y[ycut:], "r--")												#line for 2-cluster stability
plt.show()


"""
print "best guess at intersection points"
print idx
rem = 0
"""
#print "v", v
#print prc(0.512)
#return
#~ print delay
#~ print period
#~ print v[:,0]
#~ exit(0)
#x = v[:,0]
#np.around(x,3) this rounding may be useful if there is a strange interval? but removing it for now
#sys.exit()
#print "x is a", type(x), "of size", np.size(x)
"""
print "v.shape = ", v.shape
print "v", v
print "v[:,0]", v[:,0]
#for i in range(len(x)):
#	x[i] = round(x[i], 3)
print "x", x # x is basically phi, v[0], and v[1] is basically f(phi)
"""


"""
print y
print type(y)
print v
print type(v)
"""
#z = x# + 2 * delay / period
#print z
#print "testing thing"
# w is storing basin of attraction, formula f(phi) = f(1 - phi + f(phi) + 2 * delay / period)
#w = y.copy()
"""
for i in range(len(x)):
	#print round(z[i], 2), v[i][0]
	for j in range(len(x)):
		if round(z[i], 2) == v[j][0]:
			print round(z[i], 2), "adding", v[j][1], "to w"
			w[i] = v[j][1]

print "z", z
for i in range(len(x)):
	idx = np.where(v[:,0] < z[i])[0][-1] #idx is an index in v where phi smaller than you need

w = v[idx,1] + (v[idx + 1, 1] - v[idx, 1] * (z[i] - v[idx.0]) / (v[idx + 1,0] - v[idx,0]
"""
#w = v[:,1] # w should now contain f(phi)
#print "w", w
#w = v[:,1] + 2 * delay / period
# w should now contain f(phi) + 2 * delay / period
#w = 1 - v[:,0] + v[:,1] + 2 * delay / period
#print w
# w should now contain 1 - phi + f(phi)  + 2 * delay / period
"""
for i in range(len(x)):
	#print round(z[i], 2), v[i][0]
	for j in range(len(x)):
		if round(w[i], 2) == v[j][0]:
			print round(w[i], 2), "adding", v[j][1], "to w"
			w[i] = v[j][1]
"""

# w should now contain f(1 - phi + f(phi) + 2 * delay / period)
#print "w done", w
"""
sync = round(delay / period, 3)
#print "sync = ", sync
syncnum = -999
if sync / 0.002 != 0.0:
	sync = round(sync / 0.002) * 0.002
	print "sync / 0.002 = ", sync / 0.002
	print "sync changed to ", sync
for i in range(len(v[:,0])):
	if v[:,0][i] == sync:
		print "found it"
		syncnum = i
if syncnum == -999:
	print "no synchrony found on prc"
"""

#now trying to add automatic basin of attraction functionality, this should be fun

#~ while rem >= 0:
	#~ rem = int(raw_input("If there is a mistake, enter the index of a point to remove it, else enter -1\n"))
	#~ if rem >= 0:
		#~ idx = np.delete(idx, rem)
#~ for i in range(490, 495):
	#~ print i, v[:,0][i], y[i], v[:,1][i]



"""
print "Values where two cluster line intersects with PRC"
print idx

for i in range(len(v[:,0])):		#removing the ones that are off screen
	if y[i] >= 0.0:
		y = y[i:]
		x = v[:,0][i:]
		break
"""

#~ ax1.plot(x[idx], y[idx], 'ro')
#~ print x
#~ print y
#~ x = np.arange(0.0, 1.0, 0.01)
#~ print x
#y =  

#~ ax2 = ax1.twinx()
#~ p=np.genfromtxt(sys.argv[2])
#~ x=(p[:,0] - p[0,0])/(p[-1,0]-p[0,0])
#~ ax2.plot(x, p[:,1],"r-") 


#~ x = (np.linspace(0.0, 1.0, num=201)).tolist()
#~ print x
#~ ydict = {}
#~ y = []
#~ for i in range(len(x)):
	#~ ydict['x'] = x[i]
	#~ y.append(eval('(2 * x) - 1'), ydict)
#~ print ydict
#~ y = eval('2x-1')
#~ plt.plot(x, y)  
#~ plt.show()

