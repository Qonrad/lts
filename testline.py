import numpy as np
from matplotlib import pyplot as plt
import sys
print "format is python testline.py prc1.data prc2.data delay period"

f, ax1= plt.subplots()
v=np.genfromtxt(sys.argv[1]) #v is a 2 * (INTERVAL + 1) array of the prc data, [[0, 0.2], [0.01, 0.21]...[1, 0.9]]
vs=np.genfromtxt(sys.argv[2]) #vs is the same except for the prc2.data
delay = float(sys.argv[3])
period = float(sys.argv[4])

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
def secprc(phi):
	if np.any(phi) > 1 or np.any(phi) < 0:
		raise ValueError("points were given to prc < 0 or > 1")
	idx0 = np.where(vs[:,0] < phi)[0][-1]
	idx1 = idx0 + 1
	x0 = vs[idx0][0]
	x1 = vs[idx1][0]
	y0 = vs[idx0][1]
	y1 = vs[idx1][1]
	return (y0 * (x1 - phi) + y1 * (phi - x0)) / (x1 - x0)
def linterp(x, x0, y0, x1, y1):
	return (y0 * (x1 - x) + y1 * (x - x0)) / (x1 - x0)
def slp(one, two):
	return ((two[1] - one[1]) / (two[0] - one[0]))
def fullterp(preline, postline, preprc, postprc): #given 4 arbitrary x and y coordinates, returns linearly interpolated points of intersection
	#print "points are", preline, "and", preprc, "before intersection"
	#print "poins are", postline, "and", postprc, "after intersection"
	lineslp = slp(preline, postline)
	prcslp = slp(preprc, postprc)
	#print "lineslp =", lineslp
	#print "prcslp =", prcslp
	x = ((lineslp * preline[0]) - preline[1] - (prcslp * preprc[0]) + preprc[1]) / (lineslp - prcslp)
	#print "x should be", x
	interpoint = (x, linterp(x, preprc[0], preprc[1], postprc[0], postprc[1]))
	return interpoint
def intersect(xvals, yset1, yset2): #finds intersection points between two sets of y vals with the same x values, returns tuple of arrays
	intersidx = np.where(np.diff(np.sign(yset2 - yset1)) != 0)[0]
	inters = fullterp((xvals[intersidx], yset1[intersidx]), (xvals[intersidx + 1], yset1[intersidx + 1]), (xvals[intersidx], yset2[intersidx]), (xvals[intersidx + 1], yset2[intersidx + 1]))
	return np.asarray((inters[0], inters[1], slopeonprc(intersidx)))
def slopeonprc(preidx):
	x0 = v[:,0][preidx]
	x1 = v[:,0][preidx + 1]
	y0 = v[:,1][preidx]
	y1 = v[:,1][preidx + 1]
	return (y1 - y0) / (x1 - x0)

np.vectorize(prc)
np.vectorize(linterp)
np.vectorize(slopeonprc)

#making analysis sets
y = (2. * v[:,0]) - 1. - (2. * (delay / period))									#line for 2-cluster stability
ys = (2. * v[:,0]) - 1. - (2. * (delay / period)) + vs[:,1]							#line 2.0 that takes second order resetting into account
bu = prc(1 - v[:,0] + v[:,1] + (2. * delay / period))								#basin of attraction for unequal time lags, formula f(1 - phi + f(phi) + (2 * delay / period))

#printing explanations
print "\nGREEN | Basin of attraction of unequal time lags | f(1 - phi + f(phi) + (2 * delay / period))"
print "  RED | 2-cluster mode stability w/out 2nd order | (2 * phi) - 1 - (2 * (delay / period))"
print "BLUE  | 2-cluster mode stability WITH 2nd order  | (2 * phi) - 1 - (2 * (delay / period)) + f2(phi)"

#finding points of intersection
yinters = intersect(v[:,0], v[:,1], y)
ysinters = intersect(v[:,0], v[:,1], ys)
buinters = intersect(v[:,0], v[:,1], bu)

#printing them extra nicely
print "\nAntiphase line w/out 2nd order intersects at point(s)"
for i in range(len(yinters[0])):
	print yinters[:2,i], "with PRC slope", yinters[2][i]
print "\nAntiphase line with 2nd order intersects at point(s)"
for i in range(len(ysinters[0])):
	print ysinters[:2,i], "with PRC slope", ysinters[2][i]
print
#restricting analysis sets to reasonable points before graphing
ycut = np.where(y >= np.amin(v[:,1]))[0][0]											#cutting off y so that it doesn't go below the min value of f(phi)
yscut = np.where(ys >= np.amin(v[:,1]))[0][0]
"""
intersidxbu = np.where(np.diff(np.sign(bu - v[:,1])) != 0)[0]						#hopefully index of array holding line unequal time lags immediately BEFORE intersection w/ prc
intersbu = linterp((intersidxbu + intersidxbu + 1) / 2., intersidxbu, y[intersidxbu], intersidxbu + 1, y[intersidxbu + 1])
intersidx = np.where(np.diff(np.sign(y - v[:,1])) != 0)[0]							#index of array holding line 2-cluster stability immediately BEFORE intersection w/ prc
inters = fullterp((v[:,0][intersidx], y[intersidx]), (v[:,0][intersidx + 1], y[intersidx + 1]), (v[:,0][intersidx], v[:,1][intersidx]), (v[:,0][intersidx + 1], v[:,1][intersidx + 1]))
print "intersection points", inters

print "y - v[:,1]", y - v[:,1]
print "np.sign(y - v[:,1])", np.sign(y - v[:,1])
print "np.diff(np.sign(y - v[:,1]))", np.diff(np.sign(y - v[:,1]))
print "np.where(np.diff(np.sign(y - v[:,1])) != 0)", np.where(np.diff(np.sign(y - v[:,1])) != 0)
"""
"""
#inters = linterp((intersidx + intersidx + 1) / 2., intersidx, y[intersidx], intersidx + 1, y[intersidx + 1])
intersidxs = np.where(np.diff(np.sign(ys - v[:,1])) != 0)[0]
interss = linterp((intersidxs + intersidxs + 1) / 2., intersidxs, ys[intersidxs], intersidxs + 1, ys[intersidxs + 1])
print "intrinsic period is", period
#print "line for 2-cluster stability, using formula (2 * phi - 1 - (2 * (delay/period)) + f2(phi), is detected to intersect with PRC at these sets of coordinates:\n", np.dstack((v[:,0][intersidx], inters))
ycut = np.where(y >= np.amin(v[:,1]))[0][0]											#cutting off y so that it doesn't go below the min value of f(phi)
yscut = np.where(ys >= np.amin(v[:,1]))[0][0]
slope = slopeonprc(intersidx)														#approximate slope of the prc on the point that it intersects with the line for 2-cluster stability
slopes = slopeonprc(intersidxs)

print "slope of prc where it intersects with line for 2-cluster stability\n", slope
print "curve is formula f(1 - phi + f(phi) + (2 * delay / period))"
print "curve intersects with PRC at these sets of coordinates:\n", np.dstack((v[:,0][intersidxbu], intersbu))
#print "v[0][1] = f(0) is", v[0][1], "v[len(v) - 1][1] = f(1) is", v[len(v) - 1][1]
print "blue dotted thing is antiphase stabliity with second order accounted for"
print "it intersects with prc at these sets of coords", np.dstack((v[:,0][intersidxs], interss))
print "slope of prc where it intersects with line 2.0\n", slopes

print np.arange(0, len(v[:,0]) - 1, 1)
all = slopeonprc(np.arange(0, len(v[:,0]) - 1, 1))
print "all slopes\n", all
max = np.amax(all)
print "max", max


print "intersidx", intersidx
print "v[:,0][intersidx]", v[:,0][intersidx]
print "((intersidx + intersidx + 1) / 2.) / (len(y) - 1)", ((intersidx + intersidx + 1) / 2.) / (len(y) - 1)
print "y[intersidx]", y[intersidx]
print "inters", inters
"""
#synchrony calculations
sync = delay/period																	#x coordinate of the point for synchrony stability
if delay != 0:
	syncnum = prc(sync)																#y coordinate of the point for synchrony stability
else:
	syncnum = v[0][1]
syncwd = 1 - v[0][1] - v[len(v) - 1][1]
print "testing -1 < 1 - f(0) - f(1) < 1"
print "-1 <", syncwd, "< 1"
if -1 < syncwd and syncwd < 1:
	print "yes, so synchrony with delay should be stable"
else:
	print "no, so synchrony should not be stabilized by delay"
print "testing -1 < (1 - f(0)) * (1- f(1)) < 1"
syncwnd = (1 - v[0][1]) * (1 - v[len(v) - 1][1])
print "-1 <", syncwnd, "< 1"
if -1 < syncwnd and syncwnd < 1:
	print "yes, so synchrony with no delay should be stable"
else:
	print "no, so synchrony without delay should not be stable"


#plotting
ax1.plot(v[:,0], v[:,1],"k-") 														#main prc
ax1.plot(sync, syncnum, 'ro')														#point of synchrony stability on the prc
ax1.plot(v[:,0], bu, "g--") 														#drawing line for basin of attraction for unequal time lags
ax1.plot(v[yscut:,0], ys[yscut:], "b--")
ax1.plot(v[ycut:,0], y[ycut:], "r--")												#line for 2-cluster stability
ax1.plot(yinters[0], yinters[1], 'ro')												#intersection of line for 2-cluster stability with prc
ax1.plot(ysinters[0], ysinters[1], 'bo')											#intersection of line 2.0 with prc
plt.ylim(-0.2, 1.0)
plt.show()
