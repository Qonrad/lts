import numpy as np
from matplotlib import pyplot as plt
import sys
print "format is python testline.py prc.data delay period"
f, ax1= plt.subplots()
v=np.genfromtxt(sys.argv[1])
delay = float(sys.argv[2])
period = float(sys.argv[3])

#~ print delay
#~ print period
#~ print v[:,0]
#~ exit(0)



x = np.linspace(0., 1., len(v[:,0]))
for i in range(len(x)):
	x[i] = round(x[i], 3)
print x
y = (2. * x) - 1. - (2. * (delay / period))	#making line for 2-cluster stability
sync = round(delay / period, 3)
print "sync = ", sync
if sync / 0.002 != 0.0:
	sync = round(sync / 0.002) * 0.002
	print "sync / 0.002 = ", sync / 0.002
	print "sync changed to ", sync
for i in range(len(v[:,0])):
	if v[:,0][i] == sync:
		print "found it"
		syncnum = i
idx = np.argwhere(np.diff(np.sign(y - v[:,1])) != 0).reshape(-1)
print "best guess at intersection points"
print idx
rem = 0
#~ while rem >= 0:
	#~ rem = int(raw_input("If there is a mistake, enter the index of a point to remove it, else enter -1\n"))
	#~ if rem >= 0:
		#~ idx = np.delete(idx, rem)
#~ for i in range(490, 495):
	#~ print i, v[:,0][i], y[i], v[:,1][i]
ax1.plot(v[:,0], v[:,1],"k-")
ax1.plot(x[idx], y[idx], 'ro')
ax1.plot(sync, v[:,1][syncnum], 'ro')
print "Values where two cluster line intersects with PRC"
print idx

for i in range(len(x)):		#removing the ones that are off screen
	if y[i] >= 0.0:
		y = y[i:]
		x = x[i:]
		break
ax1.plot(x, y, "r--")

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

plt.show()
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
