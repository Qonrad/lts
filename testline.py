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
print x # x is basically phi, v[0], and v[1] is basically f(phi)
y = (2. * x) - 1. - (2. * (delay / period))	#making line for 2-cluster stability
print y
print type(y)
print v
print type(v)
z = x + 2 * delay / period
print z
print "testing thing"
# w is storing basin of attraction, formula f(phi) = f(1 - phi + f(phi + 2 * delay / period))
w = y.copy()
for i in range(len(x)):
	#print round(z[i], 2), v[i][0]
	for j in range(len(x)):
		if round(z[i], 2) == v[j][0]:
			print round(z[i], 2), "adding", v[j][1], "to w"
			w[i] = v[j][1]
print w
# w should now contain f(phi + 2 * delay / period)
w = 1 - x + w
print w
# w should now contain 1 - phi + f(phi  + 2 * delay / period)
for i in range(len(x)):
	#print round(z[i], 2), v[i][0]
	for j in range(len(x)):
		if round(w[i], 2) == v[j][0]:
			print round(w[i], 2), "adding", v[j][1], "to w"
			w[i] = v[j][1]
# w should now contain f(1 - phi + f(phi + 2 * delay / period))
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
#now trying to add automatic basin of attraction functionality, this should be fun

#~ while rem >= 0:
	#~ rem = int(raw_input("If there is a mistake, enter the index of a point to remove it, else enter -1\n"))
	#~ if rem >= 0:
		#~ idx = np.delete(idx, rem)
#~ for i in range(490, 495):
	#~ print i, v[:,0][i], y[i], v[:,1][i]
ax1.plot(v[:,0], v[:,1],"k-")
ax1.plot(x[idx], y[idx], 'ro')
ax1.plot(sync, v[:,1][syncnum], 'ro')
ax1.plot(x, w, "g--") #drawing line for basin of attraction
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
