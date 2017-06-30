import numpy as np
from matplotlib import pyplot as plt
import sys
f, ax1= plt.subplots()
v=np.genfromtxt(sys.argv[1])

ax1.plot(v[:,0], v[:,1],"k-") 
ax1.plot([0.5,1.],[0.,1.],"r--")

ax2 = ax1.twinx()
p=np.genfromtxt(sys.argv[2])
x=(p[:,0] - p[0,0])/(p[-1,0]-p[0,0])
ax2.plot(x, p[:,1],"r-") 

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
