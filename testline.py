import numpy as np
from matplotlib import pyplot as plt
import sys
v=np.genfromtxt(sys.argv[1])
plt.plot(v[:,0], v[:,1],"k-*") 
x = (np.linspace(0.0, 1.0, num=201)).tolist()
print x
ydict = {}
y = []
for i in range(len(x)):
	ydict['x'] = x[i]
	y.append(eval('(2 * x) - 1'), ydict)
print ydict
#~ y = eval('2x-1')
#~ plt.plot(x, y)  
#~ plt.show()
