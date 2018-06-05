from numpy import *
from matplotlib import pyplot as plt
import sys
v=genfromtxt(sys.argv[1])
plt.plot(v[:,0], v[:,1],"k-o")
plt.show()
