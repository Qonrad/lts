import numpy as np
from numpy import random as rnd
from numpy import *

n = 20 #number of total neurons in output

class random_state:
	def __init__(self, defmin, defmax, defdist, defval=0):
		self.defmin = defmin
		self.defmax = defmax
		self.defdist = defdist
		self.defval = defval
	def output(self, n, parmin=None, parmax=None, pardist=None):
		
		#did this silliness because default value of parameters can't be self.whatever for some reason
		#there is probably a much better way to do this
		if parmin is None:
			actmin = self.defmin
		else:
			actmin = parmin
		if parmax is None:
			actmax = self.defmax
		else:
			actmax = parmax
			
		if pardist is None:
			actdist = self.defdist
		else:
			actdist = pardist

		#it has been decided what parameters will be used for this output, printing them here	
		#print("min =", actmin)
		#print("max =", actmax)
		#print("dist =", actdist)

		if actdist == "normal": #if the dist is normal, we should output a random value within the range based on a normal distribution



			return ((rnd.randn(1, n) * ((actmax - actmin) / 2)) + actmin)
		elif actdist == "uniform": #we should do the same thing with standard



			return ((rnd.rand(1, n) * (actmax - actmin)) + actmin)
		elif actdist == "val":
			return np.full((1, n), self.defval)
		else:
			print("something is wrong")
			return -1

v = random_state(-100, 50, "uniform")
m = random_state(0, 1, "uniform")
h = random_state(0, 1, "uniform")
nv = random_state(0, 1, "uniform")
mn = random_state(0, 1, "val", defval=0)
s = random_state(0, 1, "uniform")
p = random_state(0, 1, "val", defval=0)

vout = v.output(n)
mout = m.output(n)
hout = h.output(n)
nvout = nv.output(n)
mnout = mn.output(n)
sout = s.output(n)
pout = p.output(n)
#print(v.defmin, v.defmax, v.defdist)
#print(v.output(20))

for i in range(n):
	print(vout[0][i])
	print(mout[0][i])
	print(hout[0][i])
	print(nvout[0][i])
	print(mnout[0][i])
	print(sout[0][i])
	print(pout[0][i])
