from numpy import random as rnd
from numpy import *

class random_state:
	def __init__(self, defmin, defmax, defdist):
		self.defmin = defmin
		self.defmax = defmax
		self.defdist = defdist
	def output(self, rmin=self.defmin, rmax=self.defmax, dist=self.defdist):
		if dist == "normal":
			return "normal"
		elif dist == "standard":
			return "standard"
		else:
			print("something is wrong")
			return -1
v = random_state(0, 100, "normal")
print(v.defmin, v.defmax, v.defdist)
