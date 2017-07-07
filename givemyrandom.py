from numpy import random as rnd
from numpy import *
import sys

if len(sys.argv) < 2:
	print "USAGE {} number_of_random variables [range] [range] ..."
	exit(1)

n = int(sys.argv[1])
#print [ x for x in sys.argv[2:] ]
vrange = [ eval(x) for x in sys.argv[2:] ]

for i in xrange(n):
	if type(vrange[i%len(vrange)] )is int or type(vrange[i%len(vrange)] )is float:
		print rnd.random()*vrange[i%len(vrange)]
	elif type(vrange[i%len(vrange)] )is list or type(vrange[i%len(vrange)] )is tuple:
		l,r = vrange[i%len(vrange)]
		print rnd.random()*(r-l)+l

