import sys
import random

randomize = False
self = False
allon = True
alloff = False

f = open(sys.argv[2], 'w')

if randomize:
	for i in range(int(sys.argv[1])):
		for j in range(int(sys.argv[1]) - 1):
			if ((not self) and i == j):
				f.write('0. ')
			elif random.choice([True, False]):
				f.write('1. ')
			else:
				f.write('0. ')
		if ((not self) and i == (int(sys.argv[1]) - 1)):
				f.write('0.\n')
		elif random.choice([True, False]):
			f.write('1.\n')
		else:
			f.write('0.\n')
else:
	if allon:
		for i in range(int(sys.argv[1])):
			for j in range(int(sys.argv[1]) - 1):
				if ((not self) and i == j):
					f.write('0. ')
				else:
					f.write('1. ')			
			if ((not self) and i == (int(sys.argv[1]) - 1)):
				f.write('0.\n')
			else:
				f.write('1.\n')
	elif alloff:
		for i in range(int(sys.argv[1])):
			for j in range(int(sys.argv[1]) - 1):
					f.write('0. ')
			f.write('0.\n')
f.close()
