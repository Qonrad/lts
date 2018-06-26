'''
this code generates weights for weights.data
It first writes weights for all-to-all connected and then optionally changes to sparsely connected, but preserving the average number of connections per neuron.
'''
n = 20
weightstr = ''
for i in range(n):
	for j in range(n):
		if i == j:
			weightstr += '0. '
		else:
			weightstr += '1. '
	weightstr += '\n'
print(weightstr)
