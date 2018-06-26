'''
this code generates weights for weights.data
It first writes weights for all-to-all connected and then optionally changes to sparsely connected, but preserving the average number of connections per neuron.
'''
n = 20
weightstr = ''
wl = [] #weights list
for i in range(n):
	for j in range(n):
		if i == j:
			weightstr += '0. '
		else:
			weightstr += '1. '
	weightstr += '\n'
print(weightstr)
def pw(wl): #prints the weights
	for i in le
	return
for i in range(n):
	for j in range(n):
		if i == j:
			wl += ['0. ']
		else:
			wl += ['1. ']. '
	wl += ['\n']
	#wl should be list of lists, simple fix, but no time now
