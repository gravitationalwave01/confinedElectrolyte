import numpy as np


for i in xrange(0,100,1):
	sig = -10*1.6e-3 + i*0.1*1.6e-3
	fin="ion_profile_s"+str(i)+".dat"
	with open(fin) as f:
		data = np.loadtxt(fin,skiprows=1)
		last = data[data.shape[0]-1,3]
		first = data[0,3]
		print(str(sig) + " " + str(last-first) + " " + str(sig/(last-first)))

