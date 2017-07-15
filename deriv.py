#this code computes the derivative of (x,y) data stored in the file given by fname.
#The code uses a simple central difference method

import numpy as np

fin = "total_energy.dat"
fout = "press.dat"
column = 1

data = np.loadtxt(fin,skiprows=1)
output = np.zeros((data.shape[0]-1,2))
for i in xrange(data.shape[0]-1):
	output[i,1] = (data[i+1,column]-data[i,column])/(data[i+1,0]-data[i,0])
	output[i,0] = 0.5*(data[i+1,0]+data[i,0])

np.savetxt(fout,output)

