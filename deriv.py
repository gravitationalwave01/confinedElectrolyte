#this code computes the derivative of (x,y) data stored in the file given by fname.
#The code uses a simple central difference method

import numpy as np

fin = "total_energy.dat"
fout = "press.dat"
column = 1
deg = 40 #degree of polynomial fit

data = np.loadtxt(fin,skiprows=1)
output = np.zeros((data.shape[0]-1,3))
myfit=np.polyfit(data[:,0],data[:,1],deg)	
print(myfit)

for i in xrange(data.shape[0]-1):
	#compute derivative manually
	output[i,1] = (data[i+1,column]-data[i,column])/(data[i+1,0]-data[i,0])
	output[i,0] = 0.5*(data[i+1,0]+data[i,0])
	
	#compute derivative from polynomial fit
	for j in xrange(deg+1):
		output[i,2] = output[i,2] + (deg-j)*myfit[j]*(output[i,0]**(deg-1-j))

np.savetxt(fout,output)

