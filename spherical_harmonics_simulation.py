import numpy as np
import matplotlib.pyplot as plt 

from pyshtools.rotate import djpi2, SHRotateRealCoef
from pyshtools.expand import SHExpandDH
from pyshtools import SHCoeffs,SHGrid

from trace_generator import generate_steps,writeout

np.random.seed(0)

def save_pattern_figures(grid,label):
	grid = SHGrid.from_array(grid)
	fig,ax = grid.plot3d()
	plt.savefig("images/patterns/pattern{}_3d.png".format(label))
	plt.close()
	fig,ax = grid.plot()
	plt.savefig("images/patterns/pattern{}_2d.png".format(label))

def pattern_generator(nx,pattern):
	ny = nx
	grid = np.zeros((nx,nx))
	if pattern == 0:
		grid[0:nx/2,0:nx/2] = 1
		grid[nx/2:,nx/2:] = 1
	elif pattern == 1:
		for i in range(0,4):
			for j in range(0,4):
				if i % 2 == 0 and j%2 == 1:
					grid[i*nx/4:(i+1)*nx/4,j*nx/4:(j+1)*nx/4] = 1
				elif i % 2 == 1 and j%2 == 0:
					grid[i*nx/4:(i+1)*nx/4,j*nx/4:(j+1)*nx/4] = 1
	elif pattern == 2:
		for i in range(0,8):
			for j in range(0,8):
				if i % 2 == 0 and j%2 == 1:
					grid[i*nx/8:(i+1)*nx/8,j*nx/8:(j+1)*nx/8] = 1
				elif i % 2 == 1 and j%2 == 0:
					grid[i*nx/8:(i+1)*nx/8,j*nx/8:(j+1)*nx/8] = 1
	elif pattern == 3:
		for i in range(0,8):
			if i % 2 == 0:
				grid[:,i*nx/8:(i+1)*nx/8] = 1
	elif pattern == 4:
		for i in range(0,8):
			if i % 2 == 0:
				grid[i*nx/8:(i+1)*nx/8,:] = 1
	return grid

def calculate_intensity(grid,theta):
	''' 
		grid - 2D scattering pattern, coordinates (theta,phi)
		theta range: [-np.pi/2:np.pi/2]
		phi range: [0:2*np.pi]

	Assume we are observing from the top of the particle.
	
	Collection (NA) is defined by our collection half-angle (theta)
		This defined the offset from the array midpoint:
			midpoint = grid.shape[0]/2 

	In the phi direction we integrate over all phi so our ROI extends from 0:grid.shape[1] in axis=1 	
	'''

	offset = int(np.round(theta*(grid.shape[0]/180.0)))
	midpoint = int(grid.shape[0]/2)

	roi = grid[midpoint-offset:midpoint+offset,:]
	total = np.sum(roi,axis=(0,1))

	return total

def rotate(coefs,tx,ty,tz,lmax):
	#Perform spherical harmonics coefficient rotation
	dj = djpi2(lmax)
	coefs = SHRotateRealCoef(cilm=coefs,x=[tx,ty,tz],dj=dj)
	return coefs

def orientation_scattering(pattern,rotation,lmax,theta_NA):
	#get rotation angles
	[tx,ty,tz] = rotation

	#expand pattern in terms of spherical harmonic coefficients
	coefs = SHExpandDH(pattern)
	#rotate pattern. lmax is the maximum order of the pattern expansion (= number of fourier components)
	coefs = rotate(coefs,tx,ty,tz,lmax)

	#convert back to a 2D array 
	grid = SHCoeffs.from_array(coefs).expand().to_array()
	#calculate intensity by summing over region of interest in 2D grid
	outp =  calculate_intensity(grid,theta_NA)

	return outp

def intensity_calculation(args):
	#compute scattered intensity by 
	# print args
	(index,scattering_pattern,angular_trace,theta) = args
	N = scattering_pattern.shape[0]
	# print angular_trace.shape
	return orientation_scattering(pattern=scattering_pattern,rotation=angular_trace[index,:],lmax=N,theta_NA=theta)

def single_particle_simulation(theta_step,theta_NA):
	N_steps = 10000
	
	trace = generate_steps(linear_stepsize=0,angular_stepsize=np.pi/20.0,N_steps=N_steps)
	
	writeout(trace,"traces/trace.npy")

	angular_trace = trace[:,3:6]
	patterns = [ np.array(pattern_generator(100,i)) for i in range(0,5)]
	for i,pat in enumerate(patterns):
		writeout(pat,"traces/pattern{}.npy".format(i))
		# save_pattern_figures(pat,i)

	intensities = np.zeros((N_steps,len(patterns)))
	from multiprocessing import Pool
	pool = Pool(3)

	for i,pattern in enumerate(patterns):
		print "Pattern",i
		args = [ (j,pattern,angular_trace,theta_NA) for j in range(N_steps)]
		# print args[0]
		# intensity_p = [intensity_calculation(a) for a in args] #for DEBUG PURPOSES
		intensity_p = np.array(pool.map(intensity_calculation,args))
		intensities[:,i] = np.array(intensity_p).T 
	
	return intensities
intensities = single_particle_simulation(theta_step=np.pi/20.0,theta_NA=50.0)
writeout(intensities,"traces/intensities.npy")
