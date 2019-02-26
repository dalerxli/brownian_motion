import matplotlib.pyplot as plt 
import numpy as np 
from trace_generator import loadin

from pyshtools.expand import SHExpandDH
from pyshtools import SHCoeffs,SHGrid


def normalise(xs):
	return (xs - np.min(xs))/(np.max(xs)-np.min(xs))

def plot_3d_pattern(pattern):
	pattern = SHGrid.from_array(pattern)
	return pattern.plot3d()
	# plt.show()

	# coefs = SHExpandDH(pattern)
	# g = coefs.expand()
	# g.plot3d(ax)

def plot_pattern(i):
	intensities = loadin("traces/intensities.npy")
	pattern = loadin("traces/pattern{}.npy".format(i))
	
	fig, ax = plot_3d_pattern(pattern)
	fig2, [ax1,ax2] = plt.subplots(2)
	ax1.imshow(pattern,cmap="gray")
	ax2.plot(normalise(intensities[:,i]))
	plt.show()


def plot_patterns(i,j):
	intensities = loadin("traces/intensities.npy")
	fig,ax = plt.subplots(1)
	ax.plot(normalise(intensities[:,i]),normalise(intensities[:,j]),'x')
	plt.show()

def plot_pattern_offset(i,offset):
	intensities = loadin("traces/intensities.npy")
	fig,ax = plt.subplots(1)
	N = intensities.shape[0]
	I = normalise(intensities[:,i])
	ax.plot(I,np.roll(I,offset),'x')
	plt.show()

# plot_pattern_offset(2,30)
plot_patterns(2,4)