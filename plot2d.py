import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from lib import rotate

def draw_particle(ax,particle,vector_length=1e-6):
	[x,y,z] = particle[0:3]
	[tx,ty,tz] = particle[3:6]

	#unit vectors in each dimension:
	rx = [1,0,0]
	ry = [0,1,0]
	rz = [0,0,1]

	rx = rotate(rx,tx,ty,tz)*vector_length
	ry = rotate(ry,tx,ty,tz)*vector_length
	rz = rotate(rz,tx,ty,tz)*vector_length

	ax.plot([x],[y],[z],"ob")
	ax.plot([x,x+rx[0]],[y,y+rx[1]],[z,z+rx[2]],"red")
	ax.plot([x,x+ry[0]],[y,y+ry[1]],[z,z+ry[2]],"green")
	ax.plot([x,x+rz[0]],[y,y+rz[1]],[z,z+rz[2]],"blue")


def draw_particles(particles,index=None):
	print "Drawing representation..."
	N_particles = particles.shape[0]
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')

	for i in range(N_particles):
		particle = particles[i,:]
		draw_particle(ax,particle)
	ax.set_title("Red: [1,0,0] vector"+"\nGreen: [0,1,0] vector"+"\nBlue: [0,0,1] vector")
	ax.set_xlim(-1e-5,1e-5)
	ax.set_ylim(-1e-5,1e-5)
	ax.set_zlim(-1e-5,1e-5)

	plt.savefig("animation/simulation{0}.png".format(index))
	plt.close("all")
	# plt.show()