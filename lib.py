import numpy as np
import matplotlib.pyplot as plt

def verify_vector_orthogonality(particles):
	#for testing that the vectors defining the particles internal coordinate system are still orthogonal after transformation
	N_particles = particles.shape[0]
	for i in range(N_particles):
		p = particles[i,:]
		[tx,ty,tz] = particles[i,3:6]

		rx = [1,0,0]
		ry = [0,1,0]
		rz = [0,0,1]

		rx = rotate(rx,tx,ty,tz)
		ry = rotate(ry,tx,ty,tz)
		rz = rotate(rz,tx,ty,tz)

		#orthogonal:
		assert(np.dot(rx,ry) < 1e-10)
		assert(np.dot(ry,rz) < 1e-10)
		assert(np.dot(rx,rz) < 1e-10)
		#unit dot product
		assert(np.dot(rx,rx)-1.0 < 1e-10)
		assert(np.dot(ry,ry)-1.0 < 1e-10)
		assert(np.dot(rz,rz)-1.0 < 1e-10)
		
	return 

def rotate(orientation,tx,ty,tz):
	#rotation of arrays of vectors by angles:
	#	theta_x - around x-axis
	#   theta_y - around y-axis
	#	theta_z - around z-axis

	orientation = np.array(orientation)
	thetas = [tx,ty,tz]
	cx,sx = np.cos(tx),np.sin(tx)
	cy,sy = np.cos(ty),np.sin(ty)
	cz,sz = np.cos(tz),np.sin(tz)

	Rx = np.array([[1,0,0],[0,cx,-sx],[0,sx,cx]])
	Ry = np.array([[cy,0,sy],[0,1,0],[-sy,0,cy]])
	Rz = np.array([[cz,-sz,0],[sz,cz,0],[0,0,1]])
	Rs = [Rx,Ry,Rz]
	for R in Rs:
		orientation = np.dot(R.T,orientation)

	return orientation


def test0(debug):
	fig, ax = plt.subplots(1)
	vx = np.array([1,0,0])
	ax.plot()
	print "before:", vx
	vx = rotate(vx,0,np.pi,0)

	plt.show()
	print "after:", vx
if __name__ == "__main__":

	test0(0)