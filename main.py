import numpy as np
np.random.seed(0)
from lib import rotate,verify_vector_orthogonality
import sys

from plot2d import draw_particles
DEBUG = 0


def initialize_particles(N_particles=10,box_dimension = 5e-6):
	particles = np.zeros((N_particles,6))
	for i in range(N_particles):
		#initialize the positions:
		particles[i,0:3] = np.random.uniform(-box_dimension,box_dimension,3)

		#orientation as direction of each unit vector (every particle is a coordinate system)
		
		# orientation = np.random.uniform(0,2*np.pi,3)
		orientation = np.random.uniform(0,0,3)
		
		particles[i,3:6] = orientation

	#Draw the initial positions of the particles:
	if DEBUG > 0: draw_particles(particles)
	return particles


def simulation_timestep(particles,step_size=10e-9,angular_step_size=(np.pi/20.0)):
	N_particles = particles.shape[0]

	#Generate unit vectors in which particle will translate
	step = np.random.uniform(-1,1,(N_particles,3))
	norms = np.linalg.norm(step,axis=1)
	for i in range(N_particles):
		step[i,:] = step[i,:]/norms[i]
	step = step*step_size

	#Step the particles:
	particles[:,0:3] = particles[:,0:3] + step
	
	#Generate the angles in which the direction vectors of each particle will be rotated:
	angular_step = np.random.uniform(-angular_step_size,angular_step_size,(N_particles,3))
	
	particles[:,3:6] = particles[:,3:6] + angular_step 
	verify_vector_orthogonality(particles)
	return particles

def phasor_sum(k,particles, wavevector_direction,detector_coordinates=np.array([0,0,10])):
	N_particles = particles.shape[0]
	#particle (x,y,z) coordinates
	positions = particles[:,0:3]
	#particle orientation (tx,ty,tz)
	angles = particles[:,3:6]

	#vector from particle to the detector:
	particle_to_detector = detector_coordinates-positions
	distances = np.linalg.norm(particle_to_detector,axis=1)
	unit_vector_to_detector = np.divide(particle_to_detector.T,distances).T

	scattering_directions = np.zeros((N_particles,3))
	for i in range(N_particles):
		scattering_directions[i,:] = rotate(wavevector_direction,angles[i,0],angles[i,1],angles[i,2])
	
	#compute the amplitude of the scattering by projecting the scattering directions from each particle onto unit vector in direction of detector
	amplitudes = np.array([np.dot(unit_vector_to_detector[i],scattering_directions[i]) for i in range(N_particles)])
	amplitudes[amplitudes < 0.0] = 0.0 #threshold to eliminate those particles whose projection is negative --> emits in different direction
	
	phasors = np.exp(-1j*np.dot(k,distances))
	phasor_sum = np.dot(amplitudes,phasors)
	intensity = np.real(np.conj(phasor_sum)*phasor_sum)
	
	if DEBUG > 0:
		print "distances:",distances.shape
		print "phasors:",phasors.shape
		print "amplitudes:",amplitudes.shape
		print "phasor_sum:", phasor_sum
		print "intensity:", intensity
		
	return intensity

def simulation_simple(wavevectors,N_steps=10000,N_particles=10,step_size=10e-9,angular_step_size=(np.pi/20.0),draw=False):
	particles = initialize_particles(N_particles) #this is our particles in the simulation
	intensity_trace = np.zeros( (len(wavevectors),N_steps) ) # this is a log of the intensity trace
	
	for i in range(N_steps):
		print "step",i
		particles = simulation_timestep(particles,step_size=step_size,angular_step_size=angular_step_size)
		if draw == True:
			draw_particles(particles,index=i)
		for j,(k,direction) in enumerate(wavevectors):
			intensity_trace[j,i] = phasor_sum(k,particles,direction)

	return intensity_trace

if __name__ == "__main__":
	N_particles = 1
	N_steps = 50
	step_size=50e-9,
	angular_step_size=0#(np.pi/10.0)
	v0 = np.array([1,0,0])
	v1 = np.array([-1,0,0])
	v1 = v1/np.linalg.norm(v1)

	wavevectors = [(2*np.pi/633e-9,v0),(2*np.pi/532e-9,v1)]
	trace = simulation(N_steps = N_steps,N_particles=N_particles,step_size=step_size,angular_step_size=angular_step_size,wavevectors=wavevectors,draw=True)
	import matplotlib.pyplot as plt
	from dls.numerics import autocorrelation,crosscorrelation

	from scipy.signal import savgol_filter

	fig,ax = plt.subplots(1)
	ax.plot(trace[0,:],color="red",label="trace0",alpha=0.1)
	ax.plot(trace[1,:],color="blue",label="trace1",alpha=0.1)
	ax.plot(np.array(trace[0])-np.array(trace[1]),color="green",label="diff")
	plt.legend()
	plt.show()

	N = 10
	window = np.ones((N,))/float(N)
	xcs = crosscorrelation(trace[0],trace[1])
	acs0 = savgol_filter(autocorrelation(trace[0]),51,3)
	acs1 = savgol_filter(autocorrelation(trace[1]),51,3)

	xcs = np.convolve(xcs,window,mode="valid")
	acs0 = np.convolve(acs0,window,mode="valid")
	acs1 = np.convolve(acs1,window,mode="valid")

	plt.semilogx(xcs[0:1000],'o-',label="Crosscorrelation (532nm v 633nm)")
	plt.semilogx(acs1[0:1000],'g+-',label="Autocorrelation (532nm)")
	plt.semilogx(acs0[0:1000],'rx-',label="Autocorrelation (633nm)")
	plt.xlabel("Simulation timestep ")
	plt.ylabel("Auto/Cross-correlation")
	plt.ylim(0,5)
	plt.minorticks_on()
	plt.grid(which="major",linewidth=0.2,linestyle="-")

	plt.grid(which="minor",linewidth=0.1,linestyle="--")
	plt.title("""
		Simulation steps:{5}
		Particle number: {0}
		Step sizes, linear: {1:} [m], angular(max,+/-): {2} [deg]
		Wavevelegths: {3} [nm]
		Scattering directions: {4}""".format(N_particles,step_size,angular_step_size*180/np.pi,[(1.0/(w[0]/(2*np.pi)))*1e9 for w in wavevectors],[list(w[1]) for w in wavevectors],N_steps))
	plt.legend()
	plt.show()
	


