import numpy as np
import matplotlib.pyplot as plt 


#simulate a single chucnk of the data
def simulate_chunk(linear_stepsize,angular_stepsize,N_steps,init=np.zeros((1,6))):
	linear_steps = np.random.uniform(-linear_stepsize,linear_stepsize,(N_steps-1,3))
	angular_steps = np.random.uniform(-angular_stepsize,angular_stepsize,(N_steps-1,3))
	steps = np.concatenate([linear_steps,angular_steps],axis=1)
	steps = np.concatenate([init,steps],axis=0)
	m = np.tril(np.ones((N_steps,N_steps)))
	outp = np.dot(m,steps)
	return outp[1:,:] #return everything but the first value 

def generate_steps(linear_stepsize,angular_stepsize,N_steps,init=np.zeros((1,6))):
	chunk_size = 10000
		
	if N_steps > chunk_size:
		chunks = N_steps/chunk_size

		outputs = []
		for i in range(chunks):
			if i > 0:
				init = outputs[-1][-1,:]
				init = init.reshape((1,6))
				print init.shape
			outputs.append(simulate_chunk(linear_stepsize,angular_stepsize,chunk_size,init))

		if float(N_steps)/chunk_size > N_steps/chunk_size:
			init = outputs[-1][-1,:]
			outputs.append(simulate_chunk(linear_stepsize,angular_stepsize,N_steps-chunk_size*chunks,init))

		outputs = [np.zeros((1,6))] + outputs
		outp = np.concatenate(outputs,axis=0)
	else:

		outp = simulate_chunk(linear_stepsize,angular_stepsize,N_steps,init=np.zeros((1,6)))
		outp = np.concatenate([init,outp],axis=0)
	return outp 

#write array to filepath
def writeout(trace_array,filepath):
	print "Writing trace: {}".format(filepath)
	np.save(filepath,trace_array)

#load array from filepath	
def loadin(filepath):
	print "Loading trace: {}".format(filepath)
	return np.load(filepath)

if __name__ == "__main__":
	# outp = generate_steps(1,1,100000)
	# writeout(outp,"traces/trace0.npy")
	outp = loadin("traces/trace0.npy")
	plt.plot(outp[:,1])
	plt.show()