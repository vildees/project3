from scitools.std import *

class graph:

	def __init__(self, filename):
		self.filename = filename
		
		self.extract_data()
	

	def extract_data(self):
			 
		infile = open(self.filename, "r")
		lines = infile.readlines()

		t_Earth = []
		r_Earth = []
		
		for line in lines:
			values = line.split()
			t_Earth.append(float(values[0]))
			r_Earth.append(float(values[1]))

		self.t_Earth = array(t_Earth)
		self.r_Earth = array(r_Earth)
		
		self.plot()

	def plot(self):
		
		plot(self.t_Earth, self.r_Earth)
		axis([0, 0.04, 0, 20])
		title("Stability of algorithm for different time steps")
		xlabel("dt")
		ylabel("Radius")
		raw_input("")
		
filename = "timestep.dat"
h = graph(filename)
