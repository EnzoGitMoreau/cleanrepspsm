import sys
import subprocess
try:
	import matplotlib.pyplot as plt
except:
	print("Requirement 'MatplotLib' not satisfied, trying to install automatically...")
	subprocess.run(["pip"]+["install","matplotlib"], capture_output = True)
	try:
		import matplotlib.pyplot as plt 
		print("Requirement 'Matplotlib' installed sucessfully")
	except:
		print("Could not install automatically, installing pip...")
		subprocess.run(["python","m","ensurepip","--upgrade"])
		print("pip installed, re-trying to install 'Matplotlib")
		subprocess.run(["pip"]+["install","matplotlib"], capture_output = True)
		try:
			import matplotlib.pyplot as plt 
			print("Requirement 'Matplotlib' installed sucessfully")
		except:
			print("Automatical installation failed")
try:
	from tqdm import tqdm
except:
	print("Requirement 'tqdm' not satisfied, trying to install automatically...")
	
	try:
		subprocess.run(["pip"]+["install","tqdm"])
		from tqdm import tqdm
		print("Requirement 'tqdm' installed sucessfully")
	except:
		print("Could not install automatically, installing pip...")
		subprocess.run(["python","-m","ensurepip","--upgrade"])
		print("pip installed, re-trying to install 'tqdm")
		subprocess.run(["pip"]+["install","tqdm"])
		try:
			from tqdm import tqdm
			print("Requirement 'tqdm' installed sucessfully")
		except:
			print("Automatical installation failed")
try: 
	import seaborn
except:
	print("Requirement 'seaborn' not satisfied, trying to install automatically...")
	subprocess.run(["pip"]+["install","seaborn"], capture_output = True)
	try:
		import seaborn
		print("Requirement 'seaborn' installed sucessfully")
	except:
		print("Could not install automatically, installing pip...")
		subprocess.run(["python","m","ensurepip","--upgrade"])
		print("pip installed, re-trying to install 'seaborn")
		subprocess.run(["pip"]+["install","seaborn"], capture_output = True)

		try:
			import seaborn
			print("Requirement 'seaborn' sucessfully installed")
		except:
			print("Automatical installation failed")
#BASE VALUES
DEFAULT_MATRIXSIZE = 20000
DEFAULT_NBMULTIPLICATION = 1000
DEFAULT_BLOCKPERCENTAGE = 0.2
DEFAULT_NBTHREADS = 4


#
#if no Internet please comment the line 7 
plt.style.use('https://github.com/dhaitz/matplotlib-stylesheets/raw/master/pitayasmoothie-light.mplstyle')


#plt.style.use('seaborn')
class Parameter():
	def __init__(self,name,default_value):
		self.value  = default_value
		self.name = name
		
	def get_default_value(self):
		return self.default_value
	
	def add_endpoint(self, start, value, nbSteps):
		self.default_value = start
		self.value = self.default_value
		self.endpoint = value 
		self.nbSteps = nbSteps
	def __str__(x):
		return x.name + ": "+str(x.value) 
	def __repr__(self):
		return self.__str__()
class IntParameter(Parameter):
	def __init__(self, name, default_value):
		self.name = name 
		self.default_value = default_value
		self.endpoint = None 
		self.it = 1
		self.value = default_value
	
	def next(self):
		if(self.endpoint == None):
			return self.default_value
		else: 
			self.it+=1
			self.value = int(self.default_value - (self.default_value- self.endpoint) * self.it / self.nbSteps)
	def log(self):
		return int(self.default_value - (self.default_value- self.endpoint) * self.it / self.nbSteps)
	def parse(self, v1,v2,v3):
		return int(v1),int(v2),int(v3)
class FloatParameter(Parameter):
	def __init__(self, name, default_value):
		self.name = name 
		self.default_value = default_value
		self.endpoint = None 
		self.it = 1 
		self.value = default_value
	def next(self):
		if(self.endpoint == None):
			return self.default_value
		else: 
			self.it+=1
			self.value = self.default_value - (self.default_value- self.endpoint) * self.it / self.nbSteps
	def log(self):
		return self.default_value - (self.default_value- self.endpoint) * self.it / self.nbSteps
	def parse(self, v1,v2,v3):
		return float(v1),float(v2),int(v3)
PARAMETERS_TEMPLATE = [IntParameter("nbThreads", -1), IntParameter("matrixSize", -1), IntParameter("matrixRep",-1), FloatParameter("blockper", -1)]
class TestInstance():
	executable_path = "tests"
	separator =" "
	
	varying_param_log = []
	
	def __init__(self,choice, start_value, max_value, nbTests):
		print(nbTests)
		self.possible_parameters = [IntParameter("nbThreads", DEFAULT_NBTHREADS), IntParameter("matrixSize", DEFAULT_MATRIXSIZE), IntParameter("matrixRep", DEFAULT_NBMULTIPLICATION), FloatParameter("blockper", DEFAULT_BLOCKPERCENTAGE)]
		self.algorithms = ["ARMPL","RSB","CYTOSIM_ORIGINAL","MatrixSymmetric", "CYTOSIM_NEW" ,"CYTOSIM_TEST"]
		self.choice = choice
		self.max_value = max_value
		self.start_value = start_value
		self.nbTests = nbTests
		self.variating_parameter = self.possible_parameters[choice]
		self.isRun = False
		self.variating_parameter.add_endpoint(self.start_value,max_value, nbTests)#adding trajectory 

		self.results = {}
		for algo in self.algorithms:
			self.results[algo] = []
		self.results[self.variating_parameter.name] = []
	def runTests(self):
		
		for i in tqdm(range(self.nbTests)):
			self.execWith(self.possible_parameters)
			for param in self.possible_parameters:
				param.next()
			self.varying_param_log.append(self.variating_parameter.log())
		
		self.isRun = True
	def execWith(self, parameters):
		print("Running with : "+str(parameters))
		args = [str(param.value) for param in self.possible_parameters]
		command = [self.executable_path] + args
		result = subprocess.run(command)
		f = open("res/compute.out")
		f.readline()#skipping first line
		res = f.readline().strip().split(self.separator)
		i = 0
		for algo in self.algorithms:
			self.results[algo].append(int(res[i]))
			i+=1
		self.results[self.variating_parameter.name].append(self.variating_parameter.log())
	def createGraph(self):
		if(self.isRun):
			plt.figure(dpi=100)
			for algo in self.algorithms:
				if(self.results[algo][0] != -1):
					plt.plot(self.results[self.variating_parameter.name], self.results[algo],label=algo)
					plt.title("Computational time in term of "+self.variating_parameter.name)
					plt.legend()
			plt.xlabel(self.variating_parameter.name)
			plt.ylabel("Computational time")
			plt.show()
			plt.savefig("lastgraph.png")
		else:
			raise Exception("Do not run createGraph before runTests")

if __name__ == "__main__":
	if(len(sys.argv)!= 5):

		print("Usage python3 testsGenerator.py parameter_choice, start_value, end_value, nbTests")
	else:
		b,c,d = PARAMETERS_TEMPLATE[int(sys.argv[1])].parse(sys.argv[2], sys.argv[3],sys.argv[4])
		
		instance = TestInstance(int(sys.argv[1]),b,c,d)
		instance.runTests()
		instance.createGraph()
	