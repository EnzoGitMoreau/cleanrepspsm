import sys
import matplotlib.pyplot as plt
from tqdm import tqdm
import subprocess
import seaborn
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
		return int(self.default_value + (self.default_value- self.endpoint) * self.it / self.nbSteps)
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
			self.value = self.default_value + (self.default_value- self.endpoint) * self.it / self.nbSteps
	def log(self):
		return self.default_value + (self.default_value- self.endpoint) * self.it / self.nbSteps
	def parse(self, v1,v2,v3):
		return float(v1),float(v2),float(v3)
PARAMETERS_TEMPLATE = [IntParameter("nbThreads", -1), IntParameter("matrixSize", -1), IntParameter("matrixRep",-1), FloatParameter("blockper", -1)]
class TestInstance():
	executable_path = "tests"
	separator =" "
	
	varying_param_log = []
	
	def __init__(self,choice, start_value, max_value, nbTests):
		self.possible_parameters = [IntParameter("nbThreads", 8), IntParameter("matrixSize", 10000), IntParameter("matrixRep",10000), FloatParameter("blockper", 0.1)]
		self.algorithms = ["ARMPL","RSB","CYTOSIM_ORIGINAL", "CYTOSIM_NEW" ,"CYTOSIM_TEST"]
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
		result = subprocess.run(command, capture_output=True, text=True)
		f = open("res/compute.out")
		f.readline()#skipping first line
		res = f.readline().strip().split(self.separator)
		i = 0
		for algo in self.algorithms:
			self.results[algo].append(int(res[i]))
			i+=1
		self.results[self.variating_parameter.name].append(int(self.variating_parameter.value))
	def createGraph(self):
		if(self.isRun):
			plt.figure(dpi=500)
			for algo in self.algorithms:
				plt.plot(self.results[self.variating_parameter.name], self.results[algo],label=algo)
				plt.title("Computational time in term of "+self.variating_parameter.name)
				plt.legend()
			plt.show()
		else:
			raise Exception("Do not run createGraph before runTests")

if __name__ == "__main__":
	if(len(sys.argv)!= 5):

		print("Usage python3 testsGenerator.py paramter_chocie, start_value, end_value, nbTests")
	else:
		b,c,d = PARAMETERS_TEMPLATE[int(sys.argv[1])].parse(sys.argv[2], sys.argv[3],sys.argv[4])
		
		instance = TestInstance(int(sys.argv[1]),b,c,d)
		instance.runTests()
		instance.createGraph()
	