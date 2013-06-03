#A General Model of Biological Invasion (GMBI)

#Copyright (C) 2012 David Savage and Michael Renton

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.


import re as re
import numpy as np
import copy as cp
import random as random

#class ParameterVariation(object):
#	def __init__(self, pname, pmin, pmax, indices, mods=[0] * len(indices)):
#		self._name = pname
#		self._min = pmin
#		self._max = pmax
#		self._indices = indices
#		self._modified_values = None
	
#	def __str__(self):
#		return "(" + str(self._name + ", " + self._min + ", " + self._max + ", " + self._indices + ")"
		
#	def modify(self, p):
#		p[indices] = self._modified_values
#		return p

class Config(object):
	def __init__(self, num_iterations, num_processors, input_dir, output_dir, display, save_time_steps, ext="", background_image=""):
		self.num_iterations = num_iterations
		self.num_processors = num_processors
		self.input_dir = input_dir
		self.output_dir = output_dir
		self.display = display
		self.save_time_steps = save_time_steps
		self.ext = ext
		self.background_image = background_image


class ParameterSet(object):
	def __init__(self):
		self.desc = None #Organism dsecription, e.g. "typical weed", "Phytophthora infestans"
		self.time_step = None	#The unit used to measure time steps, e.g. "day", "week"
		self.spatial_units = None #The unit used to measure cell dimensions, e.g. "km", "m"
		self.cell_size = None #The size of the cell expressed in terms of the spatial units
		self.max_time_steps = None #Maximum number of time steps
		self.max_distance = None #Maximum distance travelled before simulation stops
		self.landscape_size = None #The maximum x and y coordinates for the landscape
		self.psuit = None #Proportion of the landscape covered by each patch type
		self.suitability = None #Suitability of each patch type
		self.avg_cluster = None #Avg. size of a habitat cluster (i.e the number of adjacent suitable cells)
		self.landscape_type = None #The type of habitat clusters in the landscape (can be uniform, blocked or gaussian)
		self.nls = None #Number of life stages	
		self.tsfls = None #Time steps for each life stage, each life stage will be expanded into n age-classes, where tsfls = [n_1, ..., n_nls]
		self.K = None #Carrying capacity for each life stage, alternatively, an overall carrying capacity can be specified
		self.nac = None #Number of age classes, derived from nls and tsfls
		self.ncomp  = None #Number of competing age classes
		self.comp = None #Inter-age-class competition matrix, derived from nls and K
		self.ips = None #Initial population sizes, expanded using tsfls
		self.isl = None #Initial source location
		self.sentinels = None #Cells whose statistics should be output for each time step
		self.d50 = None #Median dispersal distance for each life stage, expanded using dls and dac
		self.d99 = None #99th precentile of dispersal distance for each life stage, expanded using dls and dac
		self.davg = None #Avg. dispersal direction (von mises), expanded using dls and dac
		self.dconc = None #Concentration of dispersal directions (von mises), expanded using dls and dac
		self.arad = None #Radius for active dispersal, expanded using dls and dac
		self.acrit = None #Patch selection criteria for active dispersal, ignored when arad == 0, expanded using dls and dac
		self.afa = None #Indicates whether the organism can actively move and select suitable habitat in each age class, derived from arad
		self.afls = None #Indicates whether the organism can actively move and select suitable habitat in each life stage, derived from arad
		self.allee = None #Allee threshold for each life stage, when the population size less than allee, reproduction does not occur
		self.rep = None #Reproduction rate for each life stage, expanded using tsfls
		self.survivalR = None #Probability that individual survives, but remains in the same age class
		self.survivalT = None #Probability that individual survives, and then transititons to the enxt age class
		self.hibernation_start = None #End time time step for life cycle
		self.hibernation_end = None #Start time step for life cycle
		self.tsfc = None #Time steps per cycle

	
	def copy(self):
		psc = ParameterSet()
		psc.__dict__ = cp.deepcopy(self.__dict__)
		return psc

		
	def to_list(self, keys):
		ps_list = []
		for k in keys:
			v = self.__dict__[k]
			if isinstance(v, list):
				ps_list += [vi for vi in v]
			elif not isinstance(v, str):
				ps_list += [v]
		return ps_list

		
	def header(ps, keys):
		head = []
		for k in keys:
			v = self.__dict__[k]
			if isinstance(v, list):
				head += [k + str(i + 1) for i in range(len(v))]
			elif not isinstance(v, str):
				head += [k]
		return head


	def modify_for_point(self, point, ptv):
		k = 0

		for key, pmin, pmax, ptype, indices in ptv:
			try:
				#If the parameter is vector vsalued
				for i, ind in enumerate(indices):
					self.__dict__[key][ind] = point.independents()[k + i]
			#If the parameter is not vector valued
			except TypeError:
				self.__dict__[key] = point.independents()[k]

			#Shift to the index of the next parameter
			k = k + len(indices)
				
			#if ptype == int:
				#self.__dict__[key] = np.round(self.__dict__[key]).astype(np.int)
				#	self.__dict__[key] = [int(round(v)) for v in self.__dict__[key]]
				#else:
				#	self.__dict__[key] = int(np.round(self.__dict__[key]))

			#if length > 1:
				#print self.__dict__[key]
				#self.__dict__[key] = self.__dict__[key].tolist()
			#else:
				#print self.__dict__[key]
				#self.__dict__[key] = self.__dict__[key][0]


	def convert_to_age_classes(self):
		self.nac = int(sum(self.tsfls))
	
		#Convert the initial population sizes to age-classes. Set the initial
		#population of the first age-class in each life-stage to that originally
		#specified for the life-stage
		self.ips = sum([[ips] + [0] * (ts - 1) for ips, ts in zip(self.ips, self.tsfls)], [])
		self.ips = np.array(self.ips)
	
		#Convert the proportion of surviving individuals that remain in the current age-class	
		self.survivalR = sum([[0.] * (ts - 1) + [s] for s, ts in zip(self.survivalR, self.tsfls)], [])
		self.survivalR = np.array(self.survivalR)
	
		#Convert the proportion of surviving individuals that transition to the next age-class
		self.survivalT = sum([[s] + [1] * (ts - 1) for s, ts in zip(self.survivalT, self.tsfls)], [])
		self.survivalT.pop() #The final age-class can't transition, so we drop the last term
		self.survivalT = np.array(self.survivalT)
		
		#Convert K
		self.totalK = np.sum(self.K)
		self.ncomp = len(self.tsfls) - len(self.K)
		self.K = sum([[k] * ts for k, ts in zip(self.K, self.tsfls)], [])
		self.K = self.K + [self.K[-1]] * (self.nac - len(self.K))
		self.K = np.array(self.K)
	
		#Convert the reproductive rates
		self.rep = sum([[d] * ts for d, ts in zip(self.rep, self.tsfls)], [])
		self.rep = np.array(self.rep)
	
		#Determine the age-classes where active dispersal takes place
		self.afa = np.array(sum([[v] * ts for v, ts in zip(self.afls, self.tsfls)], []))
		
		#Determine the age classes where dispersal takes place
		self.dfls = np.array([i for i, (d, r) in enumerate(zip(self.d50, self.arad)) if d > 0 or r > 0])


#Factory constructor
def parameter_set_from_file(filename):
	ps = ParameterSet()
	
	try:
		for line in open(filename):
			expr = re.sub("\s*#.*$", "", line).strip()
			if expr:
				(pname, val) = re.split("\s*=\s*", expr)
				if re.search("\[.*\]", val):
					ps.__dict__[pname] = np.array(eval(val))
				else:
					ps.__dict__[pname] = eval(val)
		ps.afls = ps.arad > 0
	
		#Tsfls must be an int, or we can't multiply sequences with it's elements
		ps.tsfls = ps.tsfls.astype(np.int)
		
		#Cell size must be a float or when we divide by
		#it later, bad things happen
		ps.cell_size = float(ps.cell_size)
	except IOError as e:
		print "Problem reading parameters from file " + filename
		print e
		exit()
	return ps


#Multiplies the array a by the scalar K when K is potentially infinity
def _mult(a, K):
	return np.array([e * K if e != 0 else 0 for e in a])


def _random_choice(x):
	return random.choice(x)

	


	