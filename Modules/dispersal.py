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


import numpy as np
import ctypes as ct
import random as random
import scipy.integrate as quad
import scipy.special as special

class Dispersal_Model(object):
	def __init__(self, ps, landscape):
		if not all([0 < d50 < d99 or (d50 == 0 and d99 == 0) for d50, d99 in zip(ps.d50, ps.d99)]):
			print "Error with parameters d50 and d99, given " + str(ps.d50)  + " and " + str(ps.d99) + \
				" but d50 must be greater than zero, and d99 must be greater than d50\n"
			exit(0)
		
		self._ps = ps
		self._landscape = landscape
		self._img = np.zeros(landscape.shape, dtype=np.int32)
		self._non_suitable = np.nonzero(1 - np.ceil(landscape))
		
		#Load the C dynamic library and initilaise the static generator
		#self._dlib = ct.cdll.LoadLibrary("Modules/CDispersal/dispersal.dylib")
		#self._dlib.init_generator(0x8f7011ee, 0xfc78ff1f, 0x3793fdff)

		#Setup a pointer to the C dispersal function. We need to specify how python types
		#should be converted to C types. We also specify the return type.
		#self._dispersalf = self._dlib.disperse
		#self._dispersalf.argtypes = [ct.c_int, ct.c_int, np.ctypeslib.ndpointer(ct.c_int),
		#						np.ctypeslib.ndpointer(ct.c_int), ct.c_float, ct.c_float, ct.c_float, ct.c_float, ct.c_float]
		#self._dispersalf.restype = None
		
		#Flag indicating if dispersal occurred outside the landscape bounds
		self._all_in_bounds = True
					
		#Calculate the Weibull shape and scale parameters
		self._shape = [_shape(d50 / ps.cell_size, d99 / ps.cell_size) for d50, d99 in zip(ps.d50, ps.d99)]
		self._scale = [_scale(d50 /  ps.cell_size, s) for d50, s in zip(ps.d50, self._shape)]
		
		#Find the probability of dispersal into the circle of radius 0.5 * cell_size,
		#which gives the probability that an individual will disperse into the same cell it is dispersing from
		self._pself = [_weibull_cdf(0.5, sh, sc) for sh, sc in zip(self._shape, self._scale)]
		
		#Determine whether dispersal is active or passive and set the dispersal
		#functions accordingly. If dispersal is active, construct a dictionary
		#that stores the list of neighbouring for each coordinate pair (i, j),
		#and the criteria function for assessing neighbouring cells

		self._neighbours = []
		self._acrit = []
		
		for arad, acrit in zip(ps.arad / ps.cell_size, ps.acrit):
			if arad > 0 and acrit: #If dispersal is active
				self._neighbours.append(_nearest_neighbours(arad, landscape))
				self._acrit.append(_criteria_function(acrit))
			else: #If dispersal is passive or does not occur for the life-stage
				self._neighbours.append(None)
				self._acrit.append(None)
	
	
	#Indicate whether all individuals for the previous dispersal event
	#dispersed within the simulation bounds
	def all_in_bounds(self):
		return self._all_in_bounds
	
	
	#Perform dispersal
	def disperse(self, N, life_stage):
		if self._ps.afls[life_stage]:
			#Perform passive dispersal then active dispersal
			self._disperse_passive(N, life_stage)
			self._disperse_active(N, life_stage)
		else:
				self._disperse_passive(N, life_stage)		
			
		#Remove any individuals in totally unsuitable population
		N[self._non_suitable] = 0


	#Perform passive dispersal
	def _disperse_passive(self, N, life_stage):
		#If the dispersal kernel is setup for the given life stage
		if self._shape[life_stage] > 0 and self._scale[life_stage] > 0:
			#Zero the intermediate immigrant matrix
			self._img[:,:] = 0
			ni, nj = self._img.shape

			disperse_C = False  ### this is where you switch to using the C code dispersal
			if disperse_C:
				self._dispersalf(ni, nj, N, self._img, self._pself[life_stage], self._scale[life_stage],
								self._shape[life_stage], self._ps.davg[life_stage], self._ps.dconc[life_stage])
			else:
				#Go through each population and passively disperse each individual,
				#mutating _img as we go
				for si in xrange(ni):
					for sj in xrange(nj):
						if N[si, sj] > 0:
							#Perform dispersal for those individuals that end up back in the source cell
							source_img = np.random.binomial(N[si, sj], self._pself[life_stage])
							self._img[si, sj] += source_img
							N[si, sj] -= source_img
	
							for k in xrange(N[si, sj]):
								try:
									if self._ps.dconc[life_stage] > 0:
										di, dj = _ptc(_rweibull(self._shape[life_stage], self._scale[life_stage], 1),
											np.random.vonmises(self._ps.davg[life_stage], self._ps.dconc[life_stage])) 
									else:
										di, dj = _ptc(_rweibull(self._shape[life_stage], self._scale[life_stage], 1),
											np.random.uniform(0, 2 * np.pi)) 
									self._img[si + di, sj + dj] += 1
								except IndexError:
									self._all_in_bounds = False

			N[:] = self._img[:]


	#Perfrom active dispersal
	def _disperse_active(self, N, life_stage):
		self._img[:] = 0
		
		#Go through each non-zero population, and have each individual
		#select a suitable neighbour, and then move to that cell
		for i, j in zip(*np.nonzero(N)):
			sns = [self._most_suitable_neighbour(i, j, N, life_stage) for k in range(int(N[i,j]))]
			
			for ind in sns:
				self._img[ind] += 1
			
			#self._img[self._most_suitable_neighbour(i, j, N, life_stage)] = N[i,j]

		N[:] = self._img

	
	#Get the most suitable neighbour based on the active radius and criteria function
	def _most_suitable_neighbour(self, i, j, N, life_stage):
		#Get the list of neighbours from the data structure
		sns = self._neighbours[life_stage][i][j]

		#If no suitable neighbours, return the given location
		if 0 == len(sns):
			return (i, j)
		#If only one suitable neighbour exists
		elif 1 == len(sns):
			return sns[0]
		#If there are multiple suitable neighbours
		else:
			#Get the neighbouring populations and suitability values
			neighbours = [N[n] for n in sns]
			suitability = [self._landscape[n] for n in sns]
			
			#Select one of the neighbouring populations that 
			#fits the specified criteria
			return sns[self._acrit[life_stage](neighbours, suitability)]


def max_suitability(n, s):
	best = max(s)
	return random.choice([i for i, v in enumerate(s) if v == best])

def max_population(n, s):
	best = max(n)
	return random.choice([i for i, v in enumerate(n) if v == best])
	
def min_suitability(n, s):
	best = min(s)
	return random.choice([i for i, v in enumerate(s) if v == best])
	
def min_population(n, s):
	best = min(n)
	return random.choice([i for i, v in enumerate(n) if v == best])
	
def random_neighbour(n, s):
	return random.choice(range(len(n)))


#Determine the criteria function for selecting
#a neighbouring cell during active dispersal
def _criteria_function(acrit):
	return eval(acrit)


#Weibull shape parameter
def _shape(d50, d99):
	if 0 < d50 < d99:
		return np.log(np.log(100.) / np.log(2.)) / np.log(d99 / d50)
	else:
		return 0.

#Weibull scale parameter
def _scale(d50, shape):
	if shape > 0:
		return d50 * np.log(2.)**(-1. / shape)
	else:
		return 0.

#Weibull cdf function
def _weibull_cdf(x, shape, scale):
	if shape > 0 and scale > 0:
		return 1. - np.exp(-(x / scale)**shape)
	else:
		return 0.


#Convert from polar to cartesian (indexed) coordinates
def _ptc(r, t):
	return (int(r * np.cos(t) + 0.5), int(r * np.sin(t) + 0.5))


def _rweibull(shape, scale, n):
	#0.5 is the radius of the largest circle contained in the source grid cell.
	#Since we've accounted for dispersal into this circle, we only draw
	#random variates from outside the circle, thus our uniform variates
	#are drawn from the interval [x,1] where x = w-cdf(0.5, shape, scale)
	us = np.random.uniform(1 - np.exp(-((0.5 / scale)**shape)), 1, n)
	return scale * ((-np.log(1 - us))**(1. / shape))

#Determine the suitable neighbours for each cell in the landscape grid.
#Return a list of lists of tuples, representing the neighbouring suitable
#habitats for each patch in the landscape. This data structure is okay for
#random access, as python lists are implemented as arrays.
def _nearest_neighbours(radius, landscape):
	dxs = range(int(-radius), int(radius + 1))
	
	nx, ny = landscape.shape

	def f(i, j):
		return 0 <= i < nx and 0 <= j < ny and landscape[i,j] > 0
	
	#Generate the set of neighbouring patches that are suitable and within the domain
	return [[[(i + dx, j + dy) for dx in dxs for dy in dxs 
			if f(i + dx, j + dy)] for j in range(ny)] for i in range(nx)]
