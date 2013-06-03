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

from population_dynamics import Population_Model
from dispersal import Dispersal_Model
from gmbiIO import population_statistics


class GMBI(object):
	def __init__(self, ps, landscape):
		self._ps = ps
		self._landscape = landscape
		self._nx, self._ny = landscape.shape
		
		#Check if the simulation should end if an individual leaves the landscape bounds
		self._break_on_leave = ps.max_distance and (
			ps.max_distance < self._nx / 2. or ps.max_distance < self._ny / 2.)
	
		#Calculate the number of individuals to seed each of the initial source locations.
		self._ips = np.floor(self._ps.ips / len(self._ps.isl))
		
		#Construct the population and dispersal sub-models
		self._pmodel = Population_Model(ps)
		self._dmodel = Dispersal_Model(ps, landscape)
		
	
	#Perform a single time step
	def step(self, N):
		#Simulate growth
		nac, nx, ny = N.shape
		for i in xrange(nx):
			for j in xrange(ny):
				if np.any(N[:, i, j] > 0):
					self._pmodel.step(N[:, i, j], self._landscape[i,j])
		
		#Simulate dispersal
		age = 0
		for life_stage in xrange(self._ps.nls):
			self._dmodel.disperse(N[age,:,:], life_stage)
			age += self._ps.tsfls[life_stage]
	
	
	#Run a single simulation
	def run(self):
                print "new run"
		#Setup the initial population array, this is a three dimensional
		#array, indexed as (k, i, j) where each element represents the number of individuals
		#in age class k, in the sub-population at landscape coordinates (i, j)
		N = np.zeros((self._ps.nac, self._nx, self._ny), dtype=np.int32)				
	
		#Seed the array with the initial populations
		for i, j in self._ps.isl:
			N[:, i, j] = self._ips
	
		#Store a summary of the initial model state in a list. In each time step, we add
		#a summary snapshot to this list and then it gets returned at the end
		Ns = [np.sum(N, 0)]
	
		#Run the simulation
		for t in range(1, self._ps.max_time_steps):
                        print "time",t
			#Perform a time step, this modifes N
			self.step(N)
			
			#Store a summary snapshot for the timestep
			Ns.append(np.sum(N, 0))
			
			#Test if the model should continue to run
			if (self._break_on_leave and not self._dmodel.all_in_bounds()) or np.sum(Ns[-1]) == 0:
				break
			
		return Ns
