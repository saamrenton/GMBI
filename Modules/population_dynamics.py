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

class Population_Model(object):
	def __init__(self, ps):
		self.ps = ps
		self.comp = self._competition_matrix(ps.tsfls, ps.K)
		(self._mature_inds,) = np.nonzero(self.ps.rep > 0)
		self._Na = np.zeros(self._mature_inds.size)
		
	#Step the given population through one time step based on the given habitat suitability
	def step(self, N, suitability):
		#Calculate the expected number of individuals given unconstrained growth 
		self._step_unconstrained(N, suitability)

		#Sample a Poisson distribution to get the stochastic
		#realisation of the unconstrained growth
		#self._apply_poisson(N)

		#Apply constraints to population growth
		self._constrain(N, suitability)

	
	#Unconstrained growth and maturation
	def _step_unconstrained(self, N, suitability):
		#This function mimics the behaviour of an Usher projection matrix.
		#The usher matrix contains a top row, a diagonal, and an off-diagonal, with
		#all other values set to zero. In this function the matrix multiplication
		#is performed by only considering the non-zero elements and performing the
		#relevant calculations to obtain the next population
		
		#Calculate the total number of mature adults (i.e those whose
		#reproduction rate is greater than zero
		total_matures = np.sum(N[self._mature_inds])

		#Calculate the total number of offspring. This is equivalent to
		#multiplying the top matrix row by the population column vector.
		#Since we're using an allee threshold though, we reduce the
		#number of individuals that actually reproduce
		offspring = 0
		if total_matures > self.ps.allee:
			self._Na[:] = (N[self._mature_inds] / total_matures) * (total_matures - self.ps.allee)
			offspring = sum(self.ps.rep[self._mature_inds] * self._Na)
		
		#Maturation. If we performed the matrix multiplication, each element
		#of the resulting population vector, except for the first one, would 
		#be of the form T_i-1 * N_i-1 + R_i * N_i. If we build up the vector
		#in reverse, we can re-use the same storage.
		#The first element (i = 0) is given by R_0 * N_0 + offspring
		for i in reversed(range(1, N.size)):
			N[i] = np.random.poisson(self.ps.survivalT[i-1] * N[i-1] + self.ps.survivalR[i] * N[i])
			#N[i] = self.ps.survivalT[i-1] * N[i-1] + self.ps.survivalR[i] * N[i]
			
		N[0] = np.random.poisson(N[0] * self.ps.survivalR[0] + offspring)
		#N[0] = N[0] * self.ps.survivalR[0] + offspring
	
	
	#Constrain a population according tho the carrying capacity and any inter-cohort competition
	def _constrain(self, N, suitability):
		#Modify the carrying capacity based on the habitat suitability
		K = suitability * self.ps.K

		#Cap the number of individuals in each age-class and 
		#determine the effective population size for competing age-classes
		capped = np.minimum(N, K)
		N[:] = np.dot(self.comp, capped)
		
		#Calculate the proportion that each age-class
		#contributes to the effective population size
		p = [c / n if n > 0 else 0. for c, n in zip(capped, N)]
		
		#Calculate the proportion of the effective population size that can remain
		p *= np.minimum([k / n if n > 0 else 0 for k, n in zip(K, N)], 1)
		
		#Adjust the population
		N[:] = np.round(N * p)
	
	
	#Setup the competition matrix based on the carrying capacity for each life stage.
	#If there are less carrying capacities specified than life stages, we assume
	#that all remaining life stages compete
	def _competition_matrix(self, tsfls, K):
		comp = np.ones((sum(tsfls), sum(tsfls)))
		ncomp = len(tsfls) - len(K) #Number of competing life-stages
		ind = 0
		
		for i, ts in enumerate(tsfls):
			for j in range(ts):
				if i < len(K) - 1:
					#Non-competing classes are marked zero
					comp[ind+j, ind+ts:] = 0
					comp[ind+ts:, ind+j] = 0
			ind += ts
	
		return comp
	
	
	def _apply_poisson(self, N):
		N[:] = np.random.poisson(N, N.shape)
