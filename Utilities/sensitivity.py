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


import argparse as ap
import numpy as np
import sys as sys
import os as os


#import rpy.robjects as robj
import rpy as rpy
sense = rpy.r.library("sensitivity")
sys.path.append("./Modules")
import parameters as pms
import gmbiIO as io

#Return the minimum and maximum parameter values
#for the parameter set using a variation factor perc
def _to_min_max(ps, perc, ptv):
	mins = []
	maxs = []

	for key, pmin, pmax, ptype, indices in ptv:
		try:
			nmins = [max(v - v * perc, pmin) for i, v in enumerate(ps.__dict__[key]) if (i in indices)]
			nmaxs = [min(v + v * perc, pmax) for i, v in enumerate(ps.__dict__[key]) if (i in indices)]
		except TypeError:
			nmins = [max(ps.__dict__[key] - ps.__dict__[key] * perc, pmin)]
			nmaxs = [min(ps.__dict__[key] + ps.__dict__[key] * perc, pmax)]
			
		if ptype is int:
			nmins = [int(np.floor(v)) for v in nmins]
			nmaxs = [int(np.ceil(v)) for v in nmaxs]

		mins += nmins
		maxs += nmaxs

	return (mins, maxs)


def _generate_oat_trajectories(r, p, levels, jumps, mins, maxs):
	#Set up the R function and vectors
	oatf = rpy.r("sensitivity:::random.oat")
	levelsr = rpy.r.FloatVector([levels] * p)
	jumpsr = rpy.r.FloatVector([jumps] * p)

	#Create a list to hold the trajectories
	ts = []
	
	#Generate r trajectories consisting of p + 1 points
	for i in range(r):
		ts.append(io.Trajectory(points=[]))
		tsr = oatf(p , levelsr, jumpsr)
		for point in np.array(tsr):
			#Scale the point and add it to the trajectory
			scaled = [v * (mx - mn) + mn for v, mn, mx in zip(point, mins, maxs)]
			ts[-1].add_point(io.Point(scaled))
	
	return ts
	

def _generate_simplex_trajectories(r, p, mins, maxs, h):
	x = np.zeros((r * (p + 1), p)) 
	
	simpf = rpy.r("sensitivity:::random.simplexes")
	minsr = rpy.r.FloatVector(mins)
	maxsr = rpy.r.FloatVector(maxs)
	x[:,:] = simpf(p, r, minsr, maxsr, h)
	return x.reshape((r, p+1, p))


#Generate a set of parameters using the oat method.
def _generate_trajectories(ps, r, levels, jumps, perc, ptv, design):
  mins, maxs = _to_min_max(ps, perc, ptv)
  if design == "oat":
    ts = _generate_oat_trajectories(r, len(mins), levels, jumps, mins, maxs)
  elif design == "simplex":
    h = design_params[0]
    ts = _generate_simplex_trajectories(r, len(mins), mins, maxs, h)
          
  return ts
    

def _calculate_elementary_effects(ts, ptv, design, scale=True):
  #Convert the trajectories to two amtrices, one holding the dependent variables (x) 
  #and one holding the mean of the independent variables in the final time step across all iterations
  x, y = _trajectories_to_matrices(ts, ptv)

  #tp =  total number of points, m = number of independent
  #variables, n = number of metrics, r = number of trajectories
  tp, m = x.shape
  _, n = y.shape
  r = tp / (m + 1)

  #Get the names of the paremters that were varied
  names = [name + "_" + str(i) for name, _, _, _, inds in ptv for i in inds]

  #Convert x to an R matrix
  #xr = rpy.r.matrix(rpy.r.c(x.reshape(x.size)), 
  #	nrow=tp, ncol=m, byrow=True, dimnames=[range(1, tp+1), names])
  rpy.r.assign('xmatrix', x.reshape(x.size))
  rpy.r("xmatrix <- matrix(xmatrix, nrow=" + str(tp) + ", ncol=" + str(m) + ", byrow=T)")

  #Calculate the elementary effects
  eef = rpy.r("sensitivity:::tell.morris") #EE function
  rpy.r("design <- list('type'='" + design + "')") 
  #designr = rpy.r.list({"type": design}) #Sample design
  ees = np.zeros((r, n, m)) #EEs for each point

  #For each metric given as a result
  for i in xrange(n):
    #If all response values are equal, don't bother calculating the ees. This
    #can happen if the model runs for a set number of time steps (i.e. all time values will be equal)
    if np.all(y[:,i] == y[0,i]):
      ees[:,i,:] = 0
    else:
      #Calculate the elementary effects
      #yr = rpy.r.c(y[:,i])
      rpy.r.assign("y", y[:,i])
      #morris_list = rpy.r.list({"X": xr, "y": yr, "design": designr, "scale": scale})
      #morris = eef(morris_list)
      morris = rpy.r(
        "sensitivity:::tell.morris(list('X'=xmatrix, 'y'=y, 'design'=design, 'scale'=T))")
      ees[:,i,:] = np.array(morris["ee"])

  return ees


def _trajectories_to_matrices(ts, ptv):
	#Get the dimensions for the matrix, (r * (p + 1) x p), where r is the number
	#of trajectories, and p is the number of points in each trajectory. The
	#function dims also returns n, the number of metrics tracked, and i, the
	#number of iterations performed for each point
	r = len(ts)
	p, n, i = ts[0].dims()
	
	#Construct the trajectory and results matrices
	x = np.zeros((r * (p + 1), p))
	y = np.zeros((r * (p + 1), n))
	
	#Copy the data from the trajectory objects to the matrices
	for i, t in enumerate(ts):
		for j, point in enumerate(t.points()):
			index = i * (p + 1) + j
			
			x[index,:] = point.independents()
			
			if np.all(x[index, :] == x[index - 1,:]):
				x[index, :] += x[index, :] / 1000.
				y[index, :] = y[index - 1, :]
			else:
				y[index,:] = point.mean()

	return x, y


if __name__ == "__main__"	:
	parser = ap.ArgumentParser()
	parser.add_argument('--config', default="config.py")
	parser.add_argument('--task', required=True)
	parser.add_argument('--design', default="oat")
	args = parser.parse_args()
	
	#Load the configuration file
	config = os.path.realpath(args.config)
	execfile(args.config)

	if "generate" == args.task:
		ps = pms.parameter_set_from_file(input_dir + "/" + organism_file)
		ts = _generate_trajectories(ps, r, levels, jumps, perc, ptv, args.design)
		io.write_trajectories(ts, ptv, [], output_dir + "/" + trajectory_file)
	elif "analyse":
		ts, ptv, sentinels = io.read_all_trajectories(output_dir + "/" + results_file)
		ees = _calculate_elementary_effects(ts, ptv, "oat")
		io.plot_ees(ees, ptv, output_dir)
