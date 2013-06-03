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

import multiprocessing as mp
import multiprocessing.pool
import subprocess as sp
import argparse as ap
import time as time
import numpy as np
import sys as sys
import os as os
import re as re
import gc as gc

#Add the GMBI modules path to the system path variable
#if ("USER" in os.environ):
 #   if(os.environ["USER"] == "mrenton" or os.environ["USER"] == "dsavage"):
#		mod_path = "/home/" + os.environ["USER"] + "/models/GMBI/Modules"
#		sys.path.append(mod_path)
#else:
sys.path.append("./Modules")

import parameters as params
import landscape as lscp
import gmbi as gmbi
import gmbiIO as io

#Some module variables that specify input and output locations
#and information for the generated PBS scripts
if "USER" in os.environ:
	pbs_header = "#!/bin/bash\n" + "#PBS -W group_list=director309\n"
	pbs_loads = "module load python-numpy\n"
	pbs_loads += "module load python-scipy\n"
	pbs_loads += "module load python-matplotlib\n"
	run_string = ("python /home/" + os.environ["USER"] + 
						"/models/GMBI/run_gmbi.py --task run_trajectories")


def run_iteration(ps, landscape):
	model = gmbi.GMBI(ps, landscape)
	return model.run()


def process_iteration(Ns, ps, landscape, config):
	io.process_iteration(Ns, ps, landscape, config)

#Run num_iteration simulations using the given parameter set and landscape.
#Each iteration is run in it's own operating system process
def run_parrallel_iterations(ps, landscape, config):
	iteration = 0

	if config.num_iterations > 1:
		means = [np.zeros(landscape.shape) for _ in xrange(ps.max_time_steps)]

	#run_iteration(ps, landscape)

	#Perform the iterations using a process pool for concurrency
	pool = mp.Pool(num_processors)
	print "running...", num_iterations, "iterations for", ps.max_time_steps,"timesteps on",num_processors,"processors"
	results = [pool.apply_async(run_iteration, args=[ps, landscape]) for _ in xrange(num_iterations)]
	pool.close()
	
	#Process iterations as they complete
	while len(results) > 0:
		completed = [i for i, r in enumerate(results) if r.ready()]
		for i in reversed(completed): #reversed so that indices aren't invalidated as we pop
			print "processing iteration " + str(iteration + 1)
			#Get the result from the list and save the iteration to file
			Ns = results.pop(i).get(None)
			config.ext = "/iteration " + str(iteration + 1)
			process_iteration(Ns, ps, landscape, config)
			iteration += 1
			
			if config.num_iterations > 1:
				#Add the population for each time step in the iteration to the total
				for t, N in enumerate(Ns):
					means[t] += N
		
		time.sleep(15)
	
	pool.close()
	pool.join()
	#run_iteration(ps, landscape)

	if config.num_iterations > 1:
		for N in means:
			N /= config.num_iterations
		config.ext = "/means"
		io.process_iteration(means, ps, landscape, config)


#Run num_iteration simulations using the given parameter set and landscape
def run_iterations_for_point(p, ps, landscape, num_iterations):
	model = gmbi.GMBI(ps, landscape)
	point = io.Point(independents=p.independents()) #no need to copy as we don't use p

	for i in xrange(num_iterations):
		res = model.run()
		point.add_iteration()

		for t, N in enumerate(res):
			point.add_time_step([t] + io.population_statistics(ps, landscape, N))
	
	return point


def run_trajectory(t, ps, landscape, ptv, num_iterations, num_processors):
	#Get the points in the trajectory
	points = t.points()

	#Determine the index of each unique point (sometimes points are equal due to rounding)
	uinds = [i for i, p in enumerate(points) if i == 0 or not p.equals(points[i-1])]

	#Create a process pool, using as many processors as are available, or
	#are required to allow each point to run concurrently
	pool = mp.Pool(processes=min(num_processors, len(points)))
	
	results = []
	for i in uinds:
		#Modify the parameter set to match the current point
		psm = ps.copy()
		psm.modify_for_point(points[i], ptv)
		psm.convert_to_age_classes()

		#Launch a process to run the simulation(s) for the point. This modifies the point in place
		args = [points[i], psm, landscape, num_iterations, num_processors]
		results.append(pool.apply_async(run_iterations_for_point, args))

	pool.close()
	pool.join()

	#Merge the unique and non-unique points back together
	for i, r in zip(uinds, results):
		points[i] = r.get(None)

	#Return a new trajectory containing the results for each point
	return io.Trajectory(points=points)


def run_trajectories_remotely(ts, processors, walltime, memory, output_dir, trj_file, config_file):
	#List to hold the PBS job ids
	ids = []
	
	#Go through each trajectory
	for i, t in enumerate(ts):
		#Create the directory to hold all the outputs
		ext = str(i + 1)
		trd = output_dir + "/Trajectory" + ext
		os.makedirs(trd)

		#Determine the names for the pbs script and the individual trajectory file
		pbsfile = trd + "/runTrajectory" + str(i + 1) + ".pbs"
		tfile = trd + "/" + trj_file + "." + str(i + 1) 
		
		#Create a file containing the current trajectory, located in the newly created directory
		io.write_trajectories([t], ptv, [], tfile)
		
		#Create a pbs script to run the trajectory
		write_pbs(pbsfile, processors, walltime, memory, config_file, ext)
		
		#Submit the script
		ids.append(sp.check_output(["qsub", "-o", trd, "-e", trd, pbsfile]).strip())

	#Wait for each of the processes to finish
	while ids != []:
		queue = sp.check_output(["qstat"]).strip()
		print "polling: " + str(len(ids)) + " trajectories remaining"
		for id_string in ids:
			if not re.search(id_string, queue):
				ids.remove(id_string)

		time.sleep(60)


def write_pbs(filename, processors, walltime, memory, config_file, ext):
		try:
			f = open(filename, "w")
			
			#Write the header and required resources
			f.write(pbs_header)
			f.write("#PBS -l walltime=" + str( walltime) + ",mem=" + 
				memory + ",nodes=1:ppn=" + str(processors) + "\n")
			
			#Modules to load
			f.write(pbs_loads)
			
			#Execution string
			f.write(run_string + " --config " + config_file + " --ext " + ext + "\n")
			
			f.close()
		except IOError as e:
			print "Error writing PBS file " + filename + "\n"
			print e
			sys.exit()
			

if __name__ == "__main__":
	parser = ap.ArgumentParser()
	parser.add_argument('--config', default="./config.py")
	parser.add_argument('--remote', action="store_true")
	parser.add_argument('--task', default="run")
	parser.add_argument('--ext', default="")
	parser.add_argument('--trajectory', default=-1)
	args = parser.parse_args()
	
	#Load the configuration file
	config_file = os.path.realpath(args.config)
	execfile(config_file)
	config = params.Config(num_iterations, num_processors, input_dir, output_dir,
								display, save_time_steps, background_image=background_image)
		
	#Construct the parameter set from the organism file specified in config
	ps = params.parameter_set_from_file(input_dir + "/" + organism_file)
	
	#Construct the landscape from the landscape file specified in config
	landscape = lscp.landscape_from_file(input_dir + "/" + landscape_file)
	ps.landscape_size = landscape.shape

	#If a single run is being performed
	if args.task == "run":
		ps.convert_to_age_classes()
		Ns = run_parrallel_iterations(ps, landscape, config)
	elif args.task == "run_trajectories":
		#Read in the trajectories from file, if ext is set, the original trajectory file has been split into multiple
		#files and each trajectory has been assigned a working sub-directory under the output_directory
		if args.ext == "":
			ts, ptv, _ = io.read_all_trajectories(input_dir + "/" + trajectory_file)
		else:	
			output_dir += "/Trajectory" + args.ext
			ts, ptv, _ = io.read_all_trajectories(output_dir + "/" + trajectory_file + "." + args.ext)
		
		if args.remote:
			run_trajectories_remotely(
				ts, num_processors, walltime, memory, output_dir, trajectory_file, config_file)
			io.cat_trajectories(len(ts), ptv, ps.sentinels, output_dir)
		else:
			#If a specific trajectory is to be run
			if args.trajectory > 0:
				ts = [run_trajectory(ts[int(args.trajectory)], ps, landscape, ptv, num_iterations, pool)]
			else:
				ts = [run_trajectory(t, ps, landscape, ptv, num_iterations, pool)  for t in ts]
			io.write_trajectories(ts, ptv, ps.sentinels, output_dir + "/results.txt")
