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

input_dir = "Organisms/Colletotrichum lupini/Inputs"
organism_file = "Colletotrichum lupini.txt"
landscape_file = "homogeneous.csv" #"heterogeneous.csv"
output_dir = "Organisms/Colletotrichum lupini/Outputs"
results_file = "results.txt"
population_file = "population.txt"
save_time_steps = True
display = True
background_image = None
num_iterations = 6
num_processors = 6

#Parameters for sensitivity analysis, these are
#used for run_sensitivity, but ignored by run_gmbi
trajectory_file = "trajectory.txt"
r = 20
levels = 10
jumps = 1
perc = 0.25

#Parameters to be varied
ptv = [("ips", 0, np.inf, int, [0]),
				("d50", 0, np.inf, float, [0]),
				("d99", 0, np.inf, float, [0]),
				("rep", 0, np.inf, float, [2]),
				("survivalT", 0, 1, float, [0,1,2]),
				("K", 1, np.inf, int, [0,1])]

#PBS parameters
walltime = "12:00:00"
memory = "4GB"
pbs_script_dir = "/scratch/director309/dsavage/scripts/"
pbs_out = "/scratch/director309/dsavage/results/pbs/"
pbs_err = "/scratch/director309/dsavage/results/pbs/"
