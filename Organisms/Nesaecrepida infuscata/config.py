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

input_dir = "Organisms/Nesaecrepida infuscata/Inputs"
organism_file = "Nesaecrepida infuscata.txt"
landscape_file = "abstract mimosa.csv"
output_dir = "Organisms/Nesaecrepida infuscata/Outputs"
results_file = "results.txt"
population_file = "population.txt"
save_time_steps = True
display = True
background_image = None
num_iterations = 1
num_processors = 1

#Parameters for sensitivity analysis, these are
#used for run_sensitivity, but ignored by run_gmbi
trajectory_file = "trajectories.txt"
r = 2#Number of trajectories
levels = 10
jumps = 1
perc = 0.5
walltime = "12:00:00"
memory = "4GB"
