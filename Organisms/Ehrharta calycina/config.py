input_dir = "Organisms/Ehrharta calycina/Inputs"
organism_file = "Ehrharta calycina.txt"
landscape_file = "Kings Park.csv"
output_dir = "Organisms/Ehrharta calycina/Outputs"
results_file = "results.txt"
population_file = "population.txt"
background_image = ""
display = True
save_time_steps = True
num_iterations = 1
num_processors = 1

#Parameters for sensitivity analysis, these are
#used for run_sensitivity, but ignored by run_gmbi
#trajectory_file = "Trajectory.trj"
trajectory_file = "small_trajectories.trj"
r = 2#Number of trajectories
levels = 10
jumps = 2
perc = 0.05
walltime = "48:00:00"
memory = "6GB"
