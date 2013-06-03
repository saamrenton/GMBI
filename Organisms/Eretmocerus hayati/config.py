input_dir = "Organisms/Eretmocerus hayati/Inputs"
organism_file = "Eretmocerus hayati.txt"
landscape_file = "Carnarvon.csv"
output_dir = "Organisms/Eretmocerus hayati/Outputs"
results_file = "results.txt"
population_file = "population.txt"
save_time_steps = True
display = True
background_image = None
num_iterations = 1
num_processors = 1

#Parameters for sensitivity analysis
trajectory_file = "trajectories.txt"

r = 20 #Number of trajectories
levels = 10 #Number of breaks between min and max values
jumps = 2  #Number of jumps for oat steps
perc = 0.20 #Percentage to vary poarameters by

#Parameters to be varied
ptv = [("ips", 0, np.inf, int, [1]),
				("d50", 0, np.inf, float, range(3, 13)),
				("d99", 0, np.inf, float, range(3, 13)),
				("rep", 0, np.inf, float, range(3, 13)),
				("allee", 0, np.inf, int, [0]),
				("survivalT", 0, 1, float, range(0, 13)),
				("K", 1, np.inf, int, range(0, 4))]

#PBS resource requests
walltime = "12:00:00"
memory = "4GB"