import sys as sys
import itertools as itertools

sys.path.append("./Modules")
import gmbiIO as io

if __name__ == "__main__":
	output_file = sys.argv[-1]

	for (t1, ptv, sentinels), (t2, _, _) in itertools.izip(io.read_trajectories_gen(sys.argv[1]), io.read_trajectories_gen(sys.argv[2])):
		t1.merge(t2)
		print "merged"
		io.write_trajectories([t1], ptv, sentinels, sys.argv[-1], append=True)