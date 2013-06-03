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
import PIL as PIL
import PIL.Image
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.collections as clt
import random as random
import numpy as np
import sys as sys
import re as re
import os as os

#Set up the colors for drawing. Populations are drawn as
#a blue trinagle, with larger populations drawn using a
#darker shade. Pyplot has a red palette that would be suitable,
#but people who are color blind have trouble distinguishing
#red from green, so we convert the red palette to blue 
red_dict = plt.cm.gist_heat_r._segmentdata
blue_dict = {'red':red_dict['blue'], 'green':red_dict['green'], 'blue':red_dict['red']}

#The landscape is draw usinf the YlGn built in color map. We
#want NaNs to be drawn as cream though, so set bad to cream
plt.cm.YlGn.set_bad("#FEFEEE", alpha=1.0)
plt.cm.gist_heat.set_bad("#FEFEEE", alpha=0.0)

#The landscape is drawn using a green palette
green_dict = {
	'red':[(0, 0, 95 / 255.),
			(1., 95 / 255., 95 / 255.)],
	'green':[(0, 0, 110 / 255.),
			(1., 180 / 255., 180 / 255.)],
	'blue':[(0, 0, 0 / 255.),
			(1., 0, 0.)]}

landscape_dict = {
	'red':[(0, 200 / 255., 200 / 255.),
			(0.1, 200 / 255., 200 / 255.),
			(0.2, 150 / 255., 150 / 255.),
			(0.3, 100 / 255., 100 / 255.),
			(0.4, 0., 0.),
			(1., 0., 0.)],
	'green':[(0, 220 / 255., 220 / 255.),
			(0.1, 220 / 255., 220 / 255.),
			(0.2, 180 / 255., 180 / 255.),
			(0.3, 180 / 255., 180 / 255.),
			(0.4, 180 / 255., 180 / 255.),
			(0.5, 160 / 255., 160 / 255.),
			(0.6, 140 / 255., 140 / 255.),
			(0.7, 120 / 255., 120 / 255.),
			(0.8, 100 / 255., 100 / 255.),
			(0.9, 60 / 255., 60 / 255.),
			(1., 60 / 255., 60 / 255.)],
	'blue':[(0, 0., 0.),
			(1., 0, 0.)]}



#Create the color maps from the dictionaries
blues = colors.LinearSegmentedColormap("blues", blue_dict, 64)
#pop_norm = colors.normalize(vmin=0, vmax=1)
#pop_norm = colors.LogNorm(vmin=0, vmax=1)
#pop_cmap = plt.cm.ScalarMappable(norm=pop_norm, cmap=plt.cm.gist_heat)
#pop_cmap = plt.cm.ScalarMappable(norm=pop_norm, cmap=blues)

greens = colors.LinearSegmentedColormap("greens", landscape_dict, 64)
greens.set_bad("#FEFEEE", alpha=1.0)


#Wrapper for trajectories and their associated results.
#Each trajectory consists of a number of points, and
#for each point a set of results.
class Trajectory(object):
	def __init__(self, points=None, num_sentinels=0):
		if points == None: points = []
		
		self._points = points
		self._num_metrics = 5 + num_sentinels

	#Add a point to the trajectory		
	def add_point(self, point):
		self._points.append(point)
	
	#Return the ith point
	def point_at_index(self, i):
		return self._points[i]

	#Return the points in the trajectory		
	def points(self):
		return self._points

	#Return the dimensions of the trajectory,
	#(num_points, num_metrics, num_iterations)			
	def dims(self):
		p, n, i = (0, 0, 0)
		if len(self._points) > 0:
			p = len(self._points) - 1
			i = self._points[0].num_iterations()
			n = self._num_metrics
			
		return (p, n, i)

	#Convert the trajectory to a matrix, where each row
	#represents a point, with the first p columns giving
	#the independent variables, and the last n columns
	#give the dependent variables
	def to_matrix(self):
		p, n, i = self.dims()
		mtrx = np.zeros((p + 1, p + n))
		for j, point in enumerate(self._points):
			mtrx[j,:p] = point.independents()
			mtrx[j,p:] = point.mean()
	
		return mtrx

	#Merge the results from another trajectory object. This
	#assumes that both trajectories have the same points
	def merge(self, t):
		for p1, p2 in zip(self._points, t.points()):
			p1._dependents += p2.dependents()


class Point(object):
	def __init__(self, independents=[], dependents=[], num_metrics=5):
		if independents == None: independents = []
		if dependents == None: dependents = []
		
		self._independents = independents
		self._dependents = dependents
		self._num_metrics = num_metrics
		
	def add_iteration(self):
		self._dependents.append([])
	
	def add_time_step(self, step):
		self._dependents[-1].append(step)
		
	def dependents(self):
		return self._dependents
	
	def independents(self):
		return self._independents
		
	def num_iterations(self):
		return len(self._dependents)

	def equals(self, p):
		return all([v1 == v2 for v1, v2 in zip(self._independents, p.independents())])

	#Return the mean for each metric across the
	#final time step of all iterations
	def mean(self):
		#The number of metrics tracked is given by the length of the first time step
		means = np.zeros(len(self._dependents[0][0]))

		#Go through the iterations, and add the 
		#final set of metrics from each to the total
		for it in self._dependents:
				means += np.array(it[-1])

		return means / self.num_iterations()


def population_statistics(ps, landscape, N):    
	nx, ny = landscape.shape	

	#Proportion of suitable habitat that's infested
	infested = np.sum(N > 0)
	infested /= float(np.sum(landscape > 0))

	#Distances travelled
	distances = []
	for i in xrange(nx):
		for j in xrange(ny):
			if N[i,j] > 0:
				dists = [np.linalg.norm([i-x, j-y]) for x, y in ps.isl]
				distances.append(min(dists))
	
	if len(distances) > 0:
		distances = sorted(distances)
		max_dist = ps.cell_size * distances[-1]
		med_dist = ps.cell_size * distances[len(distances) / 2]
	else:
		max_dist = 0
		med_dist = 0
		
	#Population totals for sentinel cells
	if ps.sentinels != None:
		sen_pops = [N[i,j] for i, j in ps.sentinels]
	else:
		sen_pops = []

	return [infested, max_dist, med_dist, np.sum(N)] + sen_pops


def process_iteration(Ns, ps, landscape, config):
	output_dir = config.output_dir + config.ext
	
	if config.background_image != None:
		background_path = config.input_dir + "/" + config.background_image
	else:
		background_path = None
	
	#Create a point to hold the iteration
	p = Point()
	p.add_iteration()
	
	#draw_population(Ns[0], landscape, ps.totalK, 0, output_dir, 2.0, background_path)
	
	if config.display:
		pool = mp.Pool(config.num_processors)

	for t in xrange(min(ps.max_time_steps, len(Ns))):
		if config.display:
			pool.apply_async(draw_population, [Ns[t], landscape, ps.totalK, t, output_dir, 2.0, background_path])
		
		p.add_time_step([t] + population_statistics(ps, landscape, Ns[t]))
	
	pool.close()

	#Write the iteration results to file as a trajectory containing a single point
	write_trajectories([Trajectory(points=[p])], None, ps.sentinels, output_dir + "/results.txt")

	if config.save_time_steps:
		np.savez(output_dir + "/populations.npz", *Ns)

	pool.join()


#Write the given point to file
def write_point(p, fout):
	#Write the point header
	fout.write(">point:[" + ",".join([str(v) for v in p.independents()]) + "]\n")
	#If there are any results for the point, write them too
	for iteration in p.dependents():
		#Iteration header
		fout.write(">iteration\n")
		
		for step in iteration:
			fout.write(",".join([str(v) for v in step]) + "\n")


#Write the given trajectory to file
def write_trajectory(t, fout):
	#Write the trajectory header
	fout.write(">trajectory\n")
	for p in t.points():
		write_point(p, fout)
	fout.write("//\n")


#Write all trajectories in the list to file
def write_trajectories(ts, ptv, sentinels, filename, append=False):
	r = len(ts)
	p, n, iters = ts[0].dims()

	try:
		#Open the results file for writing
		if append and os.path.exists(filename):
			f = open(filename, "a")
		#If this is a new file, write the header lines
		else:
			f = open(filename, "w")
			
			#Dimensions, num. trajectories, num. points - 1, max. time steps, num. iterations per point
			f.write("#r=" + str(r) + ", p=" + str(p) + ", n=" + str(n) + ", i=" + str(iters) + "\n")
			
			#Parameters that vary between points
			if ptv != None:
				f.write("#ptv=" + _ptvstr(ptv) + "\n")
			else:
				f.write("#ptv=None\n")
			
			#Sentinel cell coordinates
			if sentinels != None:
				sentinel_string = ",".join([str(tuple(s)) for s in sentinels])
				f.write("#sentinels=[" + sentinel_string + "]\n")
			else:
				f.write("#sentinels=None\n")
		
		#Now write each trajectory
		for t in ts:
			write_trajectory(t, f)
			
		f.close()
	except IOError as e:
		print "Problem writing trajectories to file " + filename
		print e
		sys.exit()		


#Read trajectores from file, and place in a list
def read_all_trajectories(filename):
	ts = []
	ptv = None
	sentinels = []
	curr_point = None
	
	try:
		fin = open(filename)
		for line in fin:
			#Header line 1, dimensions
			if re.search("#.*r=", line):
				#Read the number of trajectories, the number of parameters that
				#were varied, the number of dependent variables, and the number
				#of iterations that were performed for each point
				r, p, n, i = [eval(s) for s in re.findall("=(\d+)", line)]				
			#Header line 2, parameters varied
			elif re.search("#ptv=", line):
				expr = re.sub("#ptv=", "", line).strip()
				ptv = eval(expr)
			#Header line 3, sentinel cells
			elif re.search("#sentinels=", line):
				expr = re.sub("#sentinels=", "", line).strip()
				sentinels = eval(expr)
			#Trajectory header
			elif re.search("^>trajectory", line):
				ts.append(Trajectory(points=[], num_sentinels=len(sentinels)))
			#Point header
			elif re.search("^>point:", line):
				#Create a new point to represent the results that follow and add it to the current trajectory
				curr_point = Point(np.array(eval(re.sub(">point:", "", line).strip())), dependents=[])
				ts[-1].add_point(curr_point)
			#Iteration header
			elif re.search(">iteration", line):
				curr_point.add_iteration()
			#Time step
			elif re.search("^\d", line):
				step = tuple([eval(s) for s in line.split(",")])
				curr_point.add_time_step(step)
		fin.close()
	except IOError as e:
		print "Problem reading trajectories from file " + filename
		print e
		sys.exit()

	return (ts, ptv, sentinels)


def read_trajectories_gen(filename):
	ptv = None
	sentinels = []
	curr_point = None
	trajectory = None

	try:
		fin = open(filename)
		for line in fin:
			#Header line 1, dimensions
			if re.search("#.*r=", line):
				#Read the number of trajectories, the number of parameters that
				#were varied, the number of dependent variables, and the number
				#of iterations that were performed for each point
				r, p, n, i = [eval(s) for s in re.findall("=(\d+)", line)]				
			#Header line 2, parameters varied
			elif re.search("#ptv=", line):
				expr = re.sub("#ptv=", "", line).strip()
				ptv = eval(expr)
			#Header line 3, sentinel cells
			elif re.search("#sentinels=", line):
				expr = re.sub("#sentinels=", "", line).strip()
				sentinels = eval(expr)
			#Trajectory header
			elif re.search("^>trajectory", line):
				if trajectory == None:
					trajectory = Trajectory(points=[], num_sentinels=len(sentinels))
				else:
					yield (trajectory, ptv, sentinels)
					trajectory = Trajectory(points=[], num_sentinels=len(sentinels))
			#Point header
			elif re.search("^>point:", line):
				#Create a new point to represent the results that follow and add it to the current trajectory
				curr_point = Point(np.array(eval(re.sub(">point:", "", line).strip())), dependents=[])
				trajectory.add_point(curr_point)
			#Iteration header
			elif re.search(">iteration", line):
				curr_point.add_iteration()
			#Time step
			elif re.search("^\d", line):
				step = tuple([eval(s) for s in line.split(",")])
				curr_point.add_time_step(step)
		fin.close()
		yield (trajectory, ptv, sentinels)
	except IOError as e:
		print "Problem reading trajectories from file " + filename
		print e
		sys.exit()

	

def cat_trajectories(n, ptv, sentinels, output_dir):
	try:
		for i in xrange(n):
			rfile = output_dir + "/Trajectory" + str(i + 1) + "/results.txt"
			if os.path.isfile(rfile):
				res, _, _ = read_all_trajectories(rfile)
				write_trajectories(res, ptv, sentinels, output_dir + "/results.txt", append=(i != 0))		
	except IOError as e:
		print "Problem writing final results"
		print e
		

def draw_population(N, landscape, K, t, output_dir, scale=1, background_path=None):
	#Make sure the image directory exists
	img_dir = output_dir + "/Images"
	if not os.path.isdir(img_dir):
		os.makedirs(img_dir)
	
	#Load the background image if specified
	if background_path != None:
		background = PIL.Image.open(background_path)
		width, height = background.size
	else:
		background = None
		width, height = landscape.shape

	#Create the figure object
	figure = plt.figure(figsize=(width / 150 + 2, height / 150 + 1), dpi=150)
		
	#Draw the image and save to file
	_draw(N, landscape, K, scale, background, figure, width, height)
	plt.savefig(img_dir + "/step%.3d.png" % t)
	figure.clf()
	plt.close()
	del figure


#Display a population and landscape
def _draw(N, landscape, K, scale, background, figure, width, height): 
	pop_norm = colors.LogNorm(vmin=1, vmax=K)
	pop_cmap = plt.cm.ScalarMappable(norm=pop_norm, cmap=plt.cm.gist_heat)

	#Scale the triangles based on the landscape size
	xs = np.array([-scale / 2., 0, scale / 2.])
	ys = np.array([-scale / 2., scale / 2., -scale / 2.])
	
	plt.axis([0, width, 0, height])
	ax = figure.gca()
	
	#We want the axes to be undecorated
	ax.set_xticklabels([])
	ax.set_yticklabels([])
	ax.set_xticks([])
	ax.set_yticks([])

	#Draw the landscape, we need to rotate the landscape by 90 degrees
	#to have everything draw correctly. Taking the sqrt of the landscape
	#increases the contrast between the different habitat types, setting
	#the zero values to NaN means they get drawn in the error color (white)
	tl = np.sqrt(landscape)
	tl[tl == 0] = np.NaN;
	
	#Draw the background image is specified
	if background != None:
		ax.imshow(background, aspect='equal')
	else:
		ax.imshow(np.rot90(tl), cmap=greens, interpolation='nearest', vmin=0, vmax=1)
	
	#Draw the infested cells	using imshow
	im = ax.imshow(np.rot90(N), cmap=plt.cm.gist_heat, norm=pop_norm,
		extent=(0, width, 0, height), interpolation="nearest", aspect='auto')
	figure.colorbar(im)

	#Draw the infested cells as triangles. The way we iterate
	#through the cells is equivalent to rotating by 90 degrees
	#for i in xrange(nx):
	#	for j in xrange(ny):
	#		if N[i,j] > 0:
				#Calculate the color based on the population for the current cell
	#			col = colors.rgb2hex(pop_cmap.to_rgba(N[i,j]))
				
	#			plt.fill(xs + i, ys + (ny - j - 1), col, zorder=N[i,j], linewidth=0, antialiased=False)
	#			plt.fill(xs + i, ys + (ny - j - 1), col, zorder=N[i,j], fill=False)
		
	#Draw a color bar next to the plot
	#pop_cmap.set_clim((1, K))
	#pop_cmap.set_array(N)
	#figure.colorbar(pop_cmap)


#Plot a set of elementary effects from a sensitivity analysis
def plot_ees(ees, ptv, output_dir):
	r, m, p = ees.shape
	aprops = {'width': 1, 'frac': 0.3, 'headwidth': 3, 'shrink': 0.15}
	metrics = ["time", "infested", "max distance", "med distance", "count"]

	#Get the names of the paremters that were varied
	names = [name + " (" + str(i) + ")" for name, _, _, _, inds in ptv for i in inds]

	#Go through each metric
	for i, metric in enumerate(metrics):
		#Calculate the sensitivity statistics, considering only
		#those parameters that changed
		mu = np.mean(ees[:,i,:], axis=0)
		mu_star = np.mean(np.abs(ees[:,i,:]), axis=0)
		sigma = np.std(ees[:,i,:], axis=0)

		#Plot the figure
		figure = plt.figure(figsize=(8, 8), dpi=80)
		plt.plot(mu_star, sigma, "o", figure=figure)
		plt.xlim([-0.1, np.max(mu_star) + 0.5])
		plt.ylim([-0.1, np.max(sigma) + 0.5])
		
		for j, name in enumerate(names):
			#Get the point coordinates
			x = mu_star[j]
			y = sigma[j]

			#Annotate with an arrow which lies in a
			#random direction from the point
			plt.annotate(name, xy=(x, y), size="x-small",
				xytext=(random.choice([10, 20]), random.choice([10, 20])),
					textcoords="offset points", arrowprops=aprops)

		#Save the figure to file
		plt.savefig(output_dir + "/" + metric + ".png")


def array_to_list_string(a):
	return "[" + ",".join([str(v) for v in a]) + "]"


def _trajectories_to_matrices(ts):
	#Get the dimensions for the matrix, (r * (p + 1) x p), where r is the 
	#number of trajectories, and p + 1 is the number of points in each 
	#trajectory. The function dims also returns n, the number of metrics,
	#and i, the number of iterations performed for each point
	r = len(ts)
	p, n, i = ts[0].dims()
	
	#Construct the trajectory and results matrices
	x = np.zeros((r * (p + 1), p))
	y = np.zeros((r * (p + 1), n)) 
	
	#Copy the data from the trajectory objects to the matrices
	for i, t in enumerate(ts):
		tmtrx = t.to_matrix()
		index = i * (p + 1)
		x[index:(index+p+1),:] = tmtrx[:,:p]
		y[index:(index+p+1),:] = tmtrx[:,p:]
	
	return x, y
	

#Convert a parameters to vary dictionary to a string
def _ptvstr(ptv):
	string = "["
	for key, pmin, pmax, ptype, length in ptv:
		string += "('" + key + "', " + str(pmin) + ", "

		if np.inf == pmax:
			string += "np.inf, "
		else:
			string += str(pmax) + ", "

		if ptype == int:
			string += "int, "
		elif ptype == float: 
			string += "float, "
			
		string += str(length) + "),"
	
	return string[:-1] + "]"
