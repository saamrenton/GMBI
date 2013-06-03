import re as re
import sys as sys
import numpy as np
import argparse as ap
import random as random
from collections import Counter

neighbours4 = [(-1, 0), (1, 0), (0, -1), (0, 1)]
neighbours8 = neighbours4 + [(-1, 1), (-1, -1), (1, 1), (1, -1)]


def _river(size):
	nx, ny = size
	mny = int(ny / 2)
	landscape = 0.5 * np.ones(size)
	half_width = 5
	(b,d) = np.random.uniform(0.5,1,2)
	(a,c) = np.random.uniform(2,4,2)
	
	twopi = 2 * np.pi
	fourpi = 2 * twopi
	
	try:
		for i in range(nx):
			x = -twopi + fourpi * (i / float(nx))
			y = a * np.sin(b * x) + c * np.cos(d * x)
				
			j = int(ny / 3) + 5 * y
			#landscape[i,j-2*half_width:j+2*half_width] = 1			
			landscape[i,j-half_width:j+half_width] = 0
	except IndexError:
		print "River went too far"
	
	return landscape


def _fill_with_mrc(landscape, a, ptype, p, fill_type=0):
	mx, my = landscape.shape

	#Go through each cell in the landscape
	for i, j in np.ndindex(landscape.shape):
		if fill_type == landscape[i,j]:
			landscape[i,j] = -1 if np.random.binomial(1, 0.4) else -2
	
	if fill_type == landscape[mx / 2, my / 2]:
		landscape[mx / 2, my / 2] = -1
	_assign_cluster_types(landscape, a, ptype) #B and C	
	_fill_landscape(landscape, a, ptype) #D



def _fill_with_uniform(landscape, p, suitability, fill_type=0):
	#Go through each cell in the landscape
	for i, j in np.ndindex(landscape.shape):
		if fill_type == landscape[i,j]:
			patch_type = np.random.multinomial(1, p)
			landscape[i,j] = suitability[patch_type.tolist().index(1)]


#Generate a landscape where each pixel has a probability
#p of having one of the specifed suitabilities
def _uniform_landscape(size, p, suitability):
	nx, ny = size
	
	patch_types = np.random.multinomial(1, p, size)
	landscape = np.zeros(size)

	for i, j in np.ndindex((nx, ny)):
		landscape[i,j] = suitability[patch_types[i,j].tolist().index(1)]

	landscape[nx / 2, ny / 2] = suitability[-1]
	
	return landscape


#Add a randomly placed block to a landscape
def _add_block(x, y, area, max_y, max_x, suitability, landscape):
	if area > 9:
		width = np.round(np.random.uniform(0.8 * np.sqrt(area), 1.2 * np.sqrt(area)))
		height = np.round(area / width)
	else:
		width = np.round(np.sqrt(area))
		height = np.round(area / width)

	mx = int(width / 2)
	my = int(height / 2)

	for sx in (i % max_x for i in range(x - mx, x + mx)):
		for sy in (j % max_y for j in range(y - my, y + my)):
			landscape[sx, sy] = suitability


#Create a blocked landscape
def _blocked_landscape(landscape_size, p, area, suitability=1):
	mx, my = landscape_size
	landscape = np.zeros(landscape_size, float)

	#Generate a block around the centre
	_add_block(mx / 2, my / 2, area, mx, my, suitability, landscape)

	nb_suitable = p * mx * my
	while landscape.sum() < nb_suitable:
		i = int(round(np.random.uniform(0, mx)))
		j = int(round(np.random.uniform(0, my)))
		_add_block(i, j, area, mx, my, suitability, landscape)
		
	return landscape


#Implementation of MRC from Saura and Martinez-Milan (2000)
def _mrc_landscape(size, a, ptype, p):
	mx, my = size
	ptype = [s + 2. for s in ptype]

	landscape = np.array(np.random.binomial(1, 0.4, size), dtype=float) #A
	landscape[mx / 2, my / 2] = 1
	_assign_cluster_types(landscape, a, ptype) #B and C	
	_fill_landscape(landscape, a, ptype) #D

	return landscape - 2

	
def _assign_cluster_types(landscape, a, ptype):
	mx, my = landscape.shape
	
	#Process the middle cell, setting its type to be the
	#most suitable habitat available
	_map_cluster(landscape, mx / 2, my / 2, ptype[-1])

	#Go through the remaining cells
	for (i, j), x in np.ndenumerate(landscape):
		#If the cell is marked as suitable, but hasn't been assigned a type
		if -1 == x:
			#Assign a suitability to the current cell and then map out the cluster
			#ptype = np.random.choice(types, p=a) #only in 1.7
			t = ptype[np.random.multinomial(1, a).tolist().index(1)]
			_map_cluster(landscape, i, j, t)


#Map out a connected habitat cluster
def _map_cluster(landscape, i, j, t):
	landscape[i, j] = t
	mx, my = landscape.shape

	#Go through each of the neighbours	
	for di, dj in neighbours4:
		ni, nj = (i + di, j + dj)
		#If the neighbour is in the bounds and is suitable,
		#but hasn't hasn't yet been traversed
		if(0 <= ni < mx and 0 <= nj < my and landscape[ni, nj] == -1):
			#Map out the neighbours cluster
			_map_cluster(landscape, ni, nj, t)


#Fill in the background (non-suitable) of landscape
def _fill_landscape(landscape, a, ptype, fill_type=-2):
	mx, my = landscape.shape

	#Go through the cells in the landscape that remain unmarked
	unmarked = [ind for ind, x in np.ndenumerate(landscape) if x == fill_type]
	
	for (i, j) in np.random.permutation(unmarked):
		#Get the neighbouring cells
		ns = [landscape[i + di, j + dj] for 
			di, dj in neighbours8 if 0 <= i + di < mx and 0 <= j + dj < my]

		#If the neighbours are all unmarked, assign a random type
		if all([fill_type == n for n in ns]):
			t = ptype[np.random.multinomial(1, a).tolist().index(1)]
			landscape[i,j] = t
		#Otherwise assign the most common type
		else:
			[t] = Counter([nt for nt in ns if nt != fill_type]).most_common(1)
			landscape[i,j] = t[0]


#Generate a landscape of the specified type
def landscape(landscape_type, nx, ny, psuit, suitability, avg_cluster):
	if landscape_type.lower() == "paddocks_and_background":
		landscape = _blocked_landscape((nx, ny), psuit[-1], avg_cluster[-1], suitability[-1])
		p = [p / sum(psuit[0:-1]) for p in psuit[0:-1]]
		_fill_with_uniform(landscape, p, suitability[0:-1])
		return landscape
	elif landscape_type == "mrc":
		return _mrc_landscape((nx, ny), psuit, suitability, avg_cluster)
	elif landscape_type == "river_mrc":
		landscape = _river((nx, ny))
		_fill_with_mrc(landscape, psuit, suitability, avg_cluster, 0.5)
		landscape[nx / 2, ny / 2] = 1
		return landscape
	else:
		return _uniform_landscape((nx, ny), psuit, suitability)
		

#Load landscape from file
def landscape_from_file(filename):
	landscape = []
	
	try:
		for line in open(filename):
			landscape.append([float(s) for s in re.split("\s*,\s*", line)])
	except IOError as e:
		print "Problem reading landscape from file " + filename
		print e
		sys.exit()
	
	return np.array(landscape)


#Write a landscape to file
def write_landscape(landscape, filename):
	try:
		f = open(filename, "w")
		for row in landscape:
			f.write(",".join([str(v) for v in row]) + "\n")
	except IOError as e:
		print "Problem writing landscape to file " + filename
		print e
		sys.exit()
		
	f.close()
	

if __name__ == "__main__":
	parser = ap.ArgumentParser()
	parser.add_argument("--config", default="config.py")
	args = parser.parse_args()
	
	#Load the config file
	execfile(args.config)
	nx, ny = landscape_size
	
	ls = landscape(landscape_type, nx, ny, psuit, suitability, avg_cluster)
	write_landscape(ls, output_file)
		
	
	