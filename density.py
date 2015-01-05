from numpy import *
from math import *
from copy import *

#used to find an appropriate cut off for data that is not already binary
def find_binary(mat, bottom, top):
	mx = max(x for row in mat for x in row)

	for n in range(1, mx + 1):
		poss_mat = to_binary(mat, n)
		dense = density(poss_mat)
		if bottom <= dense <= top:	
			print "Cutoff Density: " + str(dense)
			break
		elif n == mx:
			poss_mat = find_binary(mat, bottom - .05, top + .05)

	return poss_mat

#check that a tie matrix is in binary form
def check_binary(mat):
	flat_mat = [el for row in mat for el in row]
	for i in flat_mat:
		if i != 0 or i != 1:
			return False
	return True

#based on a particular cut-off, return the resulting binary data matrix
def to_binary(mat, cutoff):
	N = len(mat)
	ones_mat = []

	for i in range(N):
		row = [0]*N
		for j in range(N):
			if mat[i][j] >= cutoff:
				row[j] = 1
		ones_mat.append(row)

	return ones_mat

#find the density of ties in a matrix (assuming agents in and out are the same)
def density(mat):
	N = len(mat)
	ties = sum(mat)

	if N == 1:
		d = 0
	else:
		d = float(ties)/(N*(N-1))

	return d

#find the density of ties between two different sets of agents
def non_sq_density(mat):
	N = len(mat)
	M = len(mat[0])

	ties = sum(mat)

	d = float(ties)/(N*M)

	return d

#used for developing the matrices that display the density between blocks
def block_dense(blocks, blocked_matrix):
	N = len(blocks)

	density_mat = []

	place1 = 0
	place2 = 0

	for i in range(N):
		row = [0]*N
		n = len(blocks[i])

		for j in range(N):
			m = len(blocks[j])

			sub_mat = []

			for x in range(n):
				xrow = [0]*m
				for y in range(m):
					xrow[y] = blocked_matrix[x + place1][y + place2]
				sub_mat.append(xrow)

			if i == j:
				d = density(sub_mat)
			else:
				d = non_sq_density(sub_mat)

			row[j] = d
			place2 += m
		density_mat.append(row)

		place2 = 0
		place1 += n

	return density_mat

