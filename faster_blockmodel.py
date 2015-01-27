from math import *
from copy import *
from density2 import *
from numpy import *

#assume input matrix is in binary, and that matrix is square
def mat_to_ties(mat):
	N = len(mat)
	ties = {}
	for x in xrange(N):
		for y in xrange(N):
			if mat[y][x]:
				if x in ties.keys():
					ties[x].append(y)
				else:
					ties[x] = [y]
	return ties

class block_mat:
	def __init__ (self, pop_size, matrix_list, num_blocks = 1):
		self.stack = matrix_list
		self.num_blocks = num_blocks
		self.ntwk_size = pop_size
		self.T = 0
		
		self.blocks = {}
		for x in xrange(num_blocks):
			self.blocks.get(x,[])

		#keep track of what block each element is in
		self.index = {}
		for x in range(self.ntwk_size):
			self.index[x] = 0

		self.nu = {}
		for name in self.stack.keys():
			mat = self.stack[name]
			nu_m = 0
			for col in mat.values():
				nu_m += len(col)
			self.nu[name] = float(nu_m) / (self.ntwk_size*(self.ntwk_size - 1))

	def set_blocks(self, block_list):
		n = len(block_list)
		if n != self.num_blocks:
			raise ValueError("There aren't the right number of blocks!")

		num_els = 0
		for block in block_list.values():
			num_els += len(block)
		if num_els != self.ntwk_size:
			raise ValueError("Not all of the individuals in your network are represented in these blocks!")

		check = False
		for x in range(self.ntwk_size):
			for block in block_list.values():
				if x in block:
					check = True
					break
			if check == False:
				raise ValueError("There's an individual missing in your block list!")

		self.blocks = block_list

		#keep track of which block each subject has been assigned to
		for name in block_list.keys():
			for el in block_list[name]:
				self.index[el] = name

		self.T = T(self.stack, self.blocks, self.nu)

	def concor(self):
		blocks = {}
		blocks[0] = [x for x in range(self.ntwk_size)]

		#keep track of the zero columns
		zeros = []
		for x in xrange(self.ntwk_size):
			xcol = []
			i = 0
			for n in self.stack.keys():
				n_mat = self.stack[n]
				if x in n_mat:
					xcol += [(el + i*self.ntwk_size) for el in n_mat[x]]
				i += 1

			if len(xcol) == 0:
				zeros.append(x)
				blocks[0].remove(x)
		
		#split blocks in the block model until there are close the desired number of blocks
		while len(blocks)*2 + 1 < self.num_blocks:
			size = len(blocks)
			block_halves = {}
			for block_index in blocks.keys():
				if len(blocks[block_index]) > 1:
					block = blocks[block_index]
					two_blocks = split_block(self.stack, block, self.ntwk_size)
					blocks[block_index] = two_blocks[0]
					block_halves[block_index + size] = two_blocks[1]
			for key in block_halves.keys():
				blocks[key] = block_halves[key]

		i = 0
		while len(blocks) + 1 < self.num_blocks:
			size = len(blocks)
			if len(blocks[i]) > 1:
				block = blocks[i]
				two_blocks = split_block(self.stack, block, self.ntwk_size)
				blocks[i] = two_blocks[0]
				blocks[i + size] = two_blocks[1]
			i += 1

		#add the block of zero columns
		if len(zeros) > 0:
			blocks[-1] = zeros

		self.set_blocks(blocks)

	def single_random_block(self):
		cut = int(self.ntwk_size/self.num_blocks)
		edge = self.num_blocks*cut

		nums = random.permutation(self.ntwk_size)
		nums.tolist()
		new_blocks = {}

		for i in range(self.num_blocks - 1):
			another_block = nums[i*cut: cut*(i+1)]
			new_blocks[i] = another_block

		last_block = nums[(self.num_blocks - 1)*cut: self.ntwk_size]
		new_blocks[self.num_blocks-1] = last_block

		return new_blocks

	def optimize(self):
		new_blocks = optimize(self.ntwk_size, self.index, self.blocks, self.stack, self.nu)
		self.set_blocks(new_blocks)

def s_ij(block_list):
	s = []
	n = len(block_list)
	for i in range(n):
		bi = len(block_list[i])
		rowi = [0]*n
		for j in range(n):
			bj = len(block_list[j])
			rowi[j] = bi*bj
			if i == j:
				rowi[j] -= bi
		s.append(rowi)

	return s

def d_ij(matrix, s_mat, block_list):
	n = len(block_list)
	d_mat = []

	for i in range(n):
		rowi = [0]*n
		for j in range(n):
			if i == j:
				rowi[j] = 0
			else:
				total = 0
				for p in block_list[i]:
					for q in block_list[j]:
						if p in matrix[q]:
							total += 1
				rowi[j] = float(total) / s_mat[i][j]
		d_mat.append(rowi)

	return d_mat

def T(matrix_list, block_list, nu):
	T_val = 0
	num_blocks = len(block_list)
	s_mat = s_ij(block_list)

	for name in matrix_list.keys():
		d_mat = d_ij(matrix_list[name], s_mat, block_list)

		for i in block_list.keys():
			for j in block_list.keys():
				T_val += s_mat[i][j]*pow((d_mat[i][j] - nu[name]), 2)

	return T_val

def optimize(ntwk_size, index_list, blocks, matrix_list, nu):
	moves = 0
	T_val = T(matrix_list, blocks, nu)
	prev_val = 0
	num_blocks = len(blocks)

	while T_val > prev_val:
		prev_val = T_val
		for i in xrange(ntwk_size):
			i_index = index_list[i]
			for j in range(num_blocks):
				if len(blocks[i_index]) > 1:
					break
				elif i_index != j:
					new_blocks = {}
					for key in blocks.keys():
						new_blocks[key] = blocks[key]
						if key == i_index:
							new_blocks[key].remove(i)
						elif key == j:
							new_blocks[key].append(i)

					new_T = T(matrix_list, new_blocks, nu)
					if new_T > T_val:
						blocks = new_blocks
						T_val = new_T
						moves += 1
						break

		print "# Moves to optimize: " + str(moves)
		return new_blocks


def split_block(mat_ties, block, ntwk_size):
	twoblocks = []
	size = len(block)

	#the M matrices here aren't going to be sparse,
	#which is why I use a list of lists to represent the matrix
	#rather than a dictionary of columns as before
	M = []

	for i in block:
		i_row = [0]*size

		#here, a column will have as elements the indices of its non-zero values
		coli = []
		x = 0
		for n in mat_ties.keys():
			n_mat = mat_ties[n]
			if i in n_mat:
				coli += [(el + x*ntwk_size) for el in n_mat[i]]
			x += 1

		for j_index in range(size):
			j = block[j_index]
			colj = []
			x = 0
			for n in mat_ties.keys():
				n_mat = mat_ties[n]
				if j in n_mat:
					colj += [(el + x*ntwk_size) for el in n_mat[j]]
				x += 1

			#new form of a column requires a different form of pearsonr
			i_row[j_index] = pearson_r2(coli, colj, ntwk_size)

		M.append(i_row)

	rep = 1

	while rep > 0:
		M_new = []

		for i in range(size):
			new_i_row = [0]*size

			coli = [row[i] for row in M]
			for j in range(size):
				colj = [row[j] for row in M]
				new_i_row[j] = pearson_r(coli, colj)

			M_new.append(new_i_row)

		M = deepcopy(M_new)

		rep = 0

		for x in range(size):
			for y in range(size):
				if abs(M_new[x][y]) < 0.99 and x != y:
					rep = 1
					break

	block1 = []
	block2 = []

	block1.append(block[0])

	for x in range(1, size):
		if M[0][x] > 0:
			block1.append(block[x])
		else:
			block2.append(block[x])

	twoblocks.append(block1)
	twoblocks.append(block2)

	return twoblocks


def pearson_r(x, y):
	if len(x) != len(y):
		print "These aren't the same length!"

	else:
		n = len(x)
		xbar = float(sum(x)) / n
    	ybar = float(sum(y)) / n
    	numerator = 0
    	xdenom = 0
    	ydenom = 0
    	for i in range(n):
    		xdiff = x[i] - xbar
    		ydiff = y[i] - ybar
    		numerator += xdiff * ydiff
    		xdenom += pow(xdiff, 2)
    		ydenom += pow(ydiff, 2)
		val = numerator / sqrt(xdenom * ydenom)

	return val

def pearson_r2(col1, col2, n):
	xbar = float(len(col1))/n
	ybar = float(len(col2))/n

	numerator = 0
	xdenom = 0
	ydenom = 0

	for i in range(n):
		xdiff = (i in col1) - xbar
		ydiff = (i in col2) - xbar
		numerator += xdiff*ydiff
		xdenom += pow(xdiff, 2)
		ydenom += pow(ydiff, 2)
	val = numerator/sqrt(xdenom*ydenom)

	return val

