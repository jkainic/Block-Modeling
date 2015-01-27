from math import *
from copy import *
from numpy import random

def sum2(mat2d):
	total = 0
	for row in mat2d:
		for x in row:
			total += x
	return total

def sum1(vec):
	total = 0
	for el in vec:
		total += el
	return total

#make the stack
def concatenate(matrix_list):
	total_mat = []
	for mat in matrix_list:
		for row in mat:
			total_mat.append(row)
	return total_mat

class block_mat:
	def __init__(self, matrix_list, num_blocks = 1):
		#list of matrices in the stack
		self.matrix_list = matrix_list
		self.blocks = []
		self.num_blocks = num_blocks
		self.num_mats = len(matrix_list)
		self.matrix_names = [i for i in range(self.num_mats)]
		#for calculating T value
		self.smat = []
		self.dmat = [0]*self.num_mats
		self.T = 0
		if len(self.matrix_list) == 0:
			self.ntwk_size = 0
		else:
			self.ntwk_size = len(matrix_list[0])
		self.index = [0]*self.ntwk_size
		
		#if only looking for one block,
		#put entire network into that block
		if num_blocks == 1:
			self.blocks.append([])
			for i in range(self.ntwk_size):
				self.blocks[0].append(i)

		self.nu = [0]*self.num_mats
		for m in range(self.num_mats):
			sum_nu = sum2(self.matrix_list[m])
			x = sum_nu/float((self.ntwk_size*(self.ntwk_size-1)))
			self.nu[m] = x

	#update block_mat object when user enters matrices in stack
	def update_info(self):
		#ensure the matrices in the stack are in the proper form
		for i in range(self.num_mats):
			if not check_binary(self.matrix_list[i]):
				self.matrix_list[i] = find_binary(self.matrix_list[i], .15, .20)
            
		self.ntwk_size = len(self.matrix_list[0])
		self.num_mats = len(self.matrix_list)
		self.dmat = [0]*self.num_mats
		self.index = [0]*self.ntwk_size

		self.nu = [0]*self.num_mats
		for m in range(self.num_mats):
			sum_nu = sum2(self.matrix_list[m])
			x = sum_nu/float((self.ntwk_size*(self.ntwk_size-1)))
			self.nu[m] = x

		if self.num_blocks == 1:
			self.blocks.append([])
			for i in range(self.ntwk_size):
				self.blocks[0].append(i)
				self.index[i] = 0
			s = pow(self.ntwk_size, 2) - self.ntwk_size
			for m in range(self.num_mats):
				d = float(sum2(self.matrix_list[m]))/s
				self.T += s*pow(d - self.nu[m], 2)

	#run CONCOR algorithm to make a block model
	def concor(self):
		blocking = []
		nonzeroblock = []
		zeroblock = []

		total_mat = concatenate(self.matrix_list)

		#separate zero columns
		for x in range(self.ntwk_size):
			colx = [row[x] for row in total_mat]
			if sum1(colx) == 0:
				zeroblock.append(x)
			else:
				nonzeroblock.append(x)
		blocking.append(nonzeroblock)

		#split (nonzero column) blocks in the block model until there are the desired number of blocks
		while len(blocking)*2 + 1 < self.num_blocks:
			size = len(blocking)
			for x in range(size):
				block = blocking[0]
				if len(block) > 1:	
					twoblocks = split_block(total_mat, block)
					blocking.remove(block)
					for b in twoblocks:
						blocking.append(b)
		if len(zeroblock) != 0:
			blocking.append(zeroblock)

		while len(blocking) < self.num_blocks:
			block = blocking[0]
			twoblocks = split_block(total_mat, block)
			blocking.remove(block)
			for b in twoblocks:
				blocking.append(b)

		self.blocking(blocking)

	#update components of T-evaluation based on a particular block model
	def blocking(self, block_list):
		n = len(block_list)

		if n != self.num_blocks:
			print "There aren't the right number of blocks!"
			return False

		block_list_els = [x for block in block_list for x in block]
		block_list_els.sort()

		if len(block_list_els) != self.ntwk_size:
			print "You don't have the right number of elements!"
			return False

		for x in range(self.ntwk_size):
			if x != block_list_els[x]:
				print "You don't have the right elements!"
				return False

		self.blocks = block_list
		
		#keep track of which block each subject has been assigned
		for i in range(n):
			for el in block_list[i]:
				self.index[el] = i

		self.smat = s_ij(block_list)
		
		for m in range(self.num_mats):
			m_mat = self.matrix_list[m]
			self.dmat[m] = d_sum(m_mat, self.blocks, self.smat)

		self.T = T(self.nu, self.smat, self.dmat)

		return True

	#generate 100 random block models and return the mode with highest T-value...
	#it may be beneficial to adjust this code so that users can pass an arbitrary 
	#function by which to compare block model effectiveness
	def random_blocks(self):
		cut = int(self.ntwk_size / self.num_blocks)
		edge = self.num_blocks * cut

		#at first, just make a block model by dividing subjects fairly equally among blocks
		i = 0
		blocks = []
		for i in range(self.num_blocks):
			new_block = [n for n in range(cut*i, cut*(i+1))]
			blocks.append(new_block)

		final_part = [n for n in range(edge, self.ntwk_size)]
		blocks[self.num_blocks - 1] += final_part

		smat = s_ij(blocks)

		dmat = []
		for m in range(self.num_mats):
			m_mat = self.matrix_list[m]
			dmat.append(d_sum(m_mat, blocks, smat))

		T_value = T(self.nu, smat, dmat)
		optimal_blocks = deepcopy(blocks)

		#then, generate the random blocks and keep track of the "best" so far
		for x in range(100):
			new_blocks = []

			nums = random.permutation(self.ntwk_size)
			allnums = [x for x in nums]

			for i in range(self.num_blocks - 1):
				another_block = allnums[i*cut: cut*(i+1)]
				new_blocks.append(another_block)

			last_block = allnums[(self.num_blocks - 1)*cut: self.ntwk_size]
			new_blocks.append(last_block)

			new_smat = s_ij(new_blocks)

			dmat = []
			for m_mat in self.matrix_list:
				dmat.append(d_sum(m_mat, new_blocks, new_smat))

			new_T = T(self.nu, new_smat, dmat)

			if new_T > T_value:
				T_value = new_T
				optimal_blocks = deepcopy(new_blocks)

		self.blocking(optimal_blocks)

	#generate a single random block, for which agent are fairly equally divided among blocks
	def single_random_block(self):
		cut = int(self.ntwk_size / self.num_blocks)
		edge = self.num_blocks * cut

		new_blocks = []

		nums = random.permutation(self.ntwk_size)
		allnums = [x for x in nums]

		for i in range(self.num_blocks - 1):
			another_block = allnums[i*cut: cut*(i+1)]
			new_blocks.append(another_block)

		last_block = allnums[(self.num_blocks - 1)*cut: self.ntwk_size]
		new_blocks.append(last_block)

		return new_blocks

	#considers moving each agent to a new block, making changes that increase the optimization function value 
	#(in this case T), going through all agents until a local maximum is found
	def optimize(self):
		moves = 0
		prev_val = 0

		while self.T > prev_val:
			prev_val = self.T
			i = 0
			while i < self.ntwk_size:
				j = 0
				iblock = self.index[i]

				while j < self.num_blocks:

					if len(self.blocks[iblock]) == 1:
						j = self.num_blocks

					elif iblock == j:
						j += 1

					else:
						new_blocks = deepcopy(self.blocks)
						new_blocks[iblock].remove(i)
						new_blocks[j].append(i)

						new_self = block_mat(self.matrix_list, self.num_blocks)
						new_self.blocking(new_blocks)

						if new_self.T > self.T:
							self.blocking(new_blocks)
							moves += 1
							j = self.num_blocks

						else:
							j += 1

				i += 1

		return moves

#calculate T value based on a particular block model
def T(nu, s_matrix, d_matrix):
	T_val = 0
	num_blocks = len(s_matrix)

	for m in range(len(d_matrix)):
		for i in range(num_blocks):
			for j in range(num_blocks):
				if i != j:
					T_val += s_matrix[i][j]*pow((d_matrix[m][i][j] - nu[m]),2)

	return T_val
		
#returns a matrix whose ij^th entry is the s value of blocks #i and #j
def s_ij(blocking):
	l = len(blocking)
	S = []

	for i in range(l):
		len1 = len(blocking[i])
		row = [0]*l
		for j in range(l):
			len2 = len(blocking[j])

			row[j] = len1*len2
			if i == j:
				row[j] -= len1
		S.append(row)

	return S

#returns a matrix whose ij^th entry is the d value of blocks #i and #j 
#(from a particular tie matrix in the stack)
def d_sum(mat, blocking, smat):
	l = len(blocking)
	D = []

	for i in range(l):
		row = [0]*l
		for j in range(l):
			sm = 0
			if i != j:
				for p in blocking[i]:
					for q in blocking[j]:
						sm += mat[p][q]
				row[j] = float(sm) / smat[i][j]
		D.append(row)

	return D

#splits one block according to the CONCOR algorithm
#by repeatedly creating matrices whose entries come from
#calculating the pearson r coefficient between columns of the prior matrix,
#beginning with the stack of tie matrices
def split_block(mat, block):
	twoblocks = []

	size = len(block)
	M = []

	for i in block:
		M_row = []
		for j in block: 
			coli = [row[i] for row in mat]
			colj = [row[j] for row in mat]
			M_row.append(pearsonr(coli, colj))
		M.append(M_row)

	rep = 1

	while rep > 0:
		M_new = []
		for i in range(size):
			new_row = [0]*size
			for j in range(size):
				coli = [row[i] for row in M]
				colj = [row[j] for row in M]
				new_row[j] = pearsonr(coli, colj)
			M_new.append(new_row)
		
		M = deepcopy(M_new)

		rep = 0

		for x in range(size):
			for y in range(size):
				if abs(M_new[x][y]) < 0.99 and x != y:
					rep += 1
					break


	block1 = []
	block2 = []

	block1.append(block[0])

	for x in range(1,size):
		if M[0][x] > 0:
			block1.append(block[x])
		else:
			block2.append(block[x])

	twoblocks.append(block1)
	twoblocks.append(block2)

	return twoblocks

#calculates the pearson r coefficient between two vectors
def pearsonr(x, y):
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
