from numpy import *
import math
from copy import *

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
			sum_nu = sum(self.matrix_list[m])
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
			sum_nu = sum(self.matrix_list[m])
			x = sum_nu/float((self.ntwk_size*(self.ntwk_size-1)))
			self.nu[m] = x

		if self.num_blocks == 1:
			self.blocks.append([])
			for i in range(self.ntwk_size):
				self.blocks[0].append(i)
				self.index[i] = 0
			s = pow(self.ntwk_size, 2) - self.ntwk_size
			for m in range(self.num_mats):
				d = float(sum(self.matrix_list[m]))/s
				self.T += s*pow(d - self.nu[m], 2)

	#make the stack
	def concatenate(matrix_list):
		total_mat = []
		for mat in matrix_list:
			for row in mat:
				total_mat.append(row)
		return total_mat

	#run CONCOR algorithm to make a block model
	def concor(self):
		blocking = []
		nonzeroblock = []
		zeroblock = []

		total_mat = concatenate(self.matrix_list)

		#separate zero columns
		for x in range(self.ntwk_size):
			colx = [row[x] for row in total_mat]
			if sum(colx) == 0:
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

	#same as optimization, but attempts to do so in fewer steps by focusing at each stage
	#only on parts of the optimization function affected by the agent being considered...
	#it is unclear at this time if this could actually be a nontrivial speed-up, because it appears
	#that there would be a proportionate number of computations
	def optimize_by_parts(self, f_diff):
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
						n = f_diff(self.matrix_list, i, iblock, j, self.nu, self.smat, self.dmat, self.blocks)
						if n > 0:
							self.blocks[iblock].remove(i)
							self.blocks[j].append(i)
							self.index[i] = j
							self.T += n

							for k in range(self.num_blocks):
								if k != iblock and k != j:
									size = len(self.blocks[k])
									self.smat[iblock][k] -= size 
									self.smat[k][iblock] -= size
									self.smat[j][k] += size
									self.smat[k][j] += size

								elif k == iblock:
									self.smat[k][k] = pow(len(self.blocks[iblock]),2) - len(self.blocks[iblock])
									self.smat[k][j] = len(self.blocks[j])*len(self.blocks[k])
									self.smat[j][k] = len(self.blocks[j])*len(self.blocks[k])
								else:
									self.smat[k][k] = pow(len(self.blocks[j]),2) - len(self.blocks[j])

							for m in range(self.num_mats):
								self.dmat[m] = d_sum(self.matrix_list[m], self.blocks, self.smat)

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
	
#calculate part of the T value affected by moving an agent from one block to another:
#returns the difference... if positive, the move is beneficial
def T_diff(matrix_list, agent, block1, block2, nu, s_matrix, d_matrix, blocks):
	num_blocks = len(blocks)
	old_part = 0
	new_part = 0

	for m in range(len(d_matrix)):
		for i in range(num_blocks):
			if i != block1 and i != block2:
				size = len(blocks[i])
				s_i1 = s_matrix[i][block1] - size
				s_i2 = s_matrix[i][block2] + size

				sum11 = 0
				sum12 = 0
				sum21 = 0
				sum22 = 0
				for j in blocks[block1]:
					sum11 += matrix_list[m][agent][j]
					sum12 += matrix_list[m][j][agent]
				for j in blocks[block2]:
					sum21 += matrix_list[m][agent][j]
					sum22 += matrix_list[m][j][agent]

				d_i1 = float(d_matrix[m][i][block1]*s_matrix[i][block1] - sum11) / s_i1
				d_1i = float(d_matrix[m][i][block1]*s_matrix[i][block1] - sum12) / s_i1
				d_i2 = float(d_matrix[m][i][block2]*s_matrix[i][block2] + sum21) / s_i2
				d_2i = float(d_matrix[m][i][block2]*s_matrix[i][block2] + sum22) / s_i2

				old_part += s_matrix[i][block1]*(pow((d_matrix[m][i][block1] - nu[m]), 2) + pow((d_matrix[m][block1][i] - nu[m]), 2))
				old_part += s_matrix[i][block2]*(pow((d_matrix[m][i][block2] - nu[m]), 2) + pow((d_matrix[m][block2][i] - nu[m]), 2))

				new_part += s_i1*(pow((d_i1 - nu[m]), 2) + pow((d_1i - nu[m]), 2))
				new_part += s_i2*(pow((d_i2 - nu[m]), 2) + pow((d_2i - nu[m]), 2))

			elif i == block1:
				if len(blocks[block1]) > 1:
					old_part += s_matrix[block1][block1]*pow((d_matrix[m][block1][block1] - nu[m]), 2)

				if len(blocks[block1]) > 2:
					s_ii = pow((len(blocks[block1]) - 1), 2) - (len(blocks[block1]) - 1)
					sum_ii = 0
					for j in blocks[i]:
						sum_ii += matrix_list[m][agent][j] + matrix_list[m][j][agent]
					d_ii = float(d_matrix[m][block1][block1]*s_matrix[block1][block1] - sum_ii) / s_ii

					new_part += s_ii*pow((d_ii - nu[m]), 2)
			else:
				if len(blocks[block2]) > 1:
					s_i1 = (len(blocks[block2]) + 1)*(len(blocks[block1]) - 1)

				sub_i1 = 0
				sub_1i = 0
				add_i1 = 0
				add_1i = 0
				for j in blocks[block2]:
					sub_i1 += matrix_list[m][agent][j]
					sub_1i += matrix_list[m][j][agent]
				for j in blocks[block1]:
					add_i1 += matrix_list[m][agent][j]
					add_1i += matrix_list[m][j][agent]
				
				d_i1 = float(d_matrix[m][block2][block1]*s_matrix[block2][block1] + add_i1 - sub_i1)/s_i1
				d_1i =  float(d_matrix[m][block1][block2]*s_matrix[block1][block2] + add_1i - sub_1i)/s_i1

				old_part += s_matrix[block2][block2]*pow((d_matrix[m][block2][block2] - nu[m]), 2)
				old_part += s_matrix[block1][block2]*(pow((d_matrix[m][block1][block2] - nu[m]), 2) + pow((d_matrix[m][block2][block1] - nu[m]), 2))

				if len(blocks[block2]) > 2:
					s_ii = pow((len(blocks[block2]) + 1), 2) - (len(blocks[block2]) + 1)
					sum_ii = 0
					for j in blocks[i]:
						sum_ii += matrix_list[m][agent][j] + matrix_list[m][j][agent]
					d_ii = float(d_matrix[m][block2][block2]*s_matrix[block2][block2] + sum_ii) / s_ii
					new_part += s_ii*pow((d_ii - nu[m]), 2)

				new_part += s_i1*(pow((d_i1 - nu[m]), 2) + pow((d_1i - nu[m]), 2))

	return new_part - old_part
		
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
				if absolute(M_new[x][y]) < 0.99 and x != y:
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

#given a particular block model, rearranges the columns and rows of matrices in the stack
#so that agents in each block are adjacent to one another
def blocked_mats(blocks, matrix_list, size):
	blocked_mats = []
	N = len(blocks)

	for mat in matrix_list:
		blocked_tie_mat = deepcopy(mat)

		for x in range(N):
			block1 = blocks[x]
			n = len(block1)
			for i in range(n):
				s = block1[i]
				for y in range(N):
					block2 = blocks[y]
					m = len(block2)
					for j in range(m):
						t = block2[j]
						blocked_tie_mat[i][j] = mat[s][t]

		blocked_mats.append(blocked_tie_mat)

	return blocked_mats
