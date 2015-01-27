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
	n = len(mat)
	m = len(mat[0])

	for i in range(n):
		for j in range(m):
			if mat[i][j] != 0 or mat[i][j] != 1:
				return False
	return True

#based on a particular cut-off, return the resulting binary data matrix
#assumes the matrix is square
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
	ties = 0
	for row in mat:
		for x in row:
			ties += x

	if N == 1:
		d = 0
	else:
		d = float(ties)/(N*(N-1))

	return d

#used for developing the matrices that display the density between blocks
def block_dense(blocks, matrix):
	density_mat = []
	N = len(blocks)

	for b1 in range(N):
		row = [0]*N
		for b2 in range(N):
			d = 0

			if b1 == b2 and len(blocks[b1]) != 1:
				block = blocks[b1]
				n = len(block)
				ties = 0
				for i in block:
					for j in block:
						ties += matrix[i][j]
				d = float(ties)/ (n*(n-1))

			elif b1 != b2:
				block1 = blocks[b1]
				block2 = blocks[b2]

				n = len(block1)
				m = len(block2)

				ties = 0

				for i in block1:
					for j in block2:
						ties += matrix[i][j] + matrix[j][i]

				d = float(ties)/(n*m)

			row[b2] = d

		density_mat.append(row)

	return density_mat





