import numpy as np

def simplex(c, A, b, Basis):
	# Primeiro passo
	B = A[:, Basis]
	Binv = np.linalg.inv(B)

	# Segundo passo
	p = np.dot(c[Basis], Binv)
	c_redu = np.zeros(len(c))
	for j in range(len(A)):
		c_redu[j] = c[j] - np.dot(p, A[:,j]) 

	return

A = np.matrix('1 2 2; 2 1 2; 0 1 2')

d,p = simplex(np.array([1,1,1]), A, None, [0,1,2])