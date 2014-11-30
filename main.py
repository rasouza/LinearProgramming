#coding: utf-8

import numpy as np
import pdb

def escalonar(A, u, l):
	

	A[l,:] = A[l,:]/u[l]
	u[l] = 1

	for i in range(len(u)):
		if i != l:
			A[i, :] = A[i, :] - A[l,:]*int(u[i])
			u[i] = u[i] - u[l]*u[i]

	return A

def simplex(c, A, b, Basis):
	# Primeiro passo
	B = A[:, Basis]
	Binv = np.linalg.inv(B)
	c = np.array(c)
	x = np.zeros(len(c))
	x[Basis] = np.linalg.solve(B,b)
	
	while True:
		# Segundo passo (Escolhe quem entra)
		p = np.dot(c[Basis], Binv)
		c_redu = np.zeros([len(A.T)]) # Inicializa o vetor de custos reduzidos
		
		for j in range(len(A.T)):
			if j not in Basis:
				c_redu[j] = c[j] - np.dot(p, A[:,j])
		
		# O algoritmo para se todos os custos são positivos
		parada = [True if x_ < 0 else False for x_ in c_redu ]

		if True not in parada:
			pdb.set_trace()
			return [np.dot(c,x), x, p]

		# Pega o primeiro j não negativo
		for k,v in enumerate(c_redu):
			if v < 0:
				break
		# no fim do loop, k é o primeiro indice não negativo
		

		# Terceiro passo
		u = np.dot(Binv, A[:,k])

		# O algoritmo para se todos os u's são negativos
		parada = [True if x_ > 0 else False for x_ in u ]		
		if True not in parada:
			# Neste caso, o custo otimo é -infinito
			return [-np.inf, None, None]

		# Quarto passo (Escolhe quem sai)
		u_b = list(zip(Basis, map(float,u)))
		theta = []
		for base,u_ in u_b:
			if u_ > 0:
				theta.append((base,x[base]/u_))
		l, t_estrela = min(theta, key = lambda x_: x_[1])
		
		# Quinto passo
		x[Basis] = np.matrix(x[Basis]).T - t_estrela*u
		x[l] = 0
		x[k] = t_estrela
		Basis[Basis.index(l)] = k
		Basis.sort()
		
		# Sexto passo
		Binv = escalonar(Binv, u, Basis.index(k))
		pdb.set_trace()

	

A = np.matrix('-1 -2 2 1; -2 1 2 0; 0 1 2 0', dtype=float)
u = np.array([2,1,3], dtype=float)

#print(simplex(np.array([1,1,1,1]), A, [0,0,1], [0,1,2]))


c = [-10, -12, -12, 0, 0, 0]
Basis = [3,4,5]
A = np.matrix('1 2 2 1 0 0; 2 1 2 0 1 0; 2 2 1 0 0 1')
b = [20, 20, 20]
print(simplex(c,A,b,Basis))