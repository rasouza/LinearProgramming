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
			return [np.dot(c,x), x, p]

		# Pega o primeiro j não negativo
		for k,v in enumerate(c_redu):
			if v < 0:
				break
		# no fim do loop, k é o primeiro indice não negativo
		

		# Terceiro passo
		u = np.dot(Binv, A[:,k])
		u[u == 0] = -np.inf # CUIDADO: tratando divisões por zero

		# O algoritmo para se todos os u's são negativos
		parada = [True if x_ > 0 else False for x_ in u ]
		
		if True not in parada:
			# Neste caso, o custo otimo é -infinito
			return [-np.inf, None, None]

		# Quarto passo (Escolhe quem sai)
		theta = np.array(x[Basis]/u.T)[0]
		l, t_estrela = min(enumerate([x_ if x_ > 0 else np.inf for x_ in theta]), key = lambda x: x[1])
		pdb.set_trace()
		l = Basis[l]
		

		u_ = iter(u)
		# Quinto passo
		u[u == -np.inf] = 0
		x[Basis] = np.matrix(x[Basis]) - t_estrela*u
		# for i in Basis:
		# 	if i != l:
		# 		x[i] = x[i] - t_estrela*u_.next()
		x[l] = 0
		x[k] = t_estrela
		Basis[Basis.index(l)] = k
		Basis.sort()
		
		
		# Sexto passo
		Binv = escalonar(Binv, u, Basis.index(k))
	

	

A = np.matrix('-1 -2 2 1; -2 1 2 0; 0 1 2 0', dtype=float)
u = np.array([2,1,3], dtype=float)

#print(simplex(np.array([1,1,1,1]), A, [0,0,1], [0,1,2]))


c = [-100, -150, -50, 0, 0]
Basis = [3,4]
A = np.matrix('5 8 3 1 0; 2 2 1 0 1')
b = [50, 20]
print(simplex(c,A,b,Basis))