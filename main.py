#coding: utf-8

import numpy as np
import pdb

def escalonar(A, u, l):
	

	A[l,:] = A[l,:]/u[l]
	u[l] = 1

	for i in range(len(u)):
		if i != l:
			A[i, :] = A[i, :] - A[l,:]*u[i]
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
		u = np.dot(Binv, A[:,j])
		u[u == 0] = -np.inf # CUIDADO: tratando divisões por zero

		# O algoritmo para se todos os u's são negativos
		parada = [True if x_ > 0 else False for x_ in u ]
		
		if True not in parada:
			# Neste caso, o custo otimo é -infinito
			return [-np.inf, None, None]

		# Quarto passo (Escolhe quem sai)
		theta = np.array(x[Basis]/u.T)[0]
		l, t_estrela = min(enumerate([x_ if x_ > 0 else np.inf for x_ in theta]))
		

		# Quinto passo
		for i in Basis:
			if i != l:
				x[i] = x[i] - t_estrela*u[i]
		x[l] = 0
		x[k] = t_estrela
		Basis[l] = k
		Basis.sort()
		
		# Sexto passo
		Binv = escalonar(Binv, u, l)
	

	

A = np.matrix('-1 -2 2 1; -2 1 2 0; 0 1 2 0', dtype=float)
u = np.array([2,1,3], dtype=float)

print(simplex(np.array([1,1,1,1]), A, [0,0,1], [0,1,2]))
