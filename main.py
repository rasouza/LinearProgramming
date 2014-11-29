#coding: utf-8

import numpy as np
import pdb

def simplex(c, A, b, Basis):
	# Primeiro passo
	B = A[:, Basis]
	Binv = np.linalg.inv(B)
	x = np.linalg.solve(B,b)

	# Segundo passo
	p = np.dot(c[Basis], Binv)
	c_redu = np.zeros([len(A)]) # Inicializa o vetor de custos reduzidos
	for j in range(len(A)):
		c_redu[j] = c[j] - np.dot(p, A[:,j])

	# O algoritmo para se todos os custos são positivos
	parada = [True if x < 0 else False for x in c_redu ]
	if True not in parada:
		return [np.dot(c,x), x, p]

	# Pega o primeiro j não negativo
	for k,v in enumerate(c_redu):
		if v > 0:
			break
	# no fim do loop, k é o primeiro indice não negativo

	# Terceiro passo
	u = np.zeros([len(A)])
	for j in range(len(A)):
		u[j] = np.dot(Binv, A[:,j])

	# O algoritmo para se todos os u's são negativos
	parada = [True if x > 0 else False for x in u ]
	if True not in parada:
		# Neste caso, o custo otimo é -infinito
		return [-np.inf, None, None]

	# Quarto passo
	theta = x/u
	l, t_estrela = min(enumerate([x if x > 0 else np.inf for x in theta]))

	# Quinto passo
	x[l] = 0
	x[k] = t_estrela
	for i in Basis:
		x[i] = x[i] - t_estrela*u[i]

	# Sexto passo


	

A = np.matrix('1 2 2; 2 1 2; 0 1 2')

simplex(np.array([1,1,1]), A, [0,0,0], [0,1,2])