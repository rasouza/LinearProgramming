#coding: utf-8

import numpy as np
import pdb

def escalonar(A, u, l, oldBasis, newBasis):

	A[l,:] = A[l,:]/u[l]
	u[l] = 1

	# Aplica a inversa [Binv | u]
	for i in range(len(u)):
		if i != l:
			A[i, :] = A[i, :] - A[l,:]*float(u[i])

	# Aplica as permutações da base
	for x in oldBasis:
		i = oldBasis.index(x)
		i_new = newBasis.index(x)
		if i_new != l:
			A[[i,i_new]] = A[[i_new,i]]

	return A

def simplex(c, A, b, Basis):
	# Primeiro passo
	B = A[:, Basis]
	Binv = np.linalg.inv(B)
	c = np.array(c)
	x = np.zeros(len(c))
	x[Basis] = np.linalg.solve(B,b)

	# print("----- PASSO 1 ------")
	# print "------ B -----------\n", B, "\n-------------------", "\n"
	# print "------ Binv -----------\n", Binv, "\n-------------------", "\n"
	# print "------ c -----------\n", c, "\n-------------------", "\n"
	# print "------ x -----------\n", x, "\n-------------------", "\n"


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


		# print("----- PASSO 2 ------")
		# print "------ Reduced cost -----------\n", c_redu, "\n-------------------", "\n"
		# print "------ J (getting in) -----------\n", k, "\n-------------------", "\n"

		# Terceiro passo
		u = np.dot(Binv, A[:,k])

		# O algoritmo para se todos os u's são negativos
		parada = [True if x_ > 0 else False for x_ in u ]
		if True not in parada:
			# Neste caso, o custo otimo é -infinito
			return [-np.inf, None, None]

		# print("----- PASSO 3 ------")
		# print "------ u -----------\n", u, "\n-------------------", "\n"

		# Quarto passo (Escolhe quem sai)
		u_b = list(zip(Basis, map(float,u)))
		theta = []
		for base,u_ in u_b:
			if u_ > 0:
				theta.append((base,x[base]/u_))
		l, t_estrela = min(theta, key = lambda x_: x_[1])

		# print("----- PASSO 4 ------")
		# print "------ U_b -----------\n", u_b, "\n-------------------", "\n"
		# print "------ L (getting out) -----------\n", l, "\n-------------------", "\n"
		# print "------ Theta -----------\n", theta, "\n-------------------", "\n"
		# print "------ Theta estrela -----------\n", t_estrela, "\n-------------------", "\n"

		# Quinto passo
		x[Basis] = np.matrix(x[Basis]).T - t_estrela*u
		x[l] = 0
		x[k] = t_estrela
		l = Basis.index(l)
		Basis[l] = k
		Basis_unsorted = Basis[:]
		Basis.sort()
		newB = A[:, Basis]

		# print("----- PASSO 5 ------")
		# print "------ xBasis -----------\n", x, "\n-------------------", "\n"
		# print "------ Basis -----------\n", Basis, "\n-------------------", "\n"
		# print "------ New B -----------\n", newB, "\n-------------------", "\n"

		# Sexto passo
		Binv = escalonar(Binv, u, l, Basis_unsorted, Basis)
		# Binv = np.linalg.inv(newB)

# c = [-10, -12, -12, 0, 0, 0]
# Basis = [3,4,5]
# A = np.matrix('1 2 2 1 0 0; 2 1 2 0 1 0; 2 2 1 0 0 1')
# b = [20, 20, 20]

c = [-100,-150,-50,0,0]
Basis = [3,4]
A = np.matrix('5 8 3 1 0; 2 2 1 0 1')
b = [50,20]

print(simplex(c,A,b,Basis))
