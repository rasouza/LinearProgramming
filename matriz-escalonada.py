import numpy

# Algoritmo para resolver sistemas lineares com matrizes escalonadas

# A.x = b
A = numpy.loadtxt("matrix.txt")
b = numpy.loadtxt("vector.txt")

print("A.x = b")
print(A)
print(b)

for i in range(b.size):
	for j in range(i-1):
		b[i] = b[i] - A[i][j]*b[j]

	if A[i][i] == 0:
		break
	b[i] = b[i]/A[i][i]

print(b)

