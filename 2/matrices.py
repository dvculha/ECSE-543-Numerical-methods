import math
import numpy as np
import time
import random

# returns True if matrix A is symmetric
def is_symmetric(A):
    n = len(A)
    if (n==1): 
        return True
    for i in range(n):
        for j in range(i + 1, n):
            return A[i][j] == A[j][i]
 
#A = [[0.2]]
#print(is_symmetric(A))

#print(if_symmetric(A))

#returns the transpose of matrix A
def transpose(A):
    
    A_transpose = [[0 for i in A] for j in A[0]]
    for i in range(len(A)):
        for j in range(len(A[0])):
            A_transpose[j][i] = A[i][j]
    return A_transpose

#A = [[2,3,0]]
#print(transpose(A))
#print(np.transpose(A))

#returns A+B
def add(A, B):
    if len(A)!= len(B) or len(A[0]) != len(B[0]):
        exit('addition not successful because of matrix size error')
    result = []
    for i in range(len(A)):
        result.append([A[i][k]+B[i][k] for k in range(len(A[0]))])
    return result

#returns A-B
def subtract(A, B):
    if len(A)!= len(B) or len(A[0]) != len(B[0]):
        exit('subtration not successful because of matrix size error')
    result = []
    for i in range(len(A)):
        result.append([A[i][k]-B[i][k] for k in range(len(A[0]))])
    return result

#returns A-B
def dot_product(A, B):
    if len(A[0])!= len(B):
        exit('dot product not successful because of matrix size error')
    result = []
    result = [[0 for col in B[0]] for row in A]
    for i in range(len(result)):
        for j in range(len(result[0])):
            result[i][j] = sum([A[i][k] * B[k][j] for k in range(len(B))])
    return result

#returns scalar*A
def scalar_product(scalar, A):    
    result = [[0 for col in range (len(A))]for row in range (len(A[0]))]
    for i in range (len(A)):
        for j in range (len(A[0])):
            result[i][j] = scalar*A[i][j]
    return result

#A = [[2,1]]
#print(scalar_product(4, A))

#returns a random symmetric positive definite square matrix of size n
def random_symmetric_positive_definite_matrix(n):
    A = np.random.randint(-10,10, size=(n,n))
    return dot_product(A,transpose(A))

#returns a random symmetric positive definite vector of size n
def random_vector(n):
    return np.random.randint(-10,10, size=(n))

