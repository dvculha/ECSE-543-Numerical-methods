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


# transforms col vectors to right format.
# [x,y,z] into [[x],[y],[z]]
def transform(A):
    try:
        cols = len(A[0])
        return A
    except TypeError:
        result = [[0] for a in range(len(A))]
        for i in range(0, len(A)):
            result[i][0] = A[i]
        return result
#A = [2,3,4]
#print(A)
#print(transform(A))

#returns the transpose of matrix A
def transpose(A):
    A = transform (A)
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
    A = transform(A)
    B = transform(B)
    if len(A)!= len(B) or len(A[0]) != len(B[0]):
        exit('addition not successful because of matrix size error')
    result = []
    for i in range(len(A)):
        result.append([A[i][k]+B[i][k] for k in range(len(A[0]))])
    return result
#A = [2,3,0]
#B = [1,3,0]
#print(add(A,B))

#returns A-B
def subtract(A, B):
    A = transform(A)
    B = transform(B)
    if len(A)!= len(B) or len(A[0]) != len(B[0]):
        exit('subtration not successful because of matrix size error')
    result = []
    for i in range(len(A)):
        result.append([A[i][k]-B[i][k] for k in range(len(A[0]))])
    return result
#A = [[2,3,0], [2,2,1]]
#B = [[1,3,0], [5,1,1]]
#print(subtract(A,B))

#returns A*B
def dot_product(A, B):
    A = transform(A)
    B = transform(B)
    if len(A[0])!= len(B):
        exit('dot product not successful because of matrix size error')
    result = []
    result = [[0 for i in range(len(B[0]))]for k in range(len(A))]
    for i in range(len(A)):
        for j in range(len(B[0])):
            for k in range(len(A[0])):
                result[i][j] += A[i][k]*B[k][j]    
    return result
#A = [[1,3], [0,1]]
#B = [3,1]
#print(dot_product(A,B))

#returns scalar*A
def scalar_product(scalar, A):    
    A = transform(A)
    result = [[0 for i in range(len(A[0]))]for k in range(len(A))]
    for i in range(len(A)):
        for j in range(len(A[0])):
            result[i][j] = scalar*A[i][j]
    return result
#A = [2,6]
#print(scalar_product(5,A))

#returns determinant of a 2x2 matrix
def determinant(A):
    det = A[0][0] * A[1][1] - A[0][1] * A[1][0]
    return det

#returns product inverse of A
def inverse (A):
    det_A = determinant(A)
    inv_A = [[None for x in range (len(A))] for y in range(len(A))]
    inv_A[0][0] = A[1][1] / det_A
    inv_A[0][1] = -1 * A[0][1] / det_A
    inv_A[1][0] = -1 * A[1][0] / det_A  
    inv_A[1][1] = A[0][0] / det_A
    return inv_A  
#A = [[2,6],[3,1]]
#print(determinant(A))
#print (inverse(A))

#returns a random symmetric positive definite square matrix of size n
def random_symmetric_positive_definite_matrix(n):
    A = np.random.randint(-10,10, size=(n,n))
    return dot_product(A,transpose(A))

#returns a random symmetric positive definite vector of size n
def random_vector(n):
    return np.random.randint(-10,10, size=(n))
