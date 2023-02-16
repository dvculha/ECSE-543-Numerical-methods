import math
import numpy as np
import time
import random
from matrices import *

def choleski(A, b, half_bandwidth=None):
    if isinstance(b[0],list):
        b = [x for r in b for x in r]
    # error flag: exit if A not symmetric
    if not (is_symmetric(A)):
        exit ("input matrix is not symmetric")

    n = len(A)

    # A is overwritten by L and b is overwritten by y
    for j in range(n - 1):
        # error flag: exit if A not positive definite
        if A[j][j] <= 0:
            exit ("input matrix is not positive definite")
        A[j][j] = math.sqrt(A[j][j])
        b[j] = b[j] / A[j][j]

        # regular approach, if half_bandwidth = none
        if (half_bandwidth == None):
            for i in range(j+1, n):
                A[i][j] = A[i][j] / A[j][j]
                b[i] = b[i] - A[i][j] * b[j]
                for k in range(j + 1, i + 1):
                    A[i][k] = A[i][k] - A[i][j] * A[k][j]

        # sparsity approach
        else:
            for i in range(j + 1, n):
                if i > j + half_bandwidth:
                    break
                
                A[i][j] = A[i][j] / A[j][j]
                b[i] = b[i] - A[i][j] * b[j]
                
                for k in range(j + 1, i + 1):
                    if (k > j + half_bandwidth):
                        break
                    A[i][k] = A[i][k] - A[i][j] * A[k][j]

    # Back substitution solving
    x = [0.0 for i in range(n)]
    for i in range(n - 1, -1, -1):
        x[i] = (b[i] - sum([A[j][i] * x[j] for j in range(i + 1, n)])) / A[i][i]
    return x

