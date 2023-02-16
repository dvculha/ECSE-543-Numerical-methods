import math
import numpy as np
import time
import random
from matrices import *

def choleski(A, b, half_bandwidth):
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
        x[i] = round((b[i] - sum([A[j][i] * x[j] for j in range(i + 1, n)])) / A[i][i], 10)
    return x

if __name__ == "__main__":

    for n in range(2, 11):
        print("\nFor N =", n)
        
        if (n==2):
            A = [[2,0],[0,8]]
            b = [-4,-16]
            #x = [-2,-2]

        if (n==3):
            A = [[8,0,2],[0,9,1],[2,1,2]]
            b = [-12,-20,-8]
            #x = [-1,-2,-2]

        if (n==4):
            A = [[9,4,-2,1],[4,6,4,-2],[-2,4,8,0],[1,-2,0,10]]
            b = [2,-4,-8,0]
            #x = [0,0,-1,0]

        if (n==5):
            A = [[11,3,-10,4,0],[3,14,-2,-1,8],[-10,-2,12,-2,2],[4,-1,-2,10,-7],[0,8,2,-7,14]]
            b = [-11,-11,10,-2,-3]
            #x = [-2,-1,-1,1,1]

        if (n==6):
            A = [[11,-4,-1,-1,-6,2],[-4,12,0,-6,-2,-2],[-1,0,15,-3,-2,2],[-1,-6,-3,5,4,0],[-6,-2,-2,4,7,0],[2,-2,2,0,0,18]]
            b = [-3,4,-31,7,8,-38]
            #x = [2,2,-1,2,2,-2]
        
        if (n==7):
            A = [[8,4,1,1,-2,-2,-5],[4,15,-9,5,-2,-10,0],[1,-9,22,2,-4,-1,6],[1,5,2,10,-3,-4,9],[-2,-2,-4,-3,10,2,-2],[-2,-10,-1,-4,2,18,-9],[-5,0,6,9,-2,-9,19]]
            b = [3,14,19,27,7,-48,51]
            #x = [2,-1,0,1,2,-2,2]

        if (n==8):
            A = [[22,13,0,-1,4,-10,-6,-1],[13,15,-6,-2,-1,0,0,3],[0,-6,17,-10,0,-7,-3,-3],[-1,-2,-10,16,5,0,-6,1],[4,-1,0,5,19,4,-1,-5],[-10,0,-7,0,4,25,8,-4],[-6,0,-3,-6,-1,8,11,-1],[-1,3,-3,1,-5,-4,-1,7]]
            b = [-45,-26,-23,31,-4,10,3,13]
            #x = [-1,-1,0,2,0,0,1,2]
        
        if (n==9):
            A = [[19,-2,-1,4,-2,2,-9,-16,7],[-2,29,-2,-12,-8,14,-12,3,1],[-1,-2,15,0,1,2,8,-6,-9],[4,-12,0,13,3,-6,2,-3,0],[-2,-8,1,3,18,-5,9,2,-7],[2,14,2,-6,-5,15,-10,2,-1],[-9,-12,8,2,9,-10,23,3,-11],[-16,3,-6,-3,2,2,3,21,-4],[7,1,-9,0,-7,-1,-11,-4,11]]
            b = [68,12,-9,7,-47,16,-69,-62,48]
            #x = [2,0,2,0,-1,-1,-2,0,2]

        if (n==10):
            A=[
            [159,53,21,37,9,-34,-13,-14,-28,15],[53,103,-16,-5,-3,44,-4,29,38,21],[21,-16,95,-5,45,-41,17,14,11,-4],
            [37,-5,-5,122,71,2,-31,12,-56,1],[9,-3,45,71,115,-7,-1,17,-5,1],
            [-34,44,-41,2,-7,85,-19,36,57,12],[-13,-4,17,-31,-1,-19,81,-26,13,23],
            [-14,29,14,12,17,36,-26,74,33,11],[-28,38,11,-56,-5,57,13,33,92,11],[15,21,-4,1,1,12,23,11,11,32]]
            
            b = [-825,-766,214,472,664,-33,-359,129,-176,-339]
            #x = [-4,-5,1,1,5,1,-4,0,0,-3]
        
        print("A =",A)
        print("b =",b)
        x = choleski(A, b, None)
        print("x =",x)

        