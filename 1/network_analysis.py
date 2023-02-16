import math
import numpy as np
import time
import os
from choleski import choleski
from matrices import *


#Textfile to network branch
#J = 1st line of file
#R = 2nd line of file
#E = 3rd line of file
#A = From first row to 4th to last row. or each new row in A, add a new row in the file too
def txt_to_network_branch(filename):
    with open(filename) as f:
        X = []
        for line in f:
            if line.strip() == '':
                continue
            X.append([float(x) for x in line.split(',')])
        J = X[0] 
        R = X[1] 
        E = X[2]
        A = X[3:]
        return A, J, R, E

# Solves for the network voltage, given A, J, R, E
# y = 1/R
# (A*y*A_transpose)*V = A*(J-y*E) 
def solve_network(A, J, R, E, half_bandwidth):

    y = [[0 for r in R] for r in R]
    for i in range(len(y)):
        y[i][i] = 1.0/R[i]

    #(A*y*A_transpose)
    left = dot_product((dot_product(A,y)),transpose(A))
    
    #A*(J-y*E)
    right = dot_product(A, subtract([[j] for j in J], dot_product(y, [[e] for e in E])))
    right = [i[0] for i in right]
    V = choleski(left, right, half_bandwidth)
    return V

if __name__ == '__main__':

    for i in range (1,6):
        print ("\nFor test circuit",i)

        if (i==1):
            filename = ('network_branch_1.txt')
        if (i==2):
            filename = ('network_branch_2.txt')
        if (i==3):
            filename = ('network_branch_3.txt')
        if (i==4):
            filename = ('network_branch_4.txt')
        if (i==5):
            filename = ('network_branch_5.txt')

        A, J, R, E = txt_to_network_branch(filename)
        print("A = ", A)
        print("J = ", J)
        print("R = ", R)
        print("E = ", E)
        v = solve_network(A, J, R, E, None)
        print ('Vout = ',v)
        
    