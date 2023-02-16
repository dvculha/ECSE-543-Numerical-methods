import numpy as np
import time
import math
from network_analysis import solve_network


def generate_network(N, R_test, E_test):

    # number of rows, columns, nodes and branches
    row = N+1
    cols = N+1
    nodes = row*cols
    branches = 2*row*cols-row-cols

    #current of each branch
    J = [0]*(branches+1) 

    #resistance of each branch
    R = [0]*(branches+1)
    for i in range (branches):
        R[i]= 10000.0 #resistance of each branch

    #find E
    E = [0]*(branches+1)

    # generate incidence matrix
    A = [[0]*(branches+1) for i in range(nodes)]

    for i in range(0, nodes):
        level = int (i/cols)
        offset = i % cols
        
        # x branches
        it_x = level * (cols - 1) + offset
        rem_x = it_x - 1

        # y branches
        it_y = nodes - row + i 
        rem_y = it_y - cols

        #middle nodes
        # if not top nodes
        if (i >= cols):
            A[i][rem_y] = -1
        # if not leftmost branches
        if (offset != 0):
            A[i][rem_x] = -1
        # if not rightmost branches
        if (offset != (cols - 1)):
            A[i][it_x] = +1
        # if not bottom nodes
        if (i <= (nodes - cols - 1)):
            A[i][it_y] = +1
        
    # test branch  
    A[0][branches] = -1
    A[nodes-1][branches] = 1
    R [branches] = R_test
    E[branches]= E_test

    A = A[:-1]
    #print ('A = ',A)
    #print ('J = ',J)
    #print ('R = ',R)
    #print ('E = ',E)
    return A, J, R, E

if __name__ == '__main__':

    for N in range(2,16):
        print ('\nFor a {}x{} mesh,' .format(N,N))

        R_test = 10000.0
        E_test = 10000.0

        A, J, R, E = generate_network(N, R_test, E_test)

        print ('When normal structure is used,')
        half_bandwidth = None
        #start timer
        time_start = time.perf_counter() 
        #solve for node voltages
        v = solve_network(A, J, R, E, half_bandwidth)
        #end timer
        time_end = time.perf_counter()
        delta_time = time_end - time_start
        #voltage at node 0
        V_node_0 = float(v[0])
        
        #calculate R_eq with voltage divider formula
        #R_eq / (R_eq + R_test) = V_node_0 / V_test
        #R_eq *V_test = V_node_0*(R_eq + R_test)
        #R_eq*(V_test-V_node_0) = V_node_0 * R_test
        R_eq = V_node_0*R_test/(E_test - V_node_0)  

        print ('Equivalent resistance =', R_eq, "ohms")
        print ('Computation time =', delta_time, "seconds")

        print ('When banded structure is used,')
        half_bandwidth = N+1
        time_start = time.perf_counter() 
        v = solve_network(A, J, R, E, half_bandwidth)
        time_end = time.perf_counter() 
        delta_time = time_end - time_start
        #print ('v=' ,v)

        V_node_0 = float(v[0])
        
        #calculate R_eq with voltage divider formula
        #R_eq / (R_eq + R_test) = V_node_0 / V_test
        #R_eq *V_test = V_node_0*(R_eq + R_test)
        #R_eq*(V_test-V_node_0) = V_node_0 * R_test
        R_eq = V_node_0*R_test/(E_test - V_node_0) 

        print ('Equivalent resistance =', R_eq, "ohms")
        print ('Computation time =', delta_time, "seconds")
     

