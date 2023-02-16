import math
from matrices import dot_product, subtract, transpose, add, scalar_product
from choleski import choleski
import numpy as np
from finite_diff import *

free_node = 19
fixed_node = 15
inner_voltage = 110
outer_voltage = 0
num_nodes_x = 6
num_nodes_y = 6

def generate_A():
    
    mesh = initial_mesh(0.02)

    S_free = [[-4 if i == j else 0 for i in range(free_node)] for j in range(free_node)]

    inner_nodes = [0,1,2,3,5,6,7,8,10,11,12,13,15,16]
    neumann_nodes_x = [4,9,14]
    neumann_nodes_y = [17,18]

    for k in range (0, free_node):
        if (k in [0,1,2,3,5,6,7,8,10,11,12,13,15,17]):
            S_free[k][k+1]=1

        if (k in [1,2,3,6,7,8,11,12,13,16,18]):
            S_free[k][k-1]=1

        if (k <= 11):
            S_free[k][k+5]=1

        if (k >= 5 and k<=16):
            S_free[k][k-5]=1

        if (k in neumann_nodes_x):
            S_free[k][k-1]=2

        if (k in [15,16]):
            S_free[k][k+2]=1

        if (k in neumann_nodes_y):
            S_free[k][k-2]=2
    
    A = S_free

    #for r in A:
    #    print(r)

    return A

#A = generate_A()

def generate_b():
    S_prescribed = [[0 for i in range(fixed_node)] for j in range(free_node)]
    
    for k in range (0, fixed_node):

        if (k <=5):
            S_prescribed[k-1][k] =1
            #S_prescribed[k][k+1] = 1

        if (k >=6 and k<=9):
            S_prescribed[5*(k-6)][k]=1

        if (k >=10 and k<=12):
            S_prescribed[k+2][k]=1
        
        if k==10:
            S_prescribed [k+6][k]=1

        if k>=13:
            S_prescribed [k+4][k]=1

    voltages_prescribed = [0,0,0,0,0,0,0,0,0,0,110,110,110,0,110]
    voltages_prescribed = [i * -1 for i in voltages_prescribed]

    #print (voltages_prescribed)

    b = np.dot(S_prescribed,voltages_prescribed).tolist()

    #for r in S_prescribed:
    #    print(r)
    print(b)

    b = np.array(b).reshape(len(b),1)

    return b

#A = generate_b()

# converts matrix A to a positive definite matrix
def convert_to_positive_definite(A, b):
    determinant = np.linalg.det(A)
    #print ('Initial determinant = ', determinant)
    positive_definite = True
    A_new = A
    b_new = b
    if determinant <=0:
        positive_definite = False  
    if positive_definite == False:
        A_transpose = transpose(A)
        A_new = dot_product(A_transpose, A)
        b_new = dot_product(A_transpose, b)
        determinant = np.linalg.det(A_new)
        #print ('New determinant = ', determinant)
    return A_new, b_new

def two_norm(vector):  
    if isinstance(vector[0],list):
        vector = [i[0] for i in vector]
    two_norm = 0
    for entry in vector:
        entry = abs(entry)
        two_norm = two_norm + entry**2
    two_norm = math.sqrt(two_norm)
    return two_norm

def inf_norm(vector):
    if isinstance(vector[0],list):
        vector = [i[0] for i in vector]
    inf_norm = 0
    for entry in vector:
        entry = abs(entry)
        if entry > inf_norm:
            inf_norm = entry
    return inf_norm

def conjugate_gradient(A, b):
    # Guess x
    x = [[0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0]]
    
    # Set r and p
    # r = b - Ax
    r = subtract(b, dot_product(A, x))
    # p = r
    p = r

    # for k = 0,1,2......,number of free nodes
    for k in range(len(x)):
        iteration = k+1
        #print ('For Iteration #',iteration,": Two norm = ",two_norm(r)," and Infinity norm = ",inf_norm(r))
        #print ('Two norm = ',two_norm(r))
        #print ('Infinity norm = ',inf_norm(r))
        #print ("r = ", r)

        # alpha_k = p_k * r_k / (p_k_transpose * A_k * p_k)
        alpha = dot_product(transpose(p),r)[0][0]/dot_product(transpose(p), dot_product(A,p))[0][0]
        
        # x_(k+1) = x_k + alpha_k * p_k
        x = add (x,np.multiply(alpha,p))

        # r_k+1 = b - A * x_(k+1)
        r = subtract(b, dot_product(A,x))

        #beta_k = - p_k_transpose * A * r_(k+1) / (p_k_transpose * A * p_k)
        beta = -1*dot_product(transpose(p), dot_product(A,r))[0][0]/dot_product(transpose(p),dot_product(A,p))[0][0]
        
        # p_(k+1) = r_(k+1) + beta_k * p_k
        #print ('p=',p)
        #print ('beta=',b)
        p = add (r,np.multiply(beta,p))

    return x

if __name__ == '__main__':

    A = generate_A()
    b = generate_b()

    # 3a - Matrix Conversion
    print('\nInitial matrices: \n')
    print("A= ")
    for r in A:
        print(r)
    print ("\nb= ", b)
    print ('\nInitial determinant =', np.linalg.det(A))
    #x = choleski(A, b)

    print('\nTransformed matrices: \n')
    A0, b0 = convert_to_positive_definite(A, b)
    print("A= ")
    for r in A0:
        print(r)
    print ("\nb= ", b0)
    print ('\nNew determinant =', np.linalg.det(A0))
    #x = choleski(A, b)

    # 3b - Choleski
    A1, b1 = convert_to_positive_definite(A, b)
    #print("A= ")
    #for r in A:
    #    print(r)
    #print ("\nb= ", b) 
    print ("\nWhen solving with Choleski\n")
    x1 = choleski(A1, b1)
    print ('x = ', x1)
    print ('\n')

    # 3b - Conjugate Gradient
    A2,b2 = convert_to_positive_definite(A,b)
    #print('A = ')
    #for r in A1:
    #    print(r)
    #print ('\nb =', b1)
    print("\nWhen solving with Conjugate Gradient\n")
    x2 = conjugate_gradient(A2, b2)
    print ('\n')
    x2 = (format([i[0] for i in x2]))
    print ("x =", x2)
    print ('\n')

'''

S_free = [
    [-4, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
    [1, -4, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
    [0, 1, -4, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
    [0, 0, 1, -4, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
    [0, 0, 0, 2, -4, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
    [1, 0, 0, 0, 0, -4, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], 
    [0, 1, 0, 0, 0, 1, -4, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0], 
    [0, 0, 1, 0, 0, 0, 1, -4, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0], 
    [0, 0, 0, 1, 0, 0, 0, 1, -4, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0], 
    [0, 0, 0, 0, 1, 0, 0, 0, 2, -4, 0, 0, 0, 0, 1, 0, 0, 0, 0], 
    [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -4, 1, 0, 0, 0, 1, 0, 0, 0], 
    [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, -4, 1, 0, 0, 0, 1, 0, 0], 
    [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, -4, 1, 0, 0, 0, 0, 0], 
    [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, -4, 1, 0, 0, 0, 0], 
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, -4, 0, 0, 0, 0], 
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -4, 1, 1, 0], 
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, -4, 0, 1], 
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, -4, 1], 
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, -4]]



S_prescribed = [
    [0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]]

voltages_prescribed = [0,0,0,0,0,0,0,0,0,0,110,110,110,0,110]

b = [[0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [-110], [-110], [-110], [0], [-110], [0], [-110]]

'''