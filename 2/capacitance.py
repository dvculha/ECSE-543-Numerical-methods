import scipy
import numpy as np
from matrices import dot_product, transpose

nodes_y = 6
nodes_x = 6
e0 = 8.854187817 * 10 ** (-12)
V = 110.0
#potentials = [0, 0, 0, 0, 0, 0], [0, 7.3254, 14.2600, 20.1675, 24.3243, 27.5000], [0, 15.0416, 29.5472, 42.0858, 49.6295, 55.0000],[0, 23.2939, 46.8015, 68.9989, 77.1081, 79.8041], [0, 31.3325, 65.3660, 110.0000, 110.0000, 110.0000], [0, 36.6700, 73.3300, 110.0000, 110, 110]
# potentials from finite element
potentials = [0, 0, 0, 0, 0, 0], [0, 7.0186, 13.6519, 19.1107, 22.2643, 23.2569], [0, 14.4223, 28.4785, 40.5265, 46.6897, 48.4989],[0, 22.1921, 45.3132, 67.8272, 75.4690, 77.3592], [0, 29.0330, 62.7550, 110.0000, 110.0000, 110.0000], [0, 31.1849, 66.6737, 110.0000, 110, 110]
#potentials from finite difference
potentials = [0, 0, 0, 0, 0, 0], [0, 7.01855349069135, 13.651930409784672, 19.11068359610681, 22.26430271342167, 23.256864494986434], [0, 14.422287599181919, 28.478476428333863, 40.52649982978976, 46.68966572649273, 48.498852862499184],[0, 22.192123760699417, 45.3131901052967, 67.82717379330018, 75.46901214903302, 77.3592182954247], [0, 29.0330096583123, 62.754978837895855, 110.0000, 110.0000, 110.0000], [0, 31.184931218436702, 66.673721421056, 110.0000, 110, 110]


S_con =  [[1, -0.5, 0, -0.5],[-0.5, 1, -0.5, 0], [0, -0.5, 1, -0.5], [-0.5, 0, -0.5, 1]]

W = 0.0
U_con = [0 for i in range (0,4)]
#print (U_con)
# nodes_y - 1 and nodes_x -1 because index nodes_y and x already get covered
for j in range (0, nodes_y - 1):
    for i in range(0, nodes_x - 1):
        U_con[0] = potentials[i][j]
        U_con[1] = potentials[i+1][j]
        U_con[2] = potentials[i+1][j+1]
        U_con[3] = potentials[i][j+1]
        #print ('U_con = ',U_con)
        #U_con_transpose = transpose([U_con])  
        U_con_transpose = np.transpose(U_con)  
        #print ('U_con_transpose = ', U_con)  
        W_new = np.dot(np.dot(U_con,S_con),U_con_transpose)
        #print ('potential = ', potential)
        
        W = W + W_new
        #print ('W = ', W)
W = 0.5*W
#print ('W = ', W)

C_unit_length = 2 * e0 * W / V**2 * 4
print ('Capacitance per unit length is', (C_unit_length), ' F/m')
print ('\n')



