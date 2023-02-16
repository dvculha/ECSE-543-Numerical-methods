import math
import csv

file = open('problem.txt','w')

inner_voltage = 110
outer_voltage = 0
h = 0.02
num_nodes_x = 6
num_nodes_y = 6

#initial node numbers mesh
mesh = [[(i + j * num_nodes_x+1) for i in range(num_nodes_x)] for j in range(num_nodes_y)]
#print (mesh)

# input node list
for j in range(num_nodes_x):
    for i in range(num_nodes_y):
        if (mesh[j][i] <= 34):
            node_number = mesh[j][i]
            x = i*h
            y = j*h
            file.write('%d %.2f %.2f\n'%(node_number, x, y))
            print('%d %.2f %.2f'%(node_number, x, y))
            
file.write('\n')
print ('\n')

# input element list
for x in range(num_nodes_x - 1):
    for y in range(num_nodes_y - 1):
        if (mesh[x][y] not in [35,36]) and (mesh[x+1][y] not in [35,36]) and (mesh[x][y+1] not in [35,36]) and (mesh[x+1][y+1] not in [35,36]):
            i = mesh[x][y]
            j = mesh[x+1][y]
            k = mesh[x][y+1]
            m = mesh[x+1][y+1]
            source = 0.0
            file.write('%d %d %d %d\n'%(i, j, k, source)) 
            print ('%d %d %d %d'%(i, j, k, source)) 
            file.write('%d %d %d %d\n'%(k, m, j, source)) 
            print('%d %d %d %d'%(k, m, j, source)) 
file.write('\n')
print('\n')


# input fixed potentials 
# dirichlet conditions
for i in range (num_nodes_x):
    for j in range (num_nodes_y):
        node_number = mesh[i][j]
        if (node_number <= 34):

            #dirichlet boundaries

            #inside inner rectangle
            if (mesh[i][j] in [28, 29, 30, 34]):
                voltage = inner_voltage
                file.write('%d %.2f\n'%(node_number,voltage))
                print('%d %.2f'%(node_number,voltage))

            # on side of outer rectangle
            if (node_number % 6 == 1) or (node_number <= 6):
                voltage = outer_voltage
                file.write('%d %.2f\n'%(node_number,voltage))
                print ('%d %.2f'%(node_number,voltage))

file.close()
'''

            #neumann boundaries
            if (node_number % 6 == 0 and node_number <= 24):
                delta_voltage_y = (outer_voltage-inner_voltage) / (4) # /(outer_width-inner_width)/h
                voltage = outer_voltage - delta_voltage_y * (node_number/6 - 1)
                file.write('%d %.2f\n'%(node_number,voltage))
                print('%d %.2f'%(node_number,voltage))
                
            if (node_number > 31 and node_number < 35):
                delta_voltage_x = (outer_voltage-inner_voltage) / (3)
                voltage = outer_voltage - delta_voltage_x * (node_number - 31)
                file.write('%d %.2f\n'%(node_number,voltage))
                print('%d %.2f'%(node_number,voltage))

'''