import math

# Constant parameters
# Divided by 2 because only using the 4th quadrant
outer_voltage = 0.0
inner_voltage = 110.0
inner_width = 0.04/2
inner_height = 0.08/2
outer_width = 0.2/2
outer_height = 0.2/2
residual_tolerance = 0.00001

def transform_target (target_x, target_y):

    if (target_x > outer_width):
        target_x = outer_width * 2 - target_x 

    if (target_x < outer_width):
        target_x = target_x

    if (target_y < outer_height):
        target_y = target_y

    if (target_y > outer_height):
        target_y = outer_height * 2 - target_y


    #target_y = outer_height - target_y
    target_y=round(target_y,5)
    target_x=round(target_x,5)

    return target_x, target_y

#print(transform_target (0.06,0.04))
#target becomes 0.04,0.06 which is correct in my coordinates :D


# Boundaries and inner rectangle, solving in the 4th quadrant only 
def initial_mesh (h):

    x = int(outer_width / h)
    y = int(outer_height / h)

    # max node indexes in inner and outer rectangles
    x_outer = int(outer_width / h)
    y_outer = int(outer_height / h)
    x_inner = int (inner_width/h)
    y_inner = int (inner_height/h)
    
    # number_of_nodes_x = outer_width / h + 1
    # number of nodes_y = outer_height / h + 1

    # Dirichlet boundary conditions

    mesh = [[inner_voltage if j <= x_inner and i <= y_inner else outer_voltage if j==x_outer and i==y_outer else 0.0 for i in range(0, x_outer+1)] for j in range(0, y_outer+1)]

    #print ('Initial mesh with Dirichlet = ' ,mesh)
    # Neuman n Boundary conditions

    # Voltage changes between boundary nodes per node
    delta_voltage_x = (outer_voltage-inner_voltage) / (x_outer - x_inner) # /(outer_width-inner_width)/h
    delta_voltage_y = (outer_voltage-inner_voltage) / (y_outer - y_inner)

    # when i = x_outer it is equal to outer_voltage
    # when i = x_inner it is equal to inner_voltage
    for i in range(x_inner+1, x_outer):
        mesh[i][0] = mesh[i-1][0] + delta_voltage_x
    
    for j in range(y_inner+1, y_outer):
        mesh[0][j] = mesh[0][j-1] + delta_voltage_y
    
    #print ('Initial mesh with Neuman = ' ,mesh)

    return mesh

def residual_tolerable(mesh, h):

    x_outer = int(outer_width / h)
    y_outer = int(outer_height / h)
    x_inner = int (inner_width/h)
    y_inner = int (inner_height/h)

    max_residual = 0

    #no need to calculate residue for boundaries
    for i in range(1, x_outer):
        for j in range(1, y_outer):
            if (i > x_inner or j > y_inner):
                current_residual = math.fabs(mesh[i-1][j] + mesh[i+1][j] + mesh[i][j-1] + mesh[i][j+1] - 4 * mesh[i][j])
                max_residual = max(max_residual,current_residual)
    
    if (max_residual < residual_tolerance):
        return True
    else: 
        return False
  

def SOR(mesh, h, w):

    x_outer = int(outer_width / h)
    y_outer = int(outer_height / h)
    x_inner = int (inner_width/h)
    y_inner = int (inner_height/h)

    # x=0, y=0s are boundaries, inner_width x inner_height is filled in too
    for j in range (1, y_outer):
        for i in range(1, x_outer):  
            if (i > x_inner) or j > int(y_inner):
                mesh[i][j] = (1 - w) * mesh[i][j] + (w/4) * (mesh[i-1][j] + mesh[i+1][j] + mesh[i][j-1] + mesh[i][j+1])    
    return mesh

def jacobi (mesh, h):

    x_outer = int(outer_width / h)
    y_outer = int(outer_height / h)
    x_inner = int (inner_width/h)
    y_inner = int (inner_height/h)

    for i in range (1, x_outer):    
        for j in range (1, y_outer):
            if (i > (x_inner) or j > (y_inner)):
                mesh[i][j] = (1/4) * (mesh[i-1][j] + mesh[i+1][j] + mesh[i][j-1] + mesh[i][j+1])
    return mesh

def SOR_solve (x_target_pos, y_target_pos, h, w):

    x_target_node = int(x_target_pos/h)
    y_target_node = int(y_target_pos/h)

    problem_mesh = initial_mesh(h)

    solution_mesh = SOR (problem_mesh, h, w)

    iteration = 1

    while (residual_tolerable(solution_mesh, h) == False):
        iteration = iteration + 1
        solution_mesh = SOR (solution_mesh, h, w)
        #print ("\nSolution mesh =",solution_mesh)
    
    #print (solution_mesh)
    
    return solution_mesh [x_target_node][y_target_node], iteration

def jacobi_solve (x_target_pos, y_target_pos, h):
    
    x_target_node = int(x_target_pos/h)
    y_target_node = int(y_target_pos/h)

    problem_mesh = initial_mesh(h)

    solution_mesh = jacobi (problem_mesh, h)

    iteration = 1

    while (residual_tolerable(solution_mesh, h) != True):
        iteration = iteration + 1
        solution_mesh = jacobi (solution_mesh, h)
        
    return solution_mesh [x_target_node][y_target_node] , iteration

#NON UNIFORM CASE

x_address = [0, 0.01,0.02, 0.03,0.04, 0.05,0.06, 0.07,0.08, 0.09,0.1]
y_address = [0, 0.01,0.02, 0.03,0.04, 0.05,0.06, 0.07,0.08, 0.09,0.1]

def initial_mesh_non_uniform ():

    # Dirichlet boundary conditions
    mesh = [[inner_voltage if j <= inner_width and i <= inner_height else outer_voltage if j==outer_width and i==outer_height else 0.0 for i in (x_address)] for j in (y_address)]
    #print ('Initial mesh with Dirichlet = ' ,mesh)
    
    # Neuman n Boundary conditions
    # Voltage changes between boundary nodes per node
    delta_voltage_x = (outer_voltage-inner_voltage) / (outer_width - inner_width) # /(outer_width-inner_width)/h
    delta_voltage_y = (outer_voltage-inner_voltage) / (outer_height - inner_height)

    # when i = x_outer it is equal to outer_voltage
    # when i = x_inner it is equal to inner_voltage
    for i in range(len(x_address)):
        if (x_address[i] > inner_width):
            mesh[i][0] = inner_voltage + delta_voltage_x * (x_address[i] - inner_width)

    for j in range(len(y_address)):
        if (y_address[j] > inner_height):
            mesh[0][j] = inner_voltage + delta_voltage_y * (y_address[j] - inner_height)

    return mesh


def residual_tolerable_non_uniform(mesh):

    max_residual = 0

    #no need to calculate residue for boundaries
    for i in range(1, len(x_address)-1):
        for j in range(1, len(y_address)-1):
            if (x_address[i] > inner_width or y_address[i] > inner_height):
                a1 = x_address[i] - x_address[i-1]
                a2 = x_address[i+1] - x_address[i]
                b1 = y_address[j+1] - y_address[j]
                b2 = y_address[j] - y_address[j-1]

                current_residual = math.fabs((mesh[i-1][j]/(a1 * (a1 + a2)) + mesh[i+1][j]/(a2 * (a1 + a2)) + mesh[i][j-1]/(b1 * (b1 + b2)) + mesh[i][j+1]/(b2 * (b1 + b2))) - mesh[i][j]*(1/(a1 * a2) + 1/(b1 * b2)))
                max_residual = max(max_residual, current_residual)
    
    if (max_residual < residual_tolerance):
        return True
    else: 
        return False

def SOR_non_uniform(mesh, w):

    problem_mesh = initial_mesh_non_uniform()

    for i in range (1, len(x_address)-1):    
        for j in range (1, len(x_address)-1):
            if (x_address[i] > inner_width or y_address[j] > inner_height):
            
                a1 = x_address[i] - x_address[i - 1]
                a2 = x_address[i + 1] - x_address[i]
                b1 = y_address[j] - y_address[j - 1]
                b2 = y_address[j + 1] - y_address[j]

                e = math.pow(a1,2) + math.pow(a2,2)
                f = math.pow(b1,2) + math.pow(b2,2)
                #mesh[i][j] = (1-w) * mesh[i][j] + w * (mesh[i-1][j]/(a1 * (a1 + a2)) + mesh[i+1][j] / (a2 * (a1 + a2)) + mesh[i][j-1] / (b1 * (b1 + b2)) + mesh[i][j+1] / (b2 * (b1 + b2)))
                mesh[i][j] = (1-w) * mesh[i][j] + w * (f/2/(e+f) * (mesh[i-1][j]+ mesh[i+1][j]) + e/2/(e+f) * (mesh[i][j-1] + mesh[i][j+1]))
    #print (mesh)

    return mesh

def SOR_non_uniform_solve (x_target_pos, y_target_pos, w):

    x_node = x_address.index(x_target_pos)
    y_node = y_address.index(y_target_pos)

    problem_mesh = initial_mesh_non_uniform()

    solution_mesh = SOR_non_uniform (problem_mesh, w)

    iteration = 1

    while (residual_tolerable_non_uniform(solution_mesh) != True):
        iteration = iteration + 1
        solution_mesh = SOR_non_uniform (solution_mesh, w)
    
    #print (solution_mesh)
    return solution_mesh [x_node][y_node] , iteration

if __name__ == "__main__":

    #since my program solves in quadrant 4 (lower right)
    target_x, target_y = 0.06, 0.04
    my_target_x, my_target_y = transform_target(target_x,target_y)

    ### part b
    h = 0.02
    
    print ('\n')
    print ("When solving with SOR with uniform spacing")
    print ('Constant spacing h is', h)
    for w in [1.0,1.1,1.2,1.25,1.3,1.35,1.4,1.5,1.6,1.7,1.8,1.9]:
        #w = x / 10
        print ("\nWhen w is ",w,":")
        solution, iteration = SOR_solve (my_target_x, my_target_y, h, w)
        print ("Voltage at ",target_x,",",target_y," is ", solution)
        print ("Iteration is" ,iteration)

    # w = 1.25 gives the least amount of iterations

    # part c
    print ('\n')
    w = 1.25
    print ("When solving with SOR with uniform spacing")
    h = 0.02
    for i in range (0,7):
        print ("\nWhen h is ",h,":")
        solution, iteration = SOR_solve (my_target_x, my_target_y, h, w)
        print ("Voltage at ",target_x,",",target_y," is ", solution)
        print ("Iteration is" ,iteration)
        h = h * 0.5

    # part d
    print ('\n')
    print ("When solving with Jacobi")
    h = 0.02
    for i in range (0,7):
        print ("\nWhen h is ",h,":")
        solution, iteration = jacobi_solve (my_target_x, my_target_y, h)
        print ("Voltage at ",target_x,",",target_y," is ", solution)
        print ("Iteration is" ,iteration)
        h = h * 0.5

    # part e

    print ('\n')
    print ("When solving with non-uniform spacing SOR")
    
    w = 1.25
    print ("\nWhen w is ",w,":")
    solution, iteration = SOR_non_uniform_solve (my_target_x, my_target_y, w)
    print ("Voltage at ",target_x,",",target_y," is ", solution)
    print ("Iteration is" ,iteration)