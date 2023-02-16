from sympy import *

def lagrange_interpolation (points_x, points_y):  
    param = symbols('B')
    n = len(points_x)
    yx = 0
    Ljx = [0.0 for j in range (n)]    
    a = [None for j in range (n)]
    # For j = 0,1,2,....,n   
    for j in range (n):
        Fjx = 1
        Fjxj = 1
        for r in range (n):
            if (r != j):
                # Fjx = PRODUCT (x-xr)
                Fjx = Fjx * (param - points_x[r])
                Fjxj = Fjxj * (points_x[j] - points_x[r])
        #Ljx = Fjx / Fjxj
        Ljx[j] = Fjx.expand()/Fjxj

    # a[j] = y(j)
    for j in range (n):
        a[j] = points_y[j]

    # yx = SUM (aj*Ljx)
    for j in range (n):
        yx = yx + a[j] * Ljx[j]             
    return yx 

def cubic_hermite(points_x,points_y):
    n = len(points_x)
    param = symbols('B')
    Ujx = []
    Vjx = []
    yx = 0
    # Lj(x)
    Ljx = [None for j in range (n)]
    #L'j(x)
    Ljx_prime = [None for j in range (n)]
    b = [None for j in range (n)] 
    a = [None for j in range (n)] 

    for j in range(n):
        Fjx = 1
        Fjxj = 1
        for r in range (n):
            if (r != j):
                # Fjx = PRODUCT (x-xr)
                Fjx = Fjx * (param - points_x[r])
                Fjxj = Fjxj * (points_x[j] - points_x[r])
        # To obtain Lj(x)
        # Lj(x) = Fjx / Fjxj
        Ljx = Fjx.expand() / Fjxj

        #To obtain L'j(x), diff Lj(x) / dx
        Ljx_prime = lambdify(param, diff(Ljx))

        U = (1 - 2 * Ljx_prime(points_x[j]) * (param - points_x[j]))*(Ljx**2)
        Ujx.append(U)
        V = (param - points_x[j])*(Ljx * Ljx)
        Vjx.append(V)

    # a[j] = y(j)
    for j in range (n):
        a[j] = points_y[j]

    # b[j] = y'(j)
    for j in range (n):
        if (j < (n-1)):    
            b[j] = (points_y[j + 1] - points_y[j]) / (points_x[j + 1] - points_x[j])
        elif (j == (n-1)):
            b[j] = points_y[j]/points_x[j]

    # y(x) = SUM (a*Ujx + b*Vjx)
    for j in range(n):
        yx = yx + a[j]*Ujx[j] + b[j]*Vjx[j]
    yx = expand(yx)
    return yx



if __name__ == "__main__":

    B = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9]
    H = [0.0, 14.7, 36.5, 71.7, 121.4, 197.4, 256.2, 348.7, 540.6, 1062.8, 2318.0, 4781.9, 8687.4, 13924.3, 22650.2]


    print ('\nPart A\n')
    points_x = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    print ('B = ', points_x)
    points_y = [0.0, 14.7, 36.5, 71.7, 121.4, 197.4]
    print ('H = ', points_y)
    a = lagrange_interpolation(points_x, points_y)
    print ('H = ',a)

    print ('\n')

    print ('Part B\n')
    points_x = [0.0, 1.3, 1.4, 1.7, 1.8, 1.9]
    print ('B = ', points_x)
    points_y = [0.0, 540.6, 1062.8, 8687.4, 13924.3, 22650.2]
    print ('H = ', points_y)
    b = lagrange_interpolation(points_x, points_y)
    print ('H = ',b)

    print ('\n')

    print ('Part C\n')
    points_x = [0.0, 1.3, 1.4, 1.7, 1.8, 1.9]
    print ('B = ', points_x)
    points_y = [0.0, 540.6, 1062.8, 8687.4, 13924.3, 22650.2]
    print ('H = ', points_y)
    c = cubic_hermite(points_x, points_y)
    print ('H = ',c)

    print ('\n')