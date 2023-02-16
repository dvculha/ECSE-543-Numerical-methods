from matrices import *
import numpy as np
import math

# constants in SI units
E = 0.2
R = 512
IsA = 0.8 * 10**(-6)
IsB = 1.1 * 10**(-6)
Vt = 0.025
tolerance = 10**(-15)

def newton_raphson (V1, V2, k):

    #initial residuals
    f1 = V1 - E + R * IsA * (math.exp((V1 - V2) / Vt) - 1)
    f2 = IsA * ((math.exp((V1 - V2) / Vt) - 1)) - IsB * (math.exp(V2 / Vt) - 1)   

    print ('\nIteration #', k)

    #Partial Derivatives
    d_f1_V1 = 1.0 + (R * IsA / Vt) * (math.exp((V1 - V2) / Vt))
    d_f1_V2 = -1.0 * (R * IsA / Vt) * (math.exp((V1 - V2) / Vt))
    d_f2_V1 = (IsA / Vt) * (math.exp((V1 - V2) / Vt))
    d_f2_V2 = -1.0 * (IsA / Vt) * (math.exp((V1 - V2) / Vt)) - (IsB / Vt) * (math.exp(V2 / Vt))    

    # Jacobian Matrix
    J = [[d_f1_V1, d_f1_V2], [d_f2_V1, d_f2_V2]]

    V = [V1, V2]
    f = [f1, f2]

    #Update voltages
    #V_new = (-f+JV)/J = -f*Jinv+V
    V = add(scalar_product(-1, dot_product(inverse(J), f)), V)
    V1 = V[0][0]
    V2 = V[1][0]
    print ('V1 = ', V1)
    print ('V2 = ', V2)

    #Update residuals
    f1 = V1 - E + R * IsA * (math.exp((V1 - V2) / Vt) - 1)
    f2 = IsA * ((math.exp((V1 - V2) / Vt) - 1)) - IsB * (math.exp(V2 / Vt) - 1)   
    
    print ('f1 = ', f1)
    print ('f2 = ', f2)

    return V1, V2, f1, f2


if __name__ == "__main__":
    # initial guesses
    V1 = 0
    V2 = 0
    #initial residuals
    f1_1 = IsA * (math.exp((V1 - V2) / Vt) - 1) - IsB * (math.exp(V2 / Vt) - 1)
    f2_1 = V1- E + R * IsB * (math.exp(V2 / Vt) - 1)

    print ('\nInitial guesses')
    print('V1 = ', V1)
    print('V2 = ', V2)
    print ('\nInitial residual')
    print('f1 = ', f1_1)
    print('f2 = ', f2_1)
    f1=1
    f2=1

    #counter
    k = 1
    while (abs(f1) >=tolerance or abs(f2)>=tolerance):
        V1, V2, f1, f2 = newton_raphson(V1, V2,k)
        k = k + 1
    print ('\n')
