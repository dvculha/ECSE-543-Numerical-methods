import math
import numpy as np

B = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9]
H = [0.0, 14.7, 36.5, 71.7, 121.4, 197.4, 256.2, 348.7, 540.6, 1062.8, 2318.0, 4781.9, 8687.4, 13924.3, 22650.2]

Area = 10**(-4)
Lc = 0.3
La = 0.005
N = 800
I = 10
perm0 = 4*math.pi*10**(-7)

Ra = La / (perm0 * Area)
#print (Ra)
MMF = N * I

tolerance = 10**(-6)

# B = flux / area
def get_B (flux):
    B = flux / Area
    return B
#print(get_B(0.3))

# B = mu * H
# mu = B / H
# H = B / mu
# H_prime = 1/mu
def get_H_and_H_prime(flux):
    Bc = flux/Area
    if Bc > B[-1]:
        mu = (B[-1] - B[-2]) / (H[-1] - H[-2])
        Hc = (Bc - B[-1]) / mu + H[-1]
        H_prime = 1/mu
        return Hc, H_prime
    
    for i in range(len(B)):
        if B[i] > Bc:
            mu = (B[i] - B[i-1]) / (H[i] - H[i-1])
            Hc = (Bc - B[i-1]) / mu + H[i-1]
            H_prime = 1/mu
            return Hc, H_prime
    
def get_f(flux):
    Hc, H_prime = get_H_and_H_prime(flux)
    f = Hc * Lc + Ra * flux - MMF
    return f

def get_f_prime(flux):
    Hc, H_prime = get_H_and_H_prime(flux)
    f_prime = H_prime * Lc / Area + Ra
    return f_prime

def newton_raphson(flux):
    print ('Initial flux = ', flux)
    i = 1
    while abs(get_f(flux)/get_f(0)) > tolerance:
        print ('\nIteration #', i)
        flux = flux - get_f(flux)/get_f_prime(flux)
        print ('Flux = ', flux)
        i += 1
    return flux

def successive_sub_xxx(flux):
    i = 1
    while abs(get_f(flux)/get_f(0)) >= tolerance:
        print ('\nIteration #', i)
        flux = get_f(flux)
        print ('Flux = ', flux)
        i += 1
    return flux

def get_Bc(flux):
    Hc = (MMF - Ra * flux) / Lc
    for i in range(len(H) - 1):
        if H[i] <= Hc < H[i + 1]:
            mu = ((B[i + 1] - B[i])/(H[i + 1] - H[i]))
            Bc = B[i] + ((Hc) - H[i]) * mu
            return Bc
        elif Hc > H[- 1]:
            mu = ((B[- 1] - B[- 2])/(H[- 1] - H[- 2]))
            Bc = B[- 1] + (Hc - H[- 1]) * mu
            return Bc

def successive_sub(flux):
    i = 1
    print ('Initial flux = ', flux)
    while abs(get_f(flux)/get_f(0)) >= tolerance:
        print ('\nIteration #', i)
        flux = Area * get_Bc(flux)
        print ('Flux = ', flux)
        i += 1
    return flux


if __name__ == "__main__":
    print ('\nSolving with Newton-Raphson')
    newton_raphson(0)
    print ('\n')

    print ('\nSolving with Successive Substitution')
    successive_sub (0)
    print('\n')
