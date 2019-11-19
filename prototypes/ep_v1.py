#/usr/bin/python3.5

import numpy as np
from time import time

N=3
MAXITER=10**3

alpha = 3.0
h=1./N

TOL = 10**-5 * h


U_old = np.zeros((N+1,N+1))
U_new = np.zeros((N+1,N+1))
F = np.zeros((N+1,N+1)) 


def initialize_F_matrix(F, N):
    for j_y in range(N+1):
        for i_x in range(N+1):
            F[i_x, j_y] = round(6. * np.pi**2 * np.cos(np.pi * (i_x*h)) * np.sin(np.pi * (j_y*h)) , 2) 

    return F


def initialize_U_edges(U_old, N):
    # Parte inferior e superior da Malha
    for i_x in range(N+1):
        U_old[i_x, 0] = round(3. * np.cos(np.pi * (h*i_x)) * np.sin(np.pi * (h*0.)), 2)
        U_old[i_x, N] = round(3. * np.cos(np.pi * (h*i_x)) * np.sin(np.pi * (h*N )), 2)

    # Parte esquerda e direita da Malha
    for j_y in range(N+1):
        U_old[0, j_y] = round(3. * np.cos(np.pi * (h*0.)) * np.sin(np.pi * (h*j_y)), 2)
        U_old[N, j_y] = round(3. * np.cos(np.pi * (h*N) ) * np.sin(np.pi * (h*j_y)), 2)

    return U_old
