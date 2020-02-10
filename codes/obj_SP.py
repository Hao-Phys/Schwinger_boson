#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : obj_SP.py
# Author            : Hao Zhang <hzhangphys@gmail.com>
# Date              : 02.07.2020
# Last Modified Date: 02.09.2020
# Last Modified By  : Hao Zhang <hzhangphys@gmail.com>

import numpy as np
import cofig as cf
import saddle_point as sp
import scipy.optimize as sciopt

"""""
Some trials show that it is difficult to obtain converging
results for all three eqns simultanously in one solve.
Therefore, we will first solve for S, then for A, B.
The self-consistent equation is equivalent to minimize
object functions defined below.

"""""
N1 = cf.N1
N2 = cf.N2
N_s = cf.N_s
spin = cf.spin
Q_half = 0.5*cf.Q_vec
Jex = cf.Jex

k1 = np.linspace(0, 1, N1)
k2 = np.linspace(0, 1, N2)

K1, K2 = np.meshgrid(k1, k2)
kk = np.zeros((2, N1, N2))
kk[0, ...] = K1 
kk[1, ...] = K2
phase_sin = np.sin(2.0*np.pi * (K1 + Q_half[0]))
phase_cos = np.cos(2.0*np.pi * (K1 + Q_half[0]))

def opt_S(v, A_in, B_in):

    Lambda_tmp = v

    Akp, Ckp = sp.saddle_values(kk,  1, Lambda_tmp, A_in, B_in, 1)
    Akm, Ckm = sp.saddle_values(kk, -1, Lambda_tmp, A_in, B_in, 1) 

    value = (Akp + Akm).sum(axis=(0, 1))/N_s

# provide the residue for scipy least square opt algorithm
    obj_S = spin - value   

    return obj_S

guess_A = cf.guess_A
guess_B = cf.guess_B
print("guessA=", guess_A)
print("guessB=", guess_B)
AB0 = np.array([guess_A, guess_B])


def opt_AB(v):

    A_delta = v[0]
    B_delta = v[1]

# In the optimization, we first solve for Lambda, because the object 
# function is steeper for Lambda (flatter for A_delta and B_delta). 
# Then we can fine tune for A_delta and B_delta

# given random A_delta, B_delta, solve Lambda
# epsilon_{0,-} = 9*J_ex/(16*N_s)
    Lambda0 = ( ((9*Jex)/(16*N_s))**2 + 0.25*27*Jex**2*A_delta**2 )**(0.5) \
            - 1.5*Jex*B_delta + 0.5/N_s
    print("Lambda0 = ", Lambda0)
    # Lambda0 = cf.guess_Lambda
    res1 = sciopt.least_squares(opt_S, Lambda0, \
                                bounds=(Lambda0-0.2, Lambda0+0.2), \
                                args=(A_delta, B_delta), verbose=2)

    Lambda_opt = res1.x
    print("the optimized cost function value for Lamda = ", res1.cost)

# use the optimized Lambda_opt as input for evaluate the cost function for
# A and B 
    Akp, Ckp = sp.saddle_values(kk,  1, Lambda_opt, A_delta, B_delta, 1)
    Akm, Ckm = sp.saddle_values(kk, -1, Lambda_opt, A_delta, B_delta, 1) 

    A_sum = ((Ckp + Ckm) * phase_sin).sum(axis=(0, 1))/N_s
    B_sum = ((Ckp + Akm) * phase_cos).sum(axis=(0, 1))/N_s

    resi_A = A_delta - A_sum
    resi_B = B_delta - B_sum
# according to sciopt.optimize documentation, return type must be an numpy 
# array with two objective residues
    resi = np.array([resi_A, resi_B])

    return resi


# res = opt_AB(AB0)
# print((res**2).sum(axis=0))
res = sciopt.least_squares(opt_AB, AB0, \
                     bounds=(AB0-0.05, AB0+0.05), verbose=1)

print("final cost function value for AB", res.cost)

x1 = res.x[0]
x2 = res.x[1]
print(x1)
print(x2)
