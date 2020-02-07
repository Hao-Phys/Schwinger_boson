#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : cofig.py
# Author            : Hao Zhang <hzhangphys@gmail.com>
# Date              : 02.07.2020
# Last Modified Date: 02.07.2020
# Last Modified By  : Hao Zhang <hzhangphys@gmail.com>
import numpy as np

spin = 0.5
N_s = 50    # number of sites
N = 2       # SP(N=2), or SU(N=2)
Jex = 1.0   # exchange constants
Q_vec = np.array([2/3, 1/3])
num_bond = 3
delta_ij = np.array([[1, 0], [0, 1], [-1, 1]])
# the relative sign of A_deltai/A_delta1 
# and B_deltai/B_delta1
A_sign = np.array([1, 1, -1])
B_sign = np.array([1, -1, -1])

# initial guesses 
guess_Lambda = 2.0 - 3.0/8.0*Jex
guess_B = -0.25 
guess_A = -np.sqrt( 4.0/(27.0*J_ex**2) * (4.0-1.0/N_s) )

# symmetry breaking field
h_SB = 1.0/N_s
