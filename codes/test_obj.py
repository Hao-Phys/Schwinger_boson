#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : test_obj.py
# Author            : Hao Zhang <hzhangphys@gmail.com>
# Date              : 02.09.2020
# Last Modified Date: 02.09.2020
# Last Modified By  : Hao Zhang <hzhangphys@gmail.com>

# test the behaviors of the objective function
# plot residuals square for lambda, A, and B
# by using the condensate values for the other two

import numpy as np
import cofig as cf
import saddle_point as sp
import matplotlib.pyplot as plt

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

guess_A = cf.guess_A
guess_B = cf.guess_B
guess_Lambda = cf.guess_Lambda

Lambda_cut = np.linspace(0.50, 1.0, 100)
A_cut = np.linspace(0.35, 0.42, 50)
B_cut = np.linspace(0.24, 0.252, 50)

Lambda_res = np.zeros(len(Lambda_cut))
A_res = np.zeros(len(A_cut))
B_res = np.zeros(len(B_cut))

for flag in range(len(Lambda_res)):

    Akp, Ckp = sp.saddle_values(kk,  1, Lambda_cut[flag], guess_A, guess_B, 1)
    Akm, Ckm = sp.saddle_values(kk, -1, Lambda_cut[flag], guess_A, guess_B, 1)
    S_sum = (Akp + Akm).sum(axis=(0, 1))/N_s
    A_sum = ( (Ckp + Ckm) * phase_sin ).sum(axis=(0,1))/N_s
    B_sum = ( (Akp + Akm) * phase_cos ).sum(axis=(0,1))/N_s

    Lambda_res[flag] = (spin-S_sum)**2 + (guess_A-A_sum)**2 + (guess_B-B_sum)**2

for flag in range(len(A_res)):
    Akp, Ckp = sp.saddle_values(kk,  1, guess_Lambda, A_cut[flag], guess_B, 1)
    Akm, Ckm = sp.saddle_values(kk, -1, guess_Lambda, A_cut[flag], guess_B, 1)
    S_sum = (Akp + Akm).sum(axis=(0, 1))/N_s
    A_sum = ( (Ckp + Ckm) * phase_sin ).sum(axis=(0,1))/N_s
    B_sum = ( (Akp + Akm) * phase_cos ).sum(axis=(0,1))/N_s

    A_res[flag] = (spin-S_sum)**2 + (guess_A-A_sum)**2 + (guess_B-B_sum)**2

for flag in range(len(B_res)):
    Akp, Ckp = sp.saddle_values(kk,  1, guess_Lambda, guess_A, B_cut[flag], 1)
    Akm, Ckm = sp.saddle_values(kk, -1, guess_Lambda, guess_A, B_cut[flag], 1)
    S_sum = (Akp + Akm).sum(axis=(0, 1))/N_s
    A_sum = ( (Ckp + Ckm) * phase_sin ).sum(axis=(0,1))/N_s
    B_sum = ( (Akp + Akm) * phase_cos ).sum(axis=(0,1))/N_s

    B_res[flag] = (spin-S_sum)**2 + (guess_A-A_sum)**2 + (guess_B-B_sum)**2

Akp, Ckp = sp.saddle_values(kk,  1, guess_Lambda, guess_A, guess_B, 1)
Akm, Ckm = sp.saddle_values(kk, -1, guess_Lambda, guess_A, guess_B, 1)

S_sum = (Akp + Akm).sum(axis=(0, 1))/N_s
A_sum = ( (Ckp + Ckm) * phase_sin ).sum(axis=(0,1))/N_s
B_sum = ( (Akp + Akm) * phase_cos ).sum(axis=(0,1))/N_s

guess_res = (spin-S_sum)**2 + (guess_A-A_sum)**2 + (guess_B-B_sum)**2

fig, ax = plt.subplots(3, 1)
ax[0].plot(Lambda_cut, Lambda_res, 'r')
ax[0].plot(guess_Lambda, guess_res, 'b*')

ax[1].plot(A_cut, A_res, 'g')
ax[1].plot(guess_A, guess_res, 'r*')

ax[2].plot(B_cut, B_res, 'b')
ax[2].plot(guess_B, guess_res, 'g*')

ax[0].set_ylim([0, 0.5])
ax[1].set_ylim([0, 0.5])
ax[2].set_ylim([0, 0.5])

plt.show()


