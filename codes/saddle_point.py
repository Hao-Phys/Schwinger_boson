#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : saddle_point.py
# Author            : Hao Zhang <hzhangphys@gmail.com>
# Date              : 02.07.2020
# Last Modified Date: 02.09.2020
# Last Modified By  : Hao Zhang <hzhangphys@gmail.com>
import numpy as np
import cofig as cf

"""""
This file calculates the saddle-point k-dependent 
values defined in the note. The unit for k is 2*pi/a
"""""

Jex = cf.Jex
A_sign = cf.A_sign
B_sign = cf.B_sign
delta_ij = cf.delta_ij
Q_half = 0.5 * cf.Q_vec 
guess_A = cf.guess_A
guess_B = cf.guess_B
guess_Lambda = cf.guess_Lambda
h_SB = cf.h_SB


def saddle_values(k, sigma, Lambda, A_delta1, B_delta1, rcase):
    """saddle_values 

    :param k: momentum, 
              can be a vector with dim=N
              or a meshgrid with dim=N
    :param sigma: flavor, -1 or 1
    :param Lambda: saddle-point values of lambda
    :param A_delta1: saddle-point values of A_delta1
    :param B_delta1: saddle-point values of B_delta1
    :param rcase: 0, return to saddle-point spinon dispersion
                  1, return to the useful quantites for solving 
                     the saddle-point equations
                  2, return to useful quantities for saddle-point 
                     Green's function
                  others, exit the program
    """

    kQ0 = k[0, ...] + Q_half[0]
    kQ1 = k[1, ...] + Q_half[1]

    mkQ0 = -k[0, ...] + Q_half[0]
    mkQ1 = -k[1, ...] + Q_half[1]

    gammaA_kQ = 0.0
    gammaA_mkQ = 0.0
    gammaB_kQ = 0.0
    gammaB_mkQ = 0.0

    for bond in range(3):

        k_dot_d = kQ0 * delta_ij[bond, 0] + kQ1 * delta_ij[bond, 1]
        mk_dot_d = mkQ0 * delta_ij[bond, 0] + mkQ1 * delta_ij[bond, 1]
        
        gammaA_kQ += Jex * (A_delta1*A_sign[bond] \
                             * np.sin(2.0*np.pi * k_dot_d))

        gammaA_mkQ += Jex * (A_delta1*A_sign[bond] \
                             * np.sin(2.0*np.pi * mk_dot_d))

        gammaB_kQ += Jex * (B_delta1*B_sign[bond] \
                             * np.cos(2.0*np.pi * k_dot_d))

        gammaB_mkQ += Jex * (B_delta1*B_sign[bond] \
                             * np.cos(2.0*np.pi * mk_dot_d))

    alpha_kQ_sq = (Lambda + gammaB_kQ)**2 - gammaA_kQ**2
    alpha_mkQ_sq = (Lambda + gammaB_mkQ)**2 - gammaA_mkQ**2 

    deltak_sq = np.sqrt( (alpha_kQ_sq-alpha_mkQ_sq)**2 + \
            ( (Lambda+gammaB_kQ+Lambda+gammaB_mkQ)**2  \
            - (gammaA_kQ-gammaA_mkQ)**2 ) * h_SB**2 )

    epsilonk = np.sqrt( 0.5*(alpha_kQ_sq + alpha_mkQ_sq + h_SB**2 \
                             + sigma * deltak_sq) )

    # notice that epsilon(k) = epsilon(-k)

    uk_sq = (Lambda+gammaB_kQ)/(2.0*epsilonk) + 0.5
    umk_sq = (Lambda+gammaB_mkQ)/(2.0*epsilonk) + 0.5

    vk_sq = (Lambda+gammaB_kQ)/(2.0*epsilonk) - 0.5
    vmk_sq = (Lambda+gammaB_mkQ)/(2.0*epsilonk) - 0.5

    zk = gammaA_kQ/(2.0*epsilonk)
    zmk = gammaA_mkQ/(2.0*epsilonk)

    delta_ksig_sq = sigma * (epsilonk**2 - alpha_kQ_sq)
    delta_mksig_sq = sigma * (epsilonk**2 - alpha_mkQ_sq)

    Ak = vk_sq * (delta_mksig_sq/deltak_sq) \
       + sigma * umk_sq * (h_SB/(4.0*deltak_sq))     

    Bk = - sigma * epsilonk * ( vk_sq*vmk_sq + zk*zmk \
            - ( h_SB/(4.0*epsilonk) )**2 ) * (h_SB/deltak_sq)

    Ck = zk * delta_mksig_sq / deltak_sq \
            - sigma * zmk * h_SB**2 / (2.0*deltak_sq)

    Dk = - sigma * epsilonk * ( vk_sq*zmk + umk_sq*zk ) * h_SB/(2.0*deltak_sq)

    Ek = uk_sq * (delta_mksig_sq/deltak_sq) + sigma*vmk_sq*h_SB**2/(4.0*deltak_sq)

    Fk = - sigma * epsilonk * ( uk_sq*umk_sq + zk*zmk \
            - ( h_SB/(4.0*epsilonk) )**2 ) * (h_SB/deltak_sq)

    if rcase == 0:
        return epsilonk
    elif rcase == 1:
        return Ak, Ck
    elif rcase == 2:
        return Ak, Bk, Ck, Dk, Ek, Fk
    else:
        print('wrong return type, exiting')
        exit()

