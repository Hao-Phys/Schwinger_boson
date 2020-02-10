#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : main_sp_eqns.py
# Author            : Hao Zhang <hzhangphys@gmail.com>
# Date              : 02.07.2020
# Last Modified Date: 02.09.2020
# Last Modified By  : Hao Zhang <hzhangphys@gmail.com>

import numpy as np
import cofig as cf
import scipy.optimize as sciopt
import obj_SP as obj

guess_A = cf.guess_A
guess_B = cf.guess_B
AB0 = np.array([guess_A, guess_B])

sciopt.least_squares(obj.opt_AB, AB0, )
