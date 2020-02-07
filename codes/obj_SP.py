#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : obj_SP.py
# Author            : Hao Zhang <hzhangphys@gmail.com>
# Date              : 02.07.2020
# Last Modified Date: 02.07.2020
# Last Modified By  : Hao Zhang <hzhangphys@gmail.com>

import numpy as np
import cofig as cf
"""""
Some trials show that it is difficult to obtain converging
results for all three eqns simultanously in one solve.
Therefore, we will first solve for S, then for A, B.

"""""
def opt_AB(v):

    A_delta = v[0]
    B_delta = v[1]

    # given random A_delta, B_delta, solve Lambda




