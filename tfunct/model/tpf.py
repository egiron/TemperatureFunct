# coding=utf-8
from __future__ import absolute_import, division, print_function

import warnings
warnings.simplefilter('ignore')
warnings.filterwarnings('ignore')

import itertools
import math
import numpy as np
import numba
from numba import njit, jit, prange, guvectorize, vectorize, cuda
from numba import vectorize, int32, int64, float32, float64
from scipy import stats


# ----------------------------------------------
# Trapezoidal Temperature Function (TPF) 
# ----------------------------------------------
def calculate_TPF(Tday, Tmin, Toptmin, Toptmax, Tmax):
    tpf = 0
    if (Toptmin > Toptmax):
        print("Min Optimum Temperature greater than Max Opt. Temperature")
        tpf = np.nan
    elif (Toptmax > Tmax):
        print("Max Optimum Temperature greater than Maximum Temperature")
        tpf = np.nan
    else:
        if ((Tday < Tmin) or (Tday > Tmax)):
            tpf = 0
        elif ((Tday >= Toptmin) and (Tday <= Toptmax)):
            tpf = 1
        elif (Tday < Toptmin):
            gradient, intercept, r_value, p_value, std_err = stats.linregress([Tmin,Toptmin],[0,1])
            tpf = Tday * gradient
        elif (Tday > Toptmax):
            gradient, intercept, r_value, p_value, std_err = stats.linregress([Toptmax, Tmax],[1,0])
            tpf = 1-((Tday-Toptmax)*abs(gradient))
    #
    return tpf

#@numba.jit()#nopython=True
def getSlope(x,y):
    A = np.vstack([x, np.ones(len(x))]).T
    m, c = np.linalg.lstsq(A, y, rcond=None)[0]
    #ret = np.linalg.lstsq(x,y)
    return m #ret[0] #coef

# Trapezoidal Temperature Function (TPF) 
@numba.jit() #nopython=True
def apply_TPF_numba(Tday, Tmin, Toptmin, Toptmax, Tmax):
    tpf = 0.0
    if ((Toptmin > Toptmax) or (Toptmax > Tmax) ):
        tpf = np.nan
    else:
        if ((Tday < Tmin) or (Tday > Tmax)):
            tpf = 0.0
        elif ((Tday >= Toptmin) and (Tday <= Toptmax)):
            tpf = 1.0
        elif (Tday < Toptmin):
            x = np.array([Tmin,Toptmin])
            y = np.array([0.0,1.0])
            slope = getSlope(x,y)
            tpf = Tday * slope
        elif (Tday > Toptmax):
            x = np.array([Toptmax,Tmax])
            y = np.array([1.0,0.0])
            slope = getSlope(x,y)
            tpf = 1-((Tday-Toptmax)*np.abs(slope))
    #
    return tpf


def apply_TPF(col_Tday, Tmin, Toptmin, Toptmax, Tmax):
    n = len(col_Tday)
    result = np.empty(n, dtype="float64")
    assert len(col_Tday) == n
    for i in prange(n):
        result[i] = apply_TPF_numba(col_Tday[i], Tmin, Toptmin, Toptmax, Tmax)
    return result

def compute_TPF(df, Tmin, Toptmin, Toptmax, Tmax):
    df["TPFTMAX"] = apply_TPF( df["Tdaymax"].to_numpy(), Tmin, Toptmin, Toptmax, Tmax )
    return df