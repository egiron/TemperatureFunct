# coding=utf-8
from __future__ import absolute_import, division, print_function

import warnings
warnings.simplefilter('ignore')
warnings.filterwarnings('ignore')

import numpy as np
import numba
from numba import cuda
from numba import vectorize, int32, int64, float32, float64
from numba import njit, jit, guvectorize


def calcPRFT(TDay, TOpt=18):
    ''' Estimate Photosynthesis reduction factor (PRFT)
        ``` python
        PRFT = 1 - 0.0025 * (TDay - TOpt)^2
        ``` 

        Warning: Deprecated.
            Stop using this function. This function estimate PRFT using only one value for each parameters.

        Parameters:
            TDay (float): Number or array of Day Temperatures
            TOpt (float): Optimum Temperature. Default value 18

        Returns: 
            (float): A number or array of PRFT
    
    '''
    PRFT = 0
    if (TDay > 0):
        PRFT = 1 - 0.0025 * (TDay - TOpt) ** 2
    return PRFT


@guvectorize([(float64[:], float64[:], float64[:])], '(n), ()->(n)')
def PRFT_gu(Tday, Topt, res):
    for i in range(Tday.shape[0]):
        res[i] = 1 - 0.0025*(Tday[i]-Topt[0])**2 if Tday[i] > 0.0 else 0.0


@numba.vectorize([float64(float64, float64)])
def _getPRFT(Tday, Topt): 
    prft = 1 - 0.0025*(Tday-Topt)**2 #if Tday > 0.0 else 0.0
    return prft

@numba.vectorize([float64(float64, float64)])
def getPRFT(Tday, Topt): 
    prft = 1 - 0.0025*(Tday-Topt)**2 if Tday > 0.0 else 0.0
    return prft

@numba.jit(parallel=True)
def apply_PRFT(Tday, Topt=18):
    n = len(Tday)
    result = np.zeros(n, dtype="float64")
    for i in range(n):
        result[i] = getPRFT(Tday[i], Topt)
    return result

def calculatePRFT(Tday, Topt=18):
    ''' Estimate Photosynthesis reduction factor (PRFT) in parallel.

        ``` python
        PRFT = 1 - 0.0025 * (TDay - TOpt)^2
        ``` 

        Parameters:
            Tday (float): Number or array of Day Temperatures
            Topt (float): Optimum Temperature. Default value 18

        Returns: 
            (float): A number or array of PRFT
    
    '''
    if (Tday is None):
        print("Day Temperature parameter is not valid")
        return
    result = []
    try:
        result = apply_PRFT(Tday, Topt)
    except:
        print("Error calculating photosynthesis reduction factor (PRFT)")
    
    return result #pd.Series(result, index=w.index, name="PRFT") int(GDD)