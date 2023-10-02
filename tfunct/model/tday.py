# coding=utf-8
from __future__ import absolute_import, division, print_function

import warnings
warnings.simplefilter('ignore')
warnings.filterwarnings('ignore')

import numpy as np
import numba
#from numba import cuda
from numba import int32, int64, float32, float64
from numba import njit, jit, guvectorize, vectorize, cuda
#import math


def calcTDay(Tmin, Tmax, tminFactor=0.25):
    ''' Calculate day time temperature. 
        TDay is a function of weighted Tmin and weighted Tmax.

        ``` python
         TDay = 0.75 * Tmax + 0.25 * TMin
        ``` 

        Warning: Deprecated.
            Stop using this function. This function estimate day time temperature (TDay)
            using only one value for each parameters.

        Parameters:
            Tmin (float): Number or array of Minimum Temperatures
            Tmax (float): Number or array of Maximum Temperatures
            tminFactor (float): Minimum Temperature factor

        Returns: 
            (float): A number of Day Temperatures

    ''' 
    if (Tmax <= Tmin):
        print ("Error: Maximum temperature is equal or lower than minimum temperature")
        return None
    tmaxFactor = 1 - tminFactor
    TDay = tmaxFactor*Tmax + tminFactor*Tmin
    return float("{:.3f}".format(TDay))


#@guvectorize([(float64[:], float64[:], float64[:], float64[:])], '(n), (n), ()->(n)')
#def tDay_gu(tn, tx, tminFactor, res):
#    tmaxFactor = 1 - tminFactor
#    for i in range(tn.shape[0]):
#        res[i] = tmaxFactor[0]*tx[i] + tminFactor[0]*tn[i]

''' Process one location at a time by estimating day time temperature'''
@guvectorize([(float64[:], float64[:], float64[:], float64[:])], '(n), (n), ()->(n)', target="cpu", nopython=True)
def tDay_gu(tn, tx, tminFactor, res):
    tmaxFactor = 1 - tminFactor[0]
    for i in range(tn.shape[0]):
        res[i] = tmaxFactor*tx[i] + tminFactor[0]*tn[i]

''' Process one location at a time by estimating day time temperature'''
@guvectorize([(float64[:], float64[:], float64[:], float64[:])], '(n), (n), ()->(n)', target="parallel", nopython=True)
def tDay_gu_parallel(tn, tx, tminFactor, res):
    tmaxFactor = 1 - tminFactor[0]
    for i in range(tn.shape[0]):
        res[i] = tmaxFactor*tx[i] + tminFactor[0]*tn[i]

''' Process one location at a time by estimating day time temperature'''
@guvectorize(['void(float64[:], float64[:], float64[:], float64[:])'], '(n), (n), ()->(n)', target='cuda' if cuda.is_available() else 'cpu')
def tDay_gu_cuda(tn, tx, tminFactor, res):
    tmaxFactor = 1 - tminFactor[0]
    for i in range(tn.shape[0]):
        res[i] = tmaxFactor*tx[i] + tminFactor[0]*tn[i]

@numba.vectorize([
            float32(float32, float32, float32),
            float64(float64, float64, float64)])
def _getTDay(Tmin, Tmax, tminFactor):
    ''' ''' 
    #if (Tmax <= Tmin):
    #    print ("Error: Maximum temperature is equal or lower than minimum temperature")
    #    return None
    tmaxFactor = 1 - tminFactor
    TDay = tmaxFactor*Tmax + tminFactor*Tmin
    return float("{:.2f}".format(TDay))

##@numba.jit(parallel=True, nopython=False)
@numba.vectorize([float32(float32, float32, float32), float64(float64, float64, float64)])
def getTDay(Tmin, Tmax, tminFactor=0.25):
    ''' ''' 
    if (Tmax <= Tmin):
        print ("Error: Maximum temperature is equal or lower than minimum temperature")
        return None
    tmaxFactor = 1 - tminFactor
    TDay = tmaxFactor*Tmax + tminFactor*Tmin
    return float("{:.3f}".format(TDay))

@numba.jit(parallel=True, nopython=False) 
def apply_TDay(Tmin, Tmax, tminFactor):
    ''' ''' 
    n = len(Tmin)
    result = np.empty(n, dtype="float64")
    assert len(Tmin) == len(Tmax) == n
    for i in range(n):
        result[i] = getTDay(Tmin[i], Tmax[i], tminFactor)
    return result

def estimate_TDay(Tmin=None, Tmax=None, tminFactor=0.25):
    ''' An optimized function to calculate day time temperature in parallel.

        ``` python
         TDay = 0.75 * Tmax + 0.25 * TMin
        ``` 

        Parameters:
            Tmin (float): Number or array of Minimum Temperatures
            Tmax (float): Number or array of Maximum Temperatures
            tminFactor (float): Minimum Temperature factor

        Returns: 
            (float): A number of Day Temperatures

    '''
    result = []
    if ( (Tmin is None) or (Tmax is None) ):
        print("Weather data not valid")
        return
    try:
        result = apply_TDay(Tmin, Tmax, tminFactor )
    except:
        print("Error calculating Day temperature")

    return result #pd.Series(result, index=w.index, name="TDay")



# ------------------------
# Optimized process
# ------------------------
def prepareTDay_dataset(sites):
    arr_tn = []
    arr_tx = []
    tminFactor = [0.25]
    maxsize = 0
    # Obtener el arreglo de mayor longitud en todas las observaciones
    for _id in range(0, len(sites)):
        tn = sites[_id].inputWPN['TMIN'].to_numpy()
        tx = sites[_id].inputWPN['TMAX'].to_numpy()
        if (len(tn)>maxsize):
            maxsize = len(tn)
        if (len(tx)>maxsize):
            maxsize = len(tx)
    # Con la matriz de mayor longitud crear los arreglos iguales de 2D 
    # con los valores de Tmin y Tmax
    for _id in range(0, len(sites)):
        tn = sites[_id].inputWPN['TMIN'].to_numpy()
        tx = sites[_id].inputWPN['TMAX'].to_numpy()
        tn2 = np.ones(maxsize, dtype=np.float64)
        tx2 = np.ones(maxsize, dtype=np.float64)
        for i in range(maxsize):
            try:
                tn2[i] = tn[i]
                tx2[i] = tx[i]
            except:
                tn2[i] = np.nan
                tx2[i] = np.nan

        arr_tn.append(tn2)
        arr_tx.append(tx2)
        #print(sites[_id].attributes['UID'], len(tn), len(tx))

    arr_tn = np.array(arr_tn)
    arr_tx = np.array(arr_tx)
    rows, cols = arr_tn.shape
    output_array = np.ones(rows*cols, dtype=np.float64).reshape(rows, cols)
    #print(output_array.shape)
    
    return arr_tn, arr_tx, output_array

''' Process all observations at the same time by estimating day time temperature'''
@guvectorize(['float64[:,:], float64[:,:], float64[:], float64[:,:]'], '(m,n), (m,n), () -> (m,n)', target='cpu')
def tday_cpu(tn, tx, tminFactor, result):
    m, n = tn.shape  # m = num of obs, n = num of days from heading to maturity
    assert len(tn) == len(tx) == m
    tmaxFactor = 1 - tminFactor[0]
    for obs in range(m):
        tmp_result = np.zeros(n, dtype=np.float64)
        for i in range(n):
            tmp_result[i] = tn[obs][i] * tminFactor[0] + tx[obs][i] * tmaxFactor
        result[obs, :] = tmp_result

''' Process all observations at the same time by estimating day time temperature'''
@guvectorize(['float64[:,:], float64[:,:], float64[:], float64[:,:]'], '(m,n), (m,n), () -> (m,n)', target='parallel')
def tday_parallel(tn, tx, tminFactor, result):
    m, n = tn.shape  # m = num of obs, n = num of days from heading to maturity
    assert len(tn) == len(tx) == m
    tmaxFactor = 1 - tminFactor[0]
    for obs in range(m):
        tmp_result = np.zeros(n, dtype=np.float64)
        for i in range(n):
            tmp_result[i] = tn[obs][i] * tminFactor[0] + tx[obs][i] * tmaxFactor
        result[obs, :] = tmp_result

''' Process all observations at the same time by estimating day time temperature'''
@guvectorize(['float64[:,:], float64[:,:], float64[:], float64[:], float64[:,:]'], '(m,n), (m,n), (), (n) -> (m,n)', target='cuda' if cuda.is_available() else 'cpu')
def tday_cuda(tn, tx, tminFactor, tmp_result, result):
    m, n = tn.shape  # m = num of obs, n = num of days from heading to maturity
    assert len(tn) == len(tx) == m
    tmaxFactor = 1 - tminFactor[0]
    for obs in range(m):
        #tmp_result = np.zeros(n, dtype=np.float64) # numpy function will not work here!!!
        for i in range(n):
            tmp_result[i] = tn[obs][i] * tminFactor[0] + tx[obs][i] * tmaxFactor
            result[obs, i] = tmp_result[i]

# ------------------------
