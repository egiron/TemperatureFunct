# coding=utf-8
from __future__ import absolute_import, division, print_function

import warnings
warnings.simplefilter('ignore')
warnings.filterwarnings('ignore')

import numpy as np
import numba
#from numba import cuda
from numba import vectorize, int32, int64, float32, float64
from numba import njit, jit, guvectorize

# --------------------------------------------------------------
# Processing iPAR - Total light interception
# --------------------------------------------------------------
@guvectorize([(float64[:], float64[:])], '(n)->(n)')
def iPAR_gu(ndvi, res):
    for i in range(ndvi.shape[0]):
        ni = ndvi[i] * 1.25 - 0.19
        ni = 0.01 if ni < 0.01 else ni
        res[i] = 0.95 if ni > 0.95 else ni

''' Total light interception - iPAR

    iPAR = NDVI * 1.25 - 0.19 # between heading and maturity (Campos et al. 2018)
    iPAR = NDVI * 1.25 - 0.21 Daughtry et al. (1992)

    - Asrar, G., Fuchs, M., Kanemasu, E.T., Hatfield, J.L., 1984. 
    Estimating absorbed photosynthetic radiation and leaf area index from spectral reflectance 
    in wheat. Agron. J. 76, 30-306.

    - Campos 2018 Remote sensing-based crop biomass with water or light-driven crop growth models in 
        wheat commercial fields.

    Warning: Deprecated.
            Stop using this function.

    Parameters:
        m (object): A tfunct model 
        ndvi (array): Array of float values

    Return: 
        An array of Total light interception values

'''
@numba.vectorize([numba.float64(numba.float64)])
def _getiPAR(ndvi):  # noqa E501
    
    ndvi = ndvi * 1.25 - 0.19
    ndvi = np.where(ndvi < 0.01, 0.01, ndvi)
    ndvi = np.where(ndvi > 0.95, 0.95, ndvi)
    return ndvi #pd.Series(ndvi)


def calcIPAR(ndvi):
    ''' Total light interception - iPAR

        Reference:
            iPAR = NDVI * 1.25 - 0.19 # between heading and maturity (Campos et al. 2018)
            
            iPAR = NDVI * 1.25 - 0.21 Daughtry et al. (1992)

            - Asrar, G., Fuchs, M., Kanemasu, E.T., Hatfield, J.L., 1984. 
            Estimating absorbed photosynthetic radiation and leaf area index from spectral reflectance 
            in wheat. Agron. J. 76, 30-306.

            - Campos 2018 Remote sensing-based crop biomass with water or light-driven crop growth models in 
                wheat commercial fields.

        Warning: Deprecated.
                Stop using this function. This function estimate iPAR using only one value for each parameters.

        Parameters:
            ndvi (array): Array of float values

        Returns: 
            (array): An array of Total light interception values

    '''
    ndvi = ndvi * 1.25 - 0.19
    ndvi = np.where(ndvi < 0.01, 0.01, ndvi)
    ndvi = np.where(ndvi > 0.95, 0.95, ndvi)
    return ndvi


@numba.vectorize([float32(float32), float64(float64)])
def getIPAR(ndvi):
    ndvi = ndvi * 1.25 - 0.19
    ndvi = 0.01 if ndvi < 0.01 else ndvi
    ndvi = 0.95 if ndvi > 0.95 else ndvi
    return float("{:.3f}".format(ndvi))

@numba.jit(parallel=True, nopython=False) 
def apply_IPAR(ndvi):
    n = len(ndvi)
    result = np.empty(n, dtype="float64")
    assert len(ndvi) == n
    for i in range(n):
        result[i] = getIPAR(ndvi[i])
    return result

def estimate_IPAR(ndvi=None):
    ''' Total light interception - iPAR.

        It is assumed that NDVI at maturity is 0.25.

        Reference:
            iPAR = NDVI * 1.25 - 0.19 # between heading and maturity (Campos et al. 2018)
            
            iPAR = NDVI * 1.25 - 0.21 Daughtry et al. (1992)

            - Asrar, G., Fuchs, M., Kanemasu, E.T., Hatfield, J.L., 1984. 
            Estimating absorbed photosynthetic radiation and leaf area index from spectral reflectance 
            in wheat. Agron. J. 76, 30-306.

            - Campos 2018 Remote sensing-based crop biomass with water or light-driven crop growth models in 
                wheat commercial fields.

        Parameters:
            ndvi (array): Array of float values

        Returns: 
            (array): An array of Total light interception values

    '''
    result = []
    if ( (ndvi is None) ):
        print("NDVI data not valid")
        return
    try:
        result = apply_IPAR(ndvi)
    except:
        print("Error calculating iPAR")

    return result 