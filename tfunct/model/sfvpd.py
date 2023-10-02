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

def calcSFvpd(VPDx, Lvpd, Uvpd, SFvpd_Lthres, SFvpd_Uthres):
    ''' '''
    sfvpd = 0
    if( VPDx <= 0 ):
        sfvpd = SFvpd_Lthres 
    elif ((VPDx > 0) and (VPDx <= Lvpd)):
        sfvpd = SFvpd_Uthres #1   
    elif ((VPDx > Lvpd) and (VPDx < Uvpd)):
        sfvpd = 1 - (VPDx-1)/(4.6-1)
    else:
        sfvpd = SFvpd_Lthres 
    
    return sfvpd

@guvectorize([(float64[:], float64[:], float64[:], float64[:], float64[:], float64[:])], 
             '(n), (), (), (), () ->(n)')
def SFvpd_gu(VPDx, Lvpd, Uvpd, SFvpd_Lthres, SFvpd_Uthres, res):
    for i in range(VPDx.shape[0]):
        sfvpd = 0
        if( VPDx[i] <= 0 ):
            sfvpd = SFvpd_Lthres[0]
        elif ((VPDx[i] > 0) and (VPDx[i] <= Lvpd[0])):
            sfvpd = SFvpd_Uthres[0] #1   
        elif ((VPDx[i] > Lvpd[0]) and (VPDx[i] < Uvpd[0])):
            sfvpd = 1 - (VPDx[i]-1)/(4.6-1)
        else:
            sfvpd = SFvpd_Lthres[0] 
        res[i] = sfvpd

@numba.vectorize([float64(float64, float64, float64, float64, float64)])
def getSFvpd(VPDx, Lvpd, Uvpd, SFvpd_Lthres, SFvpd_Uthres): 
    ''' '''
    sfvpd = 0
    if( VPDx <= 0 ):
        sfvpd = SFvpd_Lthres 
    elif ((VPDx > 0) and (VPDx <= Lvpd)):
        sfvpd = SFvpd_Uthres #1   
    elif ((VPDx > Lvpd) and (VPDx < Uvpd)):
        sfvpd = 1 - (VPDx-1)/(4.6-1)
    else:
        sfvpd = SFvpd_Lthres 
    
    return sfvpd

@numba.jit(parallel=True)
def apply_SFvpd(VPDMAX, Lvpd=1, Uvpd=4, SFvpd_Lthres=0.2, SFvpd_Uthres=1):
    ''' 
        
    '''
    n = len(VPDMAX)
    result = np.zeros(n, dtype="float64")
    for i in range(n):
        result[i] = getSFvpd(VPDMAX[i], Lvpd, Uvpd, SFvpd_Lthres, SFvpd_Uthres)
    return result

def calculateSFvpd(VPDMAX, Lvpd=1, Uvpd=4, SFvpd_Lthres=0.2, SFvpd_Uthres=1):
    ''' Calculation of Vapor pressure deficit (VPD) stress factor

        Parameters:
            VPDMAX (array): Array of daily temperature values
            Lvpd (float): A number for threshold of lower VPD. Default is 1
            Uvpd (array): A number for threshold of upper VPD. Default is 4
            SFvpd_Lthres (array): A number for threshold of stress factor of lower VPD. Default is 0.2
            SFvpd_Uthres (array): A number for threshold of stress factor of upper VPD. Default is 1
            
        Returns: 
            (array): A number or array of stressed factors of VPD
    
    '''
    if (VPDMAX is None):
        print("VPDMAX parameter is not valid")
        return
    result = []
    try:
        result = apply_SFvpd(VPDMAX, Lvpd, Uvpd, SFvpd_Lthres, SFvpd_Uthres)
    except:
        print("Error calculating VPD stress")
    
    return result  