# coding=utf-8
from __future__ import absolute_import, division, print_function

import warnings
warnings.simplefilter('ignore')
warnings.filterwarnings('ignore')

import math
import numpy as np
import numba
from numba import njit, jit, prange, guvectorize, vectorize, cuda
from numba import vectorize, int32, int64, float32, float64


# -----------------------------------------------
# Yield Estimation
# We assumed that the accumulated dry matter from heading to maturity converts to grains and final 
# yield. By subtracting 20% of dry matter as the grain moisture and respiration, 
# final yield was estimated as follows:
# -----------------------------------------------

@guvectorize([(float64[:,:,:], float64[:,:], float64[:], float64[:], float64[:,:] )], 
             '(m,n,o), (p, q), (), () -> (p,n)', nopython=True)
def estimate( data, params, is_VPDStress, tf, results):
    ''' An optimized function for estimating grain yield in one step for all observations using one of the three 
        Temperature Function (TF) such as PRFT, WETF and TPF.

        Parameters: 
            data (object): Array of arrays containing _tn, tx, ndvi, solrad, VPDx, ipar, GPP_ datasets
            params (object): An array with _RUE, tminFactor, Topt, Lvpd, Uvpd, SFvpd_Lthres, SFvpd_Uthres_ values
            is_VPDStress (array): Array of one value with T/F if Vapor Pressure Deficit stress affecting grain yield. Default is [False].
            tf (array): Array of one value representing the type of function. [1]: PRFT; [2]: WETF; [3]: TPF.
            results (array): Empty array for outputs.

        Returns:
            results (array): An array with estimated yield for each site.

    '''
    tn, tx, ndvi, solrad, VPDx, ipar, GPP = data
    m, n = params.shape
    for r in range(m):
        p = params[r]
        if (is_VPDStress[0]==True):
            if (tf[0]==1):
                RUE, tminFactor, Topt, Lvpd, Uvpd, SFvpd_Lthres,SFvpd_Uthres = p
            elif (tf[0]==2):
                RUE, Tmin, Topt, Tmax, tminFactor, Lvpd, Uvpd, SFvpd_Lthres, SFvpd_Uthres = p
            elif (tf[0]==3):
                RUE, Tmin, Toptmin, Toptmax, Tmax, tminFactor, Lvpd, Uvpd, SFvpd_Lthres, SFvpd_Uthres = p
        else:
            if (tf[0]==1):
                RUE, tminFactor, Topt = p
            elif (tf[0]==2):
                RUE, Tmin, Topt, Tmax, tminFactor = p
            elif (tf[0]==3):
                RUE, Tmin, Toptmin, Toptmax, Tmax, tminFactor = p
        tmaxFactor = 1 - tminFactor
        alpha = math.log(2) / math.log((Tmax - Tmin)/(Topt - Tmin))
        beta = 1
        gpp2 = GPP.copy()
        Tday = GPP.copy()
        temfun = GPP.copy()
        sfvpd = GPP.copy()
        for j in range(tn.shape[0]):
            Tday[j] = tmaxFactor*tx[j] + tminFactor*tn[j]
            # Estimate GPP
            if (tf[0]==1): # PRFT
                for td in range(Tday[j].size):
                    if Tday[j][td] < 0.0:
                        Tday[j][td] = 0.0
                temfun[j] = 1.0 - 0.0025*(Tday[j]-Topt)**2
                
            elif (tf[0]==2): # WETF
                for td in range(Tday[j].size):
                    if Tday[j][td] < 0.0:
                        Tday[j][td] = 0.0
                    if (Topt <= Tmin):
                        temfun[j][td] = 0.0
                    else:
                        temfun[j][td] = ( (2 * (Tday[j][td] - Tmin)**alpha * (Topt - Tmin)**alpha - (Tday[j][td] - Tmin)**(2*alpha)) / (Topt - Tmin)**(2*alpha) )**beta
                        temfun[j][td] = 0 if (Tday[j][td] < Tmin) else temfun[j][td]
                        temfun[j][td] = 0 if (Tday[j][td] > Tmax) else temfun[j][td]
                    
            elif (tf[0]==3): # TPF
                for td in range(Tday[j].size):
                    if Tday[j][td] < 0.0:
                        Tday[j][td] = 0.0
                    temfun[j][td] = 0.0
                    if ((Toptmin > Toptmax) or (Toptmax > Tmax) ):
                        temfun[j][td] = np.nan
                    else:
                        if ((Tday[j][td] < Tmin) or (Tday[j][td] > Tmax)):
                            temfun[j][td] = 0.0
                        elif ((Tday[j][td] >= Toptmin) and (Tday[j][td] <= Toptmax)):
                            temfun[j][td] = 1.0
                        elif (Tday[j][td] < Toptmin):
                            x = [Tmin,Toptmin]
                            y = [0.0,1.0]
                            A = [[ Tmin,  0.],[Toptmin,  1.]]
                            ave_x = (A[0][0] + A[1][0]) / len(A) 
                            ave_y = (A[0][1] + A[1][1]) / len(A)
                            sum_x = (A[0][0] * (A[0][0] - ave_x)) + (A[1][0] * (A[1][0] - ave_x)) 
                            sum_y = (A[0][0] * (A[0][1] - ave_y)) + (A[1][0] * (A[1][1] - ave_y))
                            slope = sum_y / sum_x
                            temfun[j][td] = Tday[j][td] * slope
                        elif (Tday[j][td] > Toptmax):
                            x = [Toptmax,Tmax] 
                            y = [1.0,0.0]
                            A = [[ Toptmax,  1.],[Tmax,  0.]]
                            ave_x = (A[0][0] + A[1][0]) / len(A) 
                            ave_y = (A[0][1] + A[1][1]) / len(A) 
                            sum_x = (A[0][0] * (A[0][0] - ave_x)) + (A[1][0] * (A[1][0] - ave_x)) 
                            sum_y = (A[0][0] * (A[0][1] - ave_y)) + (A[1][0] * (A[1][1] - ave_y))
                            slope = sum_y / sum_x
                            temfun[j][td] = 1-((Tday[j][td]-Toptmax)*np.abs(slope))
            # Stress Factor
            sfvpd[j] = 1.0
            if (is_VPDStress[0]==True):
                for k, vpd in enumerate(VPDx[j]):
                    if( vpd <= 0 ):
                        sfvpd[j][k] = SFvpd_Lthres
                    elif ((vpd > 0) and (vpd <= Lvpd )):
                        sfvpd[j][k] = SFvpd_Uthres
                    elif ((vpd > Lvpd) and (vpd < Uvpd )):
                        sfvpd[j][k] = 1 - (vpd-1.0)/(4.6-1.0)
                    else:
                        sfvpd[j][k] = SFvpd_Lthres
            gpp2[j] = solrad[j]*0.5*RUE*temfun[j]*ipar[j]*sfvpd[j]
            simYield = 0.0
            for val in gpp2[j][~np.isnan(gpp2[j])]:
                simYield += val
            results[r,j] = round(simYield * 0.008, 2) #0.8 * 0.01
#

