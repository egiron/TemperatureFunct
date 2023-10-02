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

def calcGPP(SolRad, TF, iPAR, RUE=3.0, stressFactor=1.0):
    ''' Estimate the Gross primary production

        Parameters:
            SolRad (float): Solar Radiation
            TF (float): Temperature function result (PRFT, WETF, TPF)
            iPAR (float): the photosynthetically active radiation (PAR) intercepted by a plant or crop estimated from NDVI
            RUE (float): Radiation-use efficiency. Default value is 3.0 gMJ^−1
            stressFactor (float): Stress Factor (eg. VPD stress factor)

        Warning: Deprecated.
            This function only use one value for each paramters. it is not optimized for run models in parallel.

        Returns:
            (float): the gross primary production (GPP)

            
    '''
    return SolRad*0.5*RUE*TF*iPAR*stressFactor


@guvectorize([(float64[:], float64[:], float64[:], float64[:], float64[:])], 
             '(n), (n), (n), ()->(n)')
def GPP_gu(SolRad, PRFT, iPAR, RUE, res):
    for i in range(SolRad.shape[0]):
        res[i] = SolRad[i]*0.5*RUE[0]*PRFT[i]*iPAR[i]

@guvectorize([(float64[:], float64[:], float64[:], float64[:], float64[:], float64[:])], 
             '(n), (n), (n), (n), () -> (n)')
def GPP_VPDStress_gu(SolRad, PRFT, iPAR, SFvpd, RUE, res):
    for i in range(SolRad.shape[0]):
        res[i] = SolRad[i]*0.5*RUE[0]*PRFT[i]*iPAR[i]*SFvpd[i]


@numba.vectorize([float64(float64, float64, float64, float64, float64)])
def _getGPP(SolRad, PRFT, iPAR, RUE, stressFactor): 
    return SolRad*0.5*RUE*PRFT*iPAR*stressFactor

@numba.vectorize([float64(float64, float64, float64, float64, float64)])
def getGPP(SolRad, PRFT, iPAR, RUE, stressFactor): 
    gpp = SolRad*0.5*RUE*PRFT*iPAR*stressFactor
    return gpp

@numba.jit(parallel=True)
def apply_GPP_v0(SolRad, PRFT, iPAR, RUE=3.0, stressFactor=1.0):
    ''' '''
    n = len(SolRad)
    result = np.zeros(n, dtype="float64")
    assert len(iPAR) == len(PRFT) == n
    for i in range(n):
        result[i] = getGPP(SolRad[i], PRFT[i], iPAR[i], RUE, stressFactor)
    return result

def calculateGPP(SolRad, PRFT, iPAR, RUE=3.0, stressFactor=1.0):
    ''' 
        Estimate the Gross primary production. 
        The accumulated dry matter from heading to maturity

        ``` python
        # GPP = Solar Radiation × TemFun × VPD_StressFactor × iPAR_fromNDVI × RUE 
        GPP = SolRad * 0.5 * RUE * PRFT * iPAR * stressFactor
        ```

        Parameters:
            SolRad (float): Solar Radiation
            PRFT (float): Photosynthesis reduction factor
            iPAR (float): the photosynthetically active radiation (PAR) intercepted by a plant or crop
            RUE (float): Radiation-use efficiency. Default value is 3.0 gMJ^-1
            stressFactor (float): Stress Factor

        Returns:
            (float): a number or array of values with the gross primary production (GPP)

    '''
    if (SolRad is None):
        print("Solar Radiation values are not valid")
        return
    if (PRFT is None):
        print("PRFT values are not valid")
        return
    if (iPAR is None):
        print("iPAR values are not valid")
        return
    if (RUE is None):
        print("RUE is not valid")
        return
    
    result = []
    try:
        result = apply_GPP_v0(SolRad, PRFT, iPAR, RUE, stressFactor)
    except:
        print("Error calculating GPP")
    
    return result


# ---------------------------------
# VPD stress
# ---------------------------------
@numba.jit(parallel=True)
def applyGPP_VPDStress(SolRad, PRFT, iPAR, SFvpd, RUE=3.0 ):
    ''' 
        
    '''
    n = len(SolRad)
    result = np.zeros(n, dtype="float64")
    assert len(iPAR) == len(PRFT) == len(SFvpd) == n
    for i in range(n):
        result[i] = getGPP(SolRad[i], PRFT[i], iPAR[i], RUE, SFvpd[i] )
    return result

def calculateGPP_VPDStress(SolRad, PRFT, iPAR, SFvpd, RUE=3.0 ):
    ''' Estimate the Gross primary production with stressed VPD

        ``` python
        GPP = SolRad * 0.5 * RUE * PRFT * iPAR * SFvpd
        ```

        Parameters:
            SolRad (float): Solar Radiation
            PRFT (float): Photosynthesis reduction factor
            iPAR (float): the photosynthetically active radiation (PAR) intercepted by a plant or crop
            RUE (float): Radiation-use efficiency. Default value is 3.0
            SFvpd (float): Vapor Pressure Deficit (VPD) stress factor

        Returns:
            (float): a number or array of values with the gross primary production (GPP)
    
    '''
    if (SolRad is None):
        print("Solar Radiation values are not valid")
        return
    if (PRFT is None):
        print("PRFT values are not valid")
        return
    if (iPAR is None):
        print("iPAR values are not valid")
        return
    if (SFvpd is None):
        print("SFvpd values are not valid")
        return
    if (RUE is None):
        print("RUE is not valid")
        return
    
    result = []
    try:
        result = applyGPP_VPDStress(SolRad, PRFT, iPAR, SFvpd, RUE)
    except:
        print("Error calculating GPP stressed VPD")
    
    return result


# ----------------------------------------
# Estimating GPP all observations
# ----------------------------------------
'''
    Convert Dataframe columns in vectors or arrays of 1D and fill up empty values to speed up calculations
    in parallel or cuda platforms
'''
def prepareVectors_dataset(sites):
    
    arr_tn = []
    arr_tx = []
    arr_ndvi = []
    arr_solrad = []
    arr_vpdmax = []
    arr_ipar = []
    maxsize = 0
    # Obtener el arreglo de mayor longitud en todas las observaciones
    for _id in range(0, len(sites)):
        tn = sites[_id].inputWPN['TMIN'].to_numpy()
        tx = sites[_id].inputWPN['TMAX'].to_numpy()
        ndvi = sites[_id].inputWPN['NDVI'].to_numpy()
        solrad = sites[_id].inputWPN['SolRad'].to_numpy()
        vpdmax = sites[_id].inputWPN['VPDMAX'].to_numpy()
        # Una sola evaluación sería suficiente pues todos los vectores deben de tener la misma longitud
        assert len(tn) == len(tx) == len(ndvi) == len(solrad) == len(vpdmax)
        if (len(tn)>maxsize):
            maxsize = len(tn)
        #if (len(tx)>maxsize):
        #    maxsize = len(tx)
        #if (len(ndvi)>maxsize):
        #    maxsize = len(ndvi)
    # Con la matriz de mayor longitud crear los arreglos iguales de 2D 
    # con los valores de Tmin y Tmax
    for _id in range(0, len(sites)):
        tn = sites[_id].inputWPN['TMIN'].to_numpy()
        tx = sites[_id].inputWPN['TMAX'].to_numpy()
        ndvi = sites[_id].inputWPN['NDVI'].to_numpy()
        solrad = sites[_id].inputWPN['SolRad'].to_numpy()
        vpdmax = sites[_id].inputWPN['VPDMAX'].to_numpy()
        
        # iPAR
        ipar = np.zeros(len(ndvi), dtype=np.float64)
        for i in range(0,len(ndvi)):
            ipar[i] = ndvi[i] * 1.25 - 0.19
            ipar[i] = 0.01 if ipar[i] < 0.01 else ipar[i] #ipar = np.where((ipar < 0.01), 0.01, ipar)
            ipar[i] = 0.95 if ipar[i] > 0.95 else ipar[i] #np.where((ipar > 0.95), 0.95, ipar) #
        
        tn2 = np.ones(maxsize, dtype=np.float64)
        tx2 = np.ones(maxsize, dtype=np.float64)
        ndvi2 = np.ones(maxsize, dtype=np.float64)
        solrad2 = np.ones(maxsize, dtype=np.float64)
        vpdmax2 = np.ones(maxsize, dtype=np.float64)
        ipar2 = np.ones(maxsize, dtype=np.float64)
        for i in range(maxsize):
            try:
                tn2[i] = tn[i]
                tx2[i] = tx[i]
                ndvi2[i] = ndvi[i]
                solrad2[i] = solrad[i]
                vpdmax2[i] = vpdmax[i]
                ipar2[i] = ipar[i]
            except:
                tn2[i] = np.nan
                tx2[i] = np.nan
                ndvi2[i] = np.nan
                solrad2[i] = np.nan
                vpdmax2[i] = np.nan
                ipar2[i] = np.nan

        arr_tn.append(tn2)
        arr_tx.append(tx2)
        arr_ndvi.append(ndvi2)
        arr_solrad.append(solrad2)
        arr_vpdmax.append(vpdmax2)
        arr_ipar.append(ipar2)

    arr_tn = np.array(arr_tn)
    arr_tx = np.array(arr_tx)
    arr_ndvi = np.array(arr_ndvi)
    arr_solrad = np.array(arr_solrad)
    arr_vpdmax = np.array(arr_vpdmax)
    arr_ipar = np.array(arr_ipar)
    # Create an array for outputs
    rows, cols = arr_tn.shape
    output_array = np.ones(rows*cols, dtype=np.float64).reshape(rows, cols)
    
    return arr_tn, arr_tx, arr_ndvi, arr_solrad, arr_vpdmax, arr_ipar, output_array


@guvectorize([(float64[:,:,:], float64[:,:], float64[:], float64[:], float64[:,:] )], 
             '(m,n,o), (p, q), (), () -> (p,n)', nopython=True)
def estimate( data, params, is_VPDStress, tf, results):
    ''' An optimized function for estimating GPP for all observations using one of the three 
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
            for ii in range(gpp2[j].size):
                if (gpp2[j][ii] < 0.0) or np.isnan(gpp2[j][ii]):
                    gpp2[j][ii] = 0.0
            results[r,j] = gpp2[r,j]
# ------


''' Estimating GPP in one step for all observations '''
@guvectorize([(float64[:], float64[:], float64[:], float64[:], float64[:], float64[:],
               float64[:], float64[:], float64[:], float64[:], float64[:], float64[:], float64[:], float64[:],
               float64[:])], 
             '(n), (n), (n), (n), (n), (n), (), (), (), (), (), (), (), () -> (n)')
def _GPP(tn, tx, ndvi, solrad, VPDx, ipar, tminFactor, Topt, RUE, is_VPDStress, Lvpd, Uvpd, SFvpd_Lthres, SFvpd_Uthres, GPP):
    tmaxFactor = 1 - tminFactor[0]
    for i in range(tn.shape[0]):
        # iPAR
        #ipar = ndvi[i] * 1.25 - 0.19
        #ipar = 0.01 if ipar < 0.01 else ipar #ipar = np.where((ipar < 0.01), 0.01, ipar)
        #ipar = 0.95 if ipar > 0.95 else ipar #np.where((ipar > 0.95), 0.95, ipar) #
        Tday = tmaxFactor*tx[i] + tminFactor[0]*tn[i]
        prft = 1 - 0.0025*(Tday-Topt[0])**2 if Tday > 0.0 else 0.0
        sfvpd = 1.0
        if (is_VPDStress[0]==True):
            if( VPDx[i] <= 0 ):
                sfvpd = SFvpd_Lthres[0]
            elif ((VPDx[i] > 0) and (VPDx[i] <= Lvpd[0])):
                sfvpd = SFvpd_Uthres[0] #1   
            elif ((VPDx[i] > Lvpd[0]) and (VPDx[i] < Uvpd[0])):
                sfvpd = 1 - (VPDx[i]-1)/(4.6-1)
            else:
                sfvpd = SFvpd_Lthres[0]
        #GPP[i] = solrad[i]*0.5*RUE[0]*prft*ipar[i] if (is_VPDStress[0]==False) else solrad[i]*0.5*RUE[0]*prft*ipar[i]*sfvpd
        GPP[i] = solrad[i]*0.5*RUE[0]*prft*ipar[i]*sfvpd
#
''' Estimating GPP in one step for all observations '''
@guvectorize([(float64[:], float64[:], float64[:], float64[:], float64[:], 
               float64[:], float64[:], float64[:], float64[:], float64[:], float64[:], float64[:], float64[:],
               float64[:])], 
             '(n), (n), (n), (n), (n), (), (), (), (), (), (), (), () -> (n)', nopython=True)
def GPP(tn, tx, solrad, VPDx, ipar, tminFactor, Topt, RUE, is_VPDStress, Lvpd, Uvpd, SFvpd_Lthres, SFvpd_Uthres, GPP):
    
    tmaxFactor = 1 - tminFactor[0]
    for i in range(tn.shape[0]):
        Tday = tmaxFactor*tx[i] + tminFactor[0]*tn[i]
        prft = 1 - 0.0025*(Tday-Topt[0])**2 if Tday > 0.0 else 0.0
        sfvpd = 1.0
        if (is_VPDStress[0]==True):
            if( VPDx[i] <= 0 ):
                sfvpd = SFvpd_Lthres[0]
            elif ((VPDx[i] > 0) and (VPDx[i] <= Lvpd[0])):
                sfvpd = SFvpd_Uthres[0] #1   
            elif ((VPDx[i] > Lvpd[0]) and (VPDx[i] < Uvpd[0])):
                sfvpd = 1 - (VPDx[i]-1)/(4.6-1)
            else:
                sfvpd = SFvpd_Lthres[0]
        GPP[i] = solrad[i]*0.5*RUE[0]*prft*ipar[i]*sfvpd
#
''' Estimating GPP in one step for all observations '''
@guvectorize([(float64[:], float64[:], float64[:], float64[:], float64[:], 
               float64[:], float64[:], float64[:], float64[:], float64[:], float64[:], float64[:], float64[:],
               float64[:])], 
             '(n), (n), (n), (n), (n), (), (), (), (), (), (), (), () -> (n)', target='cpu', nopython=True)
def GPP_cpu(tn, tx, ndvi, solrad, VPDx, tminFactor, Topt, RUE, is_VPDStress, Lvpd, Uvpd, SFvpd_Lthres, SFvpd_Uthres, GPP):
    
    tmaxFactor = 1 - tminFactor[0]
    for i in range(tn.shape[0]):
        # NDVI
        ipar = ndvi[i] * 1.25 - 0.19
        ipar = 0.01 if ipar < 0.01 else ipar #ipar = np.where((ipar < 0.01), 0.01, ipar)
        ipar = 0.95 if ipar > 0.95 else ipar #np.where((ipar > 0.95), 0.95, ipar) #
        Tday = tmaxFactor*tx[i] + tminFactor[0]*tn[i]
        prft = 1 - 0.0025*(Tday-Topt[0])**2 if Tday > 0.0 else 0.0
        sfvpd = 0
        if( VPDx[i] <= 0 ):
            sfvpd = SFvpd_Lthres[0]
        elif ((VPDx[i] > 0) and (VPDx[i] <= Lvpd[0])):
            sfvpd = SFvpd_Uthres[0] #1   
        elif ((VPDx[i] > Lvpd[0]) and (VPDx[i] < Uvpd[0])):
            sfvpd = 1 - (VPDx[i]-1)/(4.6-1)
        else:
            sfvpd = SFvpd_Lthres[0]
        GPP[i] = solrad[i]*0.5*RUE[0]*prft*ipar if (is_VPDStress[0]==False) else solrad[i]*0.5*RUE[0]*prft*ipar*sfvpd

#
''' Estimating GPP in one step for all observations '''
@guvectorize([(float64[:], float64[:], float64[:], float64[:], float64[:], 
               float64[:], float64[:], float64[:], float64[:], float64[:], float64[:], float64[:], float64[:],
               float64[:])], 
             '(n), (n), (n), (n), (n), (), (), (), (), (), (), (), () -> (n)', target='parallel', nopython=True)
def GPP_parallel(tn, tx, ndvi, solrad, VPDx, tminFactor, Topt, RUE, is_VPDStress, Lvpd, Uvpd, SFvpd_Lthres, SFvpd_Uthres, GPP):
    
    tmaxFactor = 1 - tminFactor[0]
    for i in range(tn.shape[0]):
        # NDVI
        ipar = ndvi[i] * 1.25 - 0.19
        ipar = 0.01 if ipar < 0.01 else ipar #ipar = np.where((ipar < 0.01), 0.01, ipar)
        ipar = 0.95 if ipar > 0.95 else ipar #np.where((ipar > 0.95), 0.95, ipar) #
        Tday = tmaxFactor*tx[i] + tminFactor[0]*tn[i]
        prft = 1 - 0.0025*(Tday-Topt[0])**2 if Tday > 0.0 else 0.0
        sfvpd = 0
        if( VPDx[i] <= 0 ):
            sfvpd = SFvpd_Lthres[0]
        elif ((VPDx[i] > 0) and (VPDx[i] <= Lvpd[0])):
            sfvpd = SFvpd_Uthres[0] #1   
        elif ((VPDx[i] > Lvpd[0]) and (VPDx[i] < Uvpd[0])):
            sfvpd = 1 - (VPDx[i]-1)/(4.6-1)
        else:
            sfvpd = SFvpd_Lthres[0]
        GPP[i] = solrad[i]*0.5*RUE[0]*prft*ipar if (is_VPDStress[0]==False) else solrad[i]*0.5*RUE[0]*prft*ipar*sfvpd

#
''' Estimating GPP in one step for all observations '''
@guvectorize([(float64[:], float64[:], float64[:], float64[:], float64[:], 
               float64[:], float64[:], float64[:], float64[:], float64[:], float64[:], float64[:], float64[:],
               float64[:])], 
             '(n), (n), (n), (n), (n), (), (), (), (), (), (), (), () -> (n)', target='cuda' if cuda.is_available() else 'cpu', nopython=True)
def GPP_cuda(tn, tx, ndvi, solrad, VPDx, tminFactor, Topt, RUE, is_VPDStress, Lvpd, Uvpd, SFvpd_Lthres, SFvpd_Uthres, GPP):
    
    tmaxFactor = 1 - tminFactor[0]
    for i in range(tn.shape[0]):
        # NDVI
        ipar = ndvi[i] * 1.25 - 0.19
        ipar = 0.01 if ipar < 0.01 else ipar #ipar = np.where((ipar < 0.01), 0.01, ipar)
        ipar = 0.95 if ipar > 0.95 else ipar #np.where((ipar > 0.95), 0.95, ipar) #
        Tday = tmaxFactor*tx[i] + tminFactor[0]*tn[i]
        prft = 1 - 0.0025*(Tday-Topt[0])**2 if Tday > 0.0 else 0.0
        sfvpd = 0
        if( VPDx[i] <= 0 ):
            sfvpd = SFvpd_Lthres[0]
        elif ((VPDx[i] > 0) and (VPDx[i] <= Lvpd[0])):
            sfvpd = SFvpd_Uthres[0] #1   
        elif ((VPDx[i] > Lvpd[0]) and (VPDx[i] < Uvpd[0])):
            sfvpd = 1 - (VPDx[i]-1)/(4.6-1)
        else:
            sfvpd = SFvpd_Lthres[0]
        GPP[i] = solrad[i]*0.5*RUE[0]*prft*ipar if (is_VPDStress[0]==False) else solrad[i]*0.5*RUE[0]*prft*ipar*sfvpd

#


# ------------------------------
# Combinations for grain yield 
# ------------------------------
''' Estimating GPP in one step for all observations '''
@guvectorize([(float64[:], float64[:], float64[:], float64[:], float64[:],  
              float64[:], float64[:], float64[:], float64[:], float64[:], float64[:], float64[:], float64[:],  
               float64[:] )], 
             '(n), (n), (n), (n), (n), (), (), (), (), (), (), (), () -> (n)', nopython=True)
def GYield_v0(tn, tx, solrad, VPDx, ipar,  tminFactor, Topt, RUE, is_VPDStress, 
           Lvpd, Uvpd, SFvpd_Lthres, SFvpd_Uthres, GPP):
    
    tmaxFactor = 1 - tminFactor[0]
    for i in range(tn.shape[0]):
        Tday = tmaxFactor*tx[i] + tminFactor[0]*tn[i]
        prft = 1 - 0.0025*(Tday-Topt[0])**2 if Tday > 0.0 else 0.0
        sfvpd = 1.0
        if (is_VPDStress[0]==True):
            if( VPDx[i] <= 0 ):
                sfvpd = SFvpd_Lthres[0]
            elif ((VPDx[i] > 0) and (VPDx[i] <= Lvpd[0])):
                sfvpd = SFvpd_Uthres[0] #1   
            elif ((VPDx[i] > Lvpd[0]) and (VPDx[i] < Uvpd[0])):
                sfvpd = 1 - (VPDx[i]-1)/(4.6-1)
            else:
                sfvpd = SFvpd_Lthres[0]
        GPP[i] = solrad[i]*0.5*RUE[0]*prft*ipar[i]*sfvpd
        
# -----------------------------------------------
# PRTF - Photosynthesis reduction factor
# -----------------------------------------------
''' An optimized function for estimating grain yield in one step for all observations using the 
    Photosynthesis reduction factor (PRFT) function.

    Parameters: 
        data (object): Array of arrays containing _tn, tx, ndvi, solrad, VPDx, ipar, GPP_ datasets
        params (object): An array with _RUE, tminFactor, Topt, Lvpd, Uvpd, SFvpd_Lthres, SFvpd_Uthres_ values
        is_VPDStress (bool): Vapor Pressure Deficit stress affecting grain yield. Default is False.
        results (array): Empty array for outputs.

    Returns:
        results (array): An array with estimated yield for each site.
'''
@guvectorize([(float64[:,:,:], float64[:,:], float64[:], float64[:,:] )], 
             '(m,n,o), (p, q), () -> (p,n)', nopython=True)
def GYield( data, params, is_VPDStress, results):
    
    tn, tx, ndvi, solrad, VPDx, ipar, GPP = data
    m, n = params.shape
    for r in range(m):
        p = params[r]
        if (is_VPDStress[0]==True):
            RUE, tminFactor, Topt, Lvpd, Uvpd, SFvpd_Lthres,SFvpd_Uthres = p
        else:
            RUE, tminFactor, Topt = p
        tmaxFactor = 1 - tminFactor
        gpp2 = GPP.copy()
        Tday = GPP.copy()
        prft = GPP.copy()
        sfvpd = GPP.copy()
        for j in range(tn.shape[0]):
            Tday[j] = tmaxFactor*tx[j] + tminFactor*tn[j] 
            for td in range(Tday[j].size):
                if Tday[j][td] < 0.0:
                    Tday[j][td] = 0.0
            prft[j] = 1.0 - 0.0025*(Tday[j]-Topt)**2
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
            gpp2[j] = solrad[j]*0.5*RUE*prft[j]*ipar[j]*sfvpd[j]
            simYield = 0.0
            for val in gpp2[j][~np.isnan(gpp2[j])]:
                simYield += val
            results[r,j] = round(simYield * 0.008, 2) #0.8 * 0.01


@numba.jit(parallel=True)
def apply_GPP(data=None, sites=None, sim_params_to_run=None, is_VPDStress=False):
    
    arr_tn, arr_tx, arr_ndvi, arr_solrad, arr_vpdmax, arr_ipar, GPP = data
    simYields = []
    for p in tqdm(sim_params_to_run):
        model.parameters = dict(list(model.parameters.items()) + list(p.items())) #{**model.parameters + **p]
        GYield_v0(arr_tn, arr_tx, arr_solrad, arr_vpdmax, arr_ipar, 
                    [model.parameters['TMIN_PERC_FACTOR']], [model.parameters['Topt']], [model.parameters['RUE']], 
                    [is_VPDStress], [model.parameters['Lvpd']], [model.parameters['Uvpd']], 
                    [model.parameters['SFvpd_Lthres']], [model.parameters['SFvpd_Uthres']], GPP)
        
        simYield = 0.0
        for val in GPP:
            simYield += np.round(val * model.parameters['YIELD_FACTOR'], 2)
            #simYield = simYield[~np.isnan(simYield)]

        simYields.append(simYield)
    return simYields


@numba.jit(parallel=True)
def createDF(df_gy=None, array_params_to_run=None, array_results=None):
    sim_res = []
    for i in range(0, len(array_params_to_run)):
        for j in range(len(df_gy)):
            sim_res.append(np.concatenate([df_gy[j], array_params_to_run[i], [array_results[i][j]]]))
    
    return sim_res


# -----------------------------------------------
# The Wang-Engel temperature function (WETF)
# -----------------------------------------------
#@numba.jit(nopython=True)
def calculate_wang_engel_temperature( Tday, Tmin, Topt, Tmax):
    if (Topt <= Tmin):
        y = 0
    else:
        alpha = math.log(2) / math.log((Tmax - Tmin)/(Topt - Tmin))
        beta = 1 
        y = ( (2 * (Tday - Tmin)**alpha * (Topt - Tmin)**alpha - (Tday - Tmin)**(2*alpha)) / (Topt - Tmin)**(2*alpha) )**beta
        #y = 0 if (Tday < Tmin) else y
        #y = 0 if (Tday > Tmax) else y
        y = np.where(((Tday < Tmin) | (Tday > Tmax) ), 0, y)
    return y

#@numba.jit(parallel=True) 
def apply_WETF(col_Tday, Tmin, Topt, Tmax):
    n = len(col_Tday)
    result = np.empty(n, dtype="float64")
    assert len(col_Tday) == n
    for i in prange(n):
        result[i] = calculate_wang_engel_temperature(col_Tday[i], Tmin, Topt, Tmax)
    return result

def compute_WETF(df, Tmin, Topt, Tmax):
    df["WETFTMAX"] = apply_WETF( df["Tdaymax"].to_numpy(), Tmin, Topt, Tmax )
    return df


''' An optimized function for estimating grain yield in one step for all observations using 
    The Wang-Engel temperature function (WETF) function.

    Parameters: 
        data (object): Array of arrays containing _tn, tx, ndvi, solrad, VPDx, ipar, GPP_ datasets
        params (object): An array with _RUE, tminFactor, Topt, Lvpd, Uvpd, SFvpd_Lthres, SFvpd_Uthres_ values
        is_VPDStress (bool): Vapor Pressure Deficit stress affecting grain yield. Default is False.
        results (array): Empty array for outputs.

    Returns:
        results (array): An array with estimated yield for each site.

'''
@guvectorize([(float64[:,:,:], float64[:,:], float64[:], float64[:,:] )], 
             '(m,n,o), (p, q), () -> (p,n)', nopython=True)
def GYield_WETF( data, params, is_VPDStress, results):
    
    tn, tx, ndvi, solrad, VPDx, ipar, GPP = data
    m, n = params.shape
    for r in range(m):
        p = params[r]
        if (is_VPDStress[0]==True):
            RUE, Tmin, Topt, Tmax, tminFactor, Lvpd, Uvpd, SFvpd_Lthres, SFvpd_Uthres = p
        else:
            RUE, Tmin, Topt, Tmax, tminFactor = p
        alpha = math.log(2) / math.log((Tmax - Tmin)/(Topt - Tmin))
        beta = 1
        tmaxFactor = 1 - tminFactor
        gpp2 = GPP.copy()
        Tday = GPP.copy()
        wetf = GPP.copy()
        sfvpd = GPP.copy()
        for j in range(tn.shape[0]):
            Tday[j] = tmaxFactor*tx[j] + tminFactor*tn[j] 
            for td in range(Tday[j].size):
                if Tday[j][td] < 0.0:
                    Tday[j][td] = 0.0
                # Calculate WETF
                if (Topt <= Tmin):
                    wetf[j][td] = 0.0
                else:
                    wetf[j][td] = ( (2 * (Tday[j][td] - Tmin)**alpha * (Topt - Tmin)**alpha - (Tday[j][td] - Tmin)**(2*alpha)) / (Topt - Tmin)**(2*alpha) )**beta
                    wetf[j][td] = 0 if (Tday[j][td] < Tmin) else wetf[j][td]
                    wetf[j][td] = 0 if (Tday[j][td] > Tmax) else wetf[j][td]
                
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
            gpp2[j] = solrad[j]*0.5*RUE*wetf[j]*ipar[j]*sfvpd[j]
            simYield = 0.0
            for val in gpp2[j][~np.isnan(gpp2[j])]:
                simYield += val
            results[r,j] = round(simYield * 0.008, 2) #0.8 * 0.01
            
#

# -----------------------------------------------
# TPF
# -----------------------------------------------
''' Estimating grain yield in one step for all observations 
    using the Trapezoidal Temperature Function (TPF)

    Warning: Deprecated.
        Stop using this function.
'''
@guvectorize([(float64[:,:,:], float64[:,:], float64[:], float64[:,:] )], 
             '(m,n,o), (p, q), () -> (p,n)', forceobj=True) #, nopython=True)
def GYield_TPF_v0( data, params, is_VPDStress, results):
    
    tn, tx, ndvi, solrad, VPDx, ipar, GPP = data
    m, n = params.shape
    for r in range(m):
        p = params[r]
        if (is_VPDStress[0]==True):
            RUE, Tmin, Toptmin, Toptmax, Tmax, tminFactor, Lvpd, Uvpd, SFvpd_Lthres, SFvpd_Uthres = p
        else:
            RUE, Tmin, Toptmin, Toptmax, Tmax, tminFactor = p
        
        tmaxFactor = 1 - tminFactor
        gpp2 = GPP.copy()
        Tday = GPP.copy()
        tpf = GPP.copy()
        sfvpd = GPP.copy()
        for j in range(tn.shape[0]):
            Tday[j] = tmaxFactor*tx[j] + tminFactor*tn[j] 
            for td in range(Tday[j].size):
                if Tday[j][td] < 0.0:
                    Tday[j][td] = 0.0
                # Calculate TPF
                tpf[j][td] = 0.0
                if ((Toptmin > Toptmax) or (Toptmax > Tmax) ):
                    tpf[j][td] = np.nan
                else:
                    if ((Tday[j][td] < Tmin) or (Tday[j][td] > Tmax)):
                        tpf[j][td] = 0.0
                    elif ((Tday[j][td] >= Toptmin) and (Tday[j][td] <= Toptmax)):
                        tpf[j][td] = 1.0
                    elif (Tday[j][td] < Toptmin):
                        x = np.array([Tmin,Toptmin]) 
                        y = np.array([0.0,1.0]) 
                        #slope = getSlope(x,y)
                        A = np.vstack([x, np.ones(len(x))]).T
                        slope, c = np.linalg.lstsq(A, y, rcond=None)[0]
                        tpf[j][td] = Tday[j][td] * slope
                    elif (Tday[j][td] > Toptmax):
                        x = np.array([Toptmax,Tmax])  
                        y = np.array([1.0,0.0]) 
                        #slope = getSlope(x,y)
                        A = np.vstack([x, np.ones(len(x))]).T
                        slope, c = np.linalg.lstsq(A, y, rcond=None)[0]
                        tpf[j][td] = 1-((Tday[j][td]-Toptmax)*np.abs(slope))
                    
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
            gpp2[j] = solrad[j]*0.5*RUE*tpf[j]*ipar[j]*sfvpd[j]
            simYield = 0.0
            for val in gpp2[j][~np.isnan(gpp2[j])]:
                simYield += val
            results[r,j] = round(simYield * 0.008, 2) 


''' An optimized function for estimating grain yield in one step for all observations using 
    the Trapezoidal Temperature Function (TPF).

    Parameters: 
        data (object): Array of arrays containing _tn, tx, ndvi, solrad, VPDx, ipar, GPP_ datasets
        params (object): An array with _RUE, tminFactor, Topt, Lvpd, Uvpd, SFvpd_Lthres, SFvpd_Uthres_ values
        is_VPDStress (bool): Vapor Pressure Deficit stress affecting grain yield. Default is False.
        results (array): Empty array for outputs.

    Returns:
        results (array): An array with estimated yield for each site.

'''
@guvectorize([(float64[:,:,:], float64[:,:], float64[:], float64[:,:] )], 
             '(m,n,o), (p, q), () -> (p,n)', nopython=True) #, forceobj=True)
def GYield_TPF( data, params, is_VPDStress, results):
    
    tn, tx, ndvi, solrad, VPDx, ipar, GPP = data
    m, n = params.shape
    for r in range(m):
        p = params[r]
        if (is_VPDStress[0]==True):
            RUE, Tmin, Toptmin, Toptmax, Tmax, tminFactor, Lvpd, Uvpd, SFvpd_Lthres, SFvpd_Uthres = p
        else:
            RUE, Tmin, Toptmin, Toptmax, Tmax, tminFactor = p
        
        tmaxFactor = 1 - tminFactor
        gpp2 = GPP.copy()
        Tday = GPP.copy()
        tpf = GPP.copy()
        sfvpd = GPP.copy()
        for j in range(tn.shape[0]):
            Tday[j] = tmaxFactor*tx[j] + tminFactor*tn[j] 
            for td in range(Tday[j].size):
                if Tday[j][td] < 0.0:
                    Tday[j][td] = 0.0
                # Calculate TPF
                tpf[j][td] = 0.0
                if ((Toptmin > Toptmax) or (Toptmax > Tmax) ):
                    tpf[j][td] = np.nan
                else:
                    if ((Tday[j][td] < Tmin) or (Tday[j][td] > Tmax)):
                        tpf[j][td] = 0.0
                    elif ((Tday[j][td] >= Toptmin) and (Tday[j][td] <= Toptmax)):
                        tpf[j][td] = 1.0
                    elif (Tday[j][td] < Toptmin):
                        x = [Tmin,Toptmin] #np.array([Tmin,Toptmin]) 
                        y = [0.0,1.0] #np.array([0.0,1.0]) 
                        #slope = getSlope(x,y)
                        #A = np.vstack([x, np.ones(len(x))]).T
                        #A = np.array(list(itertools.chain((x,np.ones(len(x)))))).T
                        #A = np.array(list(itertools.chain((x,[1.0,1.0] )))).T
                        A = [[ Tmin,  0.],[Toptmin,  1.]]
                        #slope, c = np.linalg.lstsq(A, y, rcond=None)[0]
                        #slope, c = np.linalg.lstsq(A, y)[0]
                        ave_x = (A[0][0] + A[1][0]) / len(A)  #sum([l[0] for l in A]) / len(A)  
                        ave_y = (A[0][1] + A[1][1]) / len(A)  #sum([l[1] for l in A]) / len(A)  
                        #slope = sum([l[0] * (l[1] - ave_y) for l in A]) / sum([l[0] * (l[0] - ave_x) for l in A])
                        sum_x = (A[0][0] * (A[0][0] - ave_x)) + (A[1][0] * (A[1][0] - ave_x)) 
                        sum_y = (A[0][0] * (A[0][1] - ave_y)) + (A[1][0] * (A[1][1] - ave_y))
                        slope = sum_y / sum_x
                        #b = ave_y - m * ave_x # intercept
                        tpf[j][td] = Tday[j][td] * slope
                    elif (Tday[j][td] > Toptmax):
                        x = [Toptmax,Tmax] #np.array([Toptmax,Tmax])  
                        y = [1.0,0.0] #np.array([1.0,0.0]) 
                        #slope = getSlope(x,y)
                        #A = np.vstack([x, np.ones(len(x))]).T
                        #A = np.array(list(itertools.chain((x,np.ones(len(x)))))).T
                        #A = np.array(list(itertools.chain((x,[1.0,1.0] )))).T
                        A = [[ Toptmax,  1.],[Tmax,  0.]]
                        #slope, c = np.linalg.lstsq(A, y, rcond=None)[0]
                        #slope, c = np.linalg.lstsq(A, y)[0]
                        ave_x = (A[0][0] + A[1][0]) / len(A)  #sum([l[0] for l in A]) / len(A)  
                        ave_y = (A[0][1] + A[1][1]) / len(A)  #sum([l[1] for l in A]) / len(A)  
                        #slope = sum([l[0] * (l[1] - ave_y) for l in A]) / sum([l[0] * (l[0] - ave_x) for l in A])
                        sum_x = (A[0][0] * (A[0][0] - ave_x)) + (A[1][0] * (A[1][0] - ave_x)) 
                        sum_y = (A[0][0] * (A[0][1] - ave_y)) + (A[1][0] * (A[1][1] - ave_y))
                        slope = sum_y / sum_x
                        tpf[j][td] = 1-((Tday[j][td]-Toptmax)*np.abs(slope))
                    
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
            gpp2[j] = solrad[j]*0.5*RUE*tpf[j]*ipar[j]*sfvpd[j]
            simYield = 0.0
            for val in gpp2[j][~np.isnan(gpp2[j])]:
                simYield += val
            results[r,j] = round(simYield * 0.008, 2) 

#





