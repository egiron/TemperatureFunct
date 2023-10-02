# coding=utf-8
from __future__ import absolute_import, division, print_function

import warnings
warnings.simplefilter('ignore')
warnings.filterwarnings('ignore')

from scipy import stats
#from scipy.stats import linregress
import numpy as np
import numba
from numba import vectorize, int32, int64, float32, float64


# ==============================================
# METRICS
# ==============================================
#
# ----------------------------------------------
# Mean absolute error (MAE)
# ----------------------------------------------
@numba.jit(nopython=True)
def MAE(obs,pred):
    #mae <- mean(abs(obs - pred))
    return np.mean(np.abs(obs - pred))

# ----------------------------------------------
# Mean Squared Error (MSE)
# ----------------------------------------------
@numba.jit(nopython=True)
def MSE(obs,pred):
    #mse = mean((obs - pred)^2)
    return np.mean((obs - pred)**2)

# ----------------------------------------------
# Root Mean Square Error (RMSE) 
# ----------------------------------------------
@numba.jit(nopython=True)
def RMSE(obs,pred):
    #rmse = sqrt(mean((obs - pred)^2))
    return np.sqrt(np.mean((obs - pred)**2))

# ----------------------------------------------
# Mean Absolute Percentage Error (MAPE)
# ----------------------------------------------
@numba.jit(nopython=True)
def MAPE(obs,pred):
    #mape <- mean(abs((obs - pred)/obs))*100
    return np.mean(np.abs((obs - pred)/obs))*100

# ----------------------------------------------
# R-Squared
# ----------------------------------------------
@numba.jit(nopython=True)
def R2(obs,pred):
    #r2 = 1 - (sum((pred - obs)^2) / sum((obs - mean(obs))^2))
    #acor = np.corrcoef(obs,pred)
    #acor = acor[len(acor)/2:]# use only second half
    return (np.corrcoef(obs,pred)[0, 1])**2

# ----------------------------------------------
# Root Mean Square Relative Error (RMSRE)
# ----------------------------------------------
@numba.jit(nopython=True)
def RMSRE(obs,pred):
    #rmsre = 100 * sqrt(mean(((obs - pred) / obs)^2))
    return 100 * np.sqrt(np.mean(((obs - pred) / obs)**2))

# ----------------------------------------------
# Nash–Sutcliffe model efficiency (EF)
# ----------------------------------------------
# EF is a distance measure that compares model MSE with the MSE of using 
# the average of measured values as an estimator. Therefore, EF is useful 
# for making statements about the skill of a model relative to this simple 
# reference estimator. For a model that simulates perfectly,
# EF = 1, and for a model that has the same squared error of simulation as 
# the mean of the measurements, EF = 0. EF is positive for a model that has 
# a smaller squared error than the mean of the measurements.
@numba.jit(nopython=True)
def EF(obs,pred):
    #ef = 1 - ( sum((obs - pred)^2) / sum((obs - mean(obs))^2) )
    return 1 - ( np.sum((obs - pred)**2) / np.sum((obs - np.mean(obs))**2) )


# ----------------------------------------------
# Pearson Correlation Coefficient
# ----------------------------------------------
''' Pearson Correlation Coefficient '''
def pearson_CC(x,y):
    sxy = np.sum((x - x.mean())*(y - y.mean())) / x.shape[0]
    rho = sxy / (np.std(x)*np.std(y))
    return rho

# ----------------------------------------------
# Lin’s Concordance Correlation Coefficient (CCC)
# Computes Lin's (1989, 2000) concordance correlation coefficient for 
# agreement on a continuous measure obtained by two methods. The 
# concordance correlation coefficient combines measures of both precision 
# and accuracy to determine how far the observed data deviate from the 
# line of perfect concordance (that is, the line at 45 degrees on a square 
# scatter plot).
# ----------------------------------------------
''' Concordance Correlation Coefficient '''
def _ccc(x,y):
    sxy = np.sum((x - x.mean())*(y - y.mean())) / x.shape[0]
    rhoc = 2*sxy / (np.var(x) + np.var(y) + (x.mean() - y.mean())**2)
    return rhoc

'''Lin's Concordance correlation coefficient'''
@numba.jit(parallel=True)
def CCC(y_true, y_pred):
    # covariance between y_true and y_pred
    s_xy = np.cov([y_true, y_pred])[0,1]
    # means
    x_m = np.mean(y_true)
    y_m = np.mean(y_pred)
    # variances
    s_x_sq = np.var(y_true)
    s_y_sq = np.var(y_pred)
    
    # condordance correlation coefficient
    ccc = (2.0*s_xy) / (s_x_sq + s_y_sq + (x_m-y_m)**2)
    
    return ccc

'''
    Variable C.b is a bias correction factor that measures how far the best-fit line deviates 
    from a line at 45 degrees (a measure of accuracy). 
    
    No deviation from the 45 degree line occurs when C.b = 1. See Lin (1989 page 258).
'''
@numba.jit(parallel=True)
def Cb(x,y):
    k = len(y)
    yb = np.mean(y)
    sy2 = np.var(y) * (k - 1) / k
    sd1 = np.std(y)
    #print(k, yb, sy2, sd1)
    xb = np.mean(x)
    sx2 = np.var(x) * (k - 1) / k
    sd2 = np.std(x)
    
    r = np.corrcoef(x, y)[0,1] ## same as pearson CC
    sl = r * sd1 / sd2
    sxy = r * np.sqrt(sx2 * sy2)
    p = 2 * sxy / (sx2 + sy2 + (yb - xb)**2)
    
    # The following taken from the Stata code for function "concord" (changed 290408):
    bcf = p / r
    return bcf


# ----------------------------------------------
@numba.jit(parallel=True) #nopython=True
def _calculateMetrics(obs, pred):
    m_mae = MAE(obs,pred)
    m_mse = MSE(obs,pred)
    m_rmse = RMSE(obs,pred)
    m_rmsre = RMSRE(obs,pred)
    m_R2 = R2(obs,pred)
    m_EF = EF(obs,pred)
    m_MAPE = MAPE(obs,pred)
    m_accuracy = 100 - m_MAPE
    m_Cb = Cb(obs,pred)
    m_CCC = CCC(obs,pred)
    #Create the linear regression
    #A = np.vstack([obs, np.ones(len(obs))]).T
    #a, res, rank, s = np.linalg.lstsq(A, pred, rcond=None)
    #m_R2a = (1-res/(obs.size*np.var(pred,axis=0)))[0]
    #print(a, res, rank, s)
    #m_slope, m_intercept = np.linalg.lstsq(A, pred, rcond=None)[0]
    m_slope, m_intercept, r_value, m_pvalue, std_err = stats.linregress(obs,pred)
    #print(slope, intercept, r_value, p_value, std_err)
    #m_R2a, m_pvalue = pearsonr(obs,pred)
    #print("STDERR:", std_err, "\n", "MAE:", m_mae, "\n", "MSE:", m_mse, "\n",
    #     "RMSE:", m_rmse, "\n", "RMSRE:", m_rmsre, "\n",
    #     "MAPE:", m_MAPE, "\n",
    #     "p-value:", m_pvalue, "\n", 
    #     "R-squared:", m_R2, "\n",
    #     #"Adj-R-squared:", m_AdjR2, "\n", 
    #      "EF:", m_EF, "\n",
    #     "Accuracy: ", m_accuracy)

    return [m_mae, m_mse, m_rmse, m_rmsre, m_MAPE, m_pvalue, 
              m_R2, m_EF, m_intercept, m_slope, std_err, m_Cb, m_CCC, m_accuracy]
#

@numba.jit(parallel=True)
def calculateMetrics(obs, pred):
    nd = 3
    m_mae = round(MAE(obs,pred), nd)
    m_mse = round(MSE(obs,pred), nd)
    m_rmse = round(RMSE(obs,pred), nd)
    m_rmsre = round(RMSRE(obs,pred), nd)
    m_R2 = round(R2(obs,pred), nd)
    m_EF = round(EF(obs,pred), nd)
    m_MAPE = round(MAPE(obs,pred), nd)
    m_accuracy = round(100 - m_MAPE, 2)
    m_Cb = round(Cb(obs,pred), 2)
    m_CCC = round(CCC(obs,pred), 2)
    
    m_slope, m_intercept, r_value, m_pvalue, std_err = stats.linregress(obs,pred)

    return [m_mae, m_mse, m_rmse, m_rmsre, m_MAPE, m_pvalue, 
              m_R2, m_EF, round(m_intercept, 4), round(m_slope, 4), m_Cb, m_CCC, m_accuracy]




