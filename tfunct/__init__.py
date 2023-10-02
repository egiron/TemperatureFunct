# coding=utf-8
# Load libraries and existing datasets
from __future__ import absolute_import, division, print_function

import warnings
warnings.simplefilter('ignore')
warnings.filterwarnings('ignore')

__version__ = "tfunct version 1.0.0"
__author__ = "Azam Lashkari, Urs Christoph Schulthess, Ernesto Giron Echeverry"
__copyright__ = "Copyright (c) 2023 CIMMYT-Henan Collaborative Innovation Center"
__license__ = "Public Domain"

import gc
import numpy as np
from numba import cuda
import pandas as pd
from datetime import date, datetime
from tqdm import tqdm

from . import *
from .data import *
from .util import *
from .model import *

from .model import tday, ipar, prft, sfvpd, gpp, metrics

# 
class Site(object):
    ''' Site class

        Object containing attributes and functions related to the nursery site.

        Attributes:
            uid (integer): The unique identifier for the site.
            attributes (object): The default attributes for each location in IWIN dataset.
            params (dictionary): The parameters to use during calculations.
            pheno_dates (array): The phenology dates of the trial.
            weather (array): Table with daily weather data for each location.
            raw_ndvi (array): NDVI values during growing period.
            inputWPN (object): A dataframe with phenology, NDVI and weather data for each site.
            errors (dictionary): Dictionary with different errors presented during the model processing.
        
        Methods:
            getAttr (function): Get the attibutes of the site

    '''
    def __init__(self, uid, attributes, params=None):
        self.uid = uid
        self.attributes = attributes
        self.pheno_dates = None
        self.weather = None
        self.raw_ndvi = None
        self.inputWPN = None
        self.errors = []
        self.params = dict()
        if (params is not None):
            self.params = {**self.params, **params}
            
    def getAttr(self):
        self.attributes['errors'] = self.errors
        return self.attributes
    
    def __str__(self):
        return f"{self.uid}"
    
    
    # ---------------------------------------
    # Calculate day time temperature - TDay
    # ---------------------------------------
    def getTDay(self, m=None, tminFactor=None):
        '''Calculate day time temperature for the selected site.
        
        Parameters:
            m (object): Model with information to estimate grain yield
            tminFactor (float): Minimum Temperature factor

        Returns: 
            (array): A number or array of Day Temperature
        
        ''' 
        result = []
        if (self.weather is None):
            print("Weather data not valid")
            return
        if (tminFactor is None and m is not None):
            tminFactor = m.parameters["TMIN_PERC_FACTOR"]
        elif (tminFactor is None):
            tminFactor = self.params["TMIN_PERC_FACTOR"]
        try:
            if (('TMIN' in list(self.weather)) and ('TMAX' in list(self.weather)) ) :
                #result = tday.estimate_TDay(self.weather['TMIN'].to_numpy(), self.weather['TMAX'].to_numpy(), tminFactor )
                #result = tday._getTDay(self.weather['TMIN'].to_numpy(), self.weather['TMAX'].to_numpy(), tminFactor )
                # if (cuda.is_available()):
                result = np.zeros((self.weather['TMIN'].to_numpy().shape[0]))
                tday.tDay_gu(self.weather['TMIN'].to_numpy(), self.weather['TMAX'].to_numpy(), tminFactor, result)
                # else:
                #     inputWPN['Tdaymax'] = estimate_TDay(Tmin=self.weather['TMIN'].to_numpy(), 
                #                                         Tmax=self.weather['TMAX'].to_numpy(), 
                #                                         tminFactor=tminFactor)
                    
            else:
                print("Values for TMIN and TMAX were not found")
        except Exception as err:
            print("Error calculating Day temperature {}. Error: {}".format(self.uid, err))
            self.errors.append({"uid": self.uid, 
                                "error": "Calculating Day temperature. Error: {}".format(err)})
        
        self.inputWPN['Tdaymax'] = result
        return result
    
    
    
    # ------------------------------------------------
    # Estimate Photosynthesis reduction factor - PRFT
    # ------------------------------------------------
    def getPRFT(self, m=None, TDay=None, TOpt=None):
        ''' Estimate Photosynthesis reduction factor (PRFT) for each site.

            Parameters:
                m (object): A tfunct model
                TDay (float): Number or array of Day Temperatures
                TOpt (float): Optimum Temperature. Default value 18

            Returns: 
                (float): A number or array of PRFT

        '''
        
        if (TDay is None):
            #print("Day Temperature parameter is not valid")
            if ('Tdaymax' in self.inputWPN):
                TDay = self.inputWPN['Tdaymax'].to_numpy()
            else:
                TDay = self.getTDay(m)
        if (TOpt is None and m is not None):
            TOpt = m.parameters["Topt"]
        elif (TOpt is None):
            TOpt = self.params["Topt"]
            
        result = []
        try:
            #result = prft.calculatePRFT(TDay, TOpt )
            #result = prft._getPRFT(TDay, TOpt )
            
            result = np.zeros((TDay.shape[0]))
            prft.PRFT_gu(TDay, TOpt, result )
        except Exception as err:
            print("Error calculating photosynthesis reduction factor {}. Error: {}".format(self.uid, err))
            self.errors.append({"uid": self.uid,  
                                "error": "Calculating photosynthesis reduction factor. Error: {}"
                                .format(err)})
        
        self.inputWPN['PRFT'] = result
        return result
    
    # ------------------------------------------------
    # Calculation VPD stress
    # ------------------------------------------------
    def getSFvpd(self, m=None, VPDMAX=None, Lvpd=None, Uvpd=None, SFvpd_Lthres=None, SFvpd_Uthres=None):
        ''' Calculation of Vapor pressure deficit (VPD) stress factor for each site

            Parameters:
                m (object): A tfunct model
                VPDMAX (array): Array of daily temperature values
                Lvpd (float): A number for threshold of lower VPD. Default is 1
                Uvpd (array): A number for threshold of upper VPD. Default is 4
                SFvpd_Lthres (array): A number for threshold of stress factor of lower VPD. Default is 0.2
                SFvpd_Uthres (array): A number for threshold of stress factor of upper VPD. Default is 1
                
            Returns: 
                (array): A number or array of stressed factors of VPD

        '''
        if (VPDMAX is None):
            if ('VPDMAX' in self.inputWPN):
                VPDMAX = self.inputWPN['VPDMAX'].to_numpy()
        
        if (Lvpd is None and m is not None):
            Lvpd = m.parameters["Lvpd"]
        elif (Lvpd is None):
            Lvpd = self.params["Lvpd"]
            
        if (Uvpd is None and m is not None):
            Uvpd = m.parameters["Uvpd"]
        elif (Uvpd is None):
            Uvpd = self.params["Uvpd"]
        
        if (SFvpd_Lthres is None and m is not None):
            SFvpd_Lthres = m.parameters["SFvpd_Lthres"]
        elif (SFvpd_Lthres is None):
            SFvpd_Lthres = self.params["SFvpd_Lthres"]
            
        if (SFvpd_Uthres is None and m is not None):
            SFvpd_Uthres = m.parameters["SFvpd_Uthres"]
        elif (SFvpd_Uthres is None):
            SFvpd_Uthres = self.params["SFvpd_Uthres"]
            
        result = []
        try:
            #result = sfvpd.calculateSFvpd(VPDMAX, Lvpd, Uvpd, SFvpd_Lthres, SFvpd_Uthres)
            #result = sfvpd.getSFvpd(VPDMAX, Lvpd, Uvpd, SFvpd_Lthres, SFvpd_Uthres)
            
            result = np.zeros((VPDMAX.shape[0]))
            sfvpd.SFvpd_gu(VPDMAX, Lvpd, Uvpd, SFvpd_Lthres, SFvpd_Uthres, result)
            
            self.inputWPN['SFvpd'] = result
        except Exception as err:
            print("Error calculating SFvpd - {}. Error: {}".format(self.uid, err))
            self.errors.append({"uid": self.uid,  
                                "error": "Calculating SFvpd. Error: {}"
                                .format(err)})
        
        return result
    
    # ------------------------------------------------
    # GPP
    # ------------------------------------------------
    def getGPP(self, m=None, SolRad=None, PRFT=None, iPAR=None, RUE=3.0, stressFactor=1.0, SFvpd=None, is_VPDStress=False):
        ''' Estimate the Gross primary production (GPP) for each site

            Parameters:
                m (object): A tfunct model
                SolRad (float): Solar Radiation
                PRFT (float): Photosynthesis reduction factor
                iPAR (float): the photosynthetically active radiation (PAR) intercepted by a plant or crop
                RUE (float): Radiation-use efficiency. Default value is 3.0
                stressFactor (float): Stress Factor
                SFvpd (float): Stress Factor for VPD
                is_VPDStress (bool): Vapor pressure deficit stress. Default is `False`.

            Returns:
                (float): the gross primary production (GPP)
        '''
        if (SolRad is None):
            if ('SolRad' in self.inputWPN):
                SolRad = self.inputWPN['SolRad'].to_numpy()
        if (PRFT is None):
            if ('PRFT' in self.inputWPN):
                PRFT = self.inputWPN['PRFT'].to_numpy()
        if (iPAR is None):
            if ('iPAR' in self.inputWPN):
                iPAR = self.inputWPN['iPAR'].to_numpy()
        if (RUE is None and m is not None):
            RUE = m.parameters["RUE"]
        elif (RUE is None):
            RUE = self.params["RUE"]
            
        if (is_VPDStress is True):
            if (SFvpd is None):
                if ('SFvpd' in self.inputWPN):
                    SFvpd = self.inputWPN['SFvpd'].to_numpy()
                else:
                    SFvpd = self.getSFvpd(m)
            
        result = []
        try:
            if (is_VPDStress is True):
                #result = gpp.calculateGPP_VPDStress(SolRad, PRFT, iPAR, SFvpd, RUE)
                #result = gpp._getGPP(SolRad, PRFT, iPAR, RUE, SFvpd)
                
                result = np.zeros((SolRad.shape[0]))
                gpp.GPP_VPDStress_gu(SolRad, PRFT, iPAR, SFvpd, RUE, result )

                self.inputWPN['GPP_SFvpd'] = result
            else:
                #result = gpp.calculateGPP(SolRad, PRFT, iPAR, RUE, stressFactor)
                #result = gpp._getGPP(SolRad, PRFT, iPAR, RUE, stressFactor)
                
                result = np.zeros((SolRad.shape[0]))
                gpp.GPP_gu(SolRad, PRFT, iPAR, RUE, result )
                
                self.inputWPN['GPP'] = result
        except Exception as err:
            print("Error calculating GPP - {}. Error: {}".format(self.uid, err))
            self.errors.append({"uid": self.uid,  
                                "error": "Calculating GPP. Error: {}"
                                .format(err)})
        
        return result
    
    
    # ---------------------------------------
    # Estimate grain yield
    # ---------------------------------------
    ''' Estimate grain yield for each site
    
        Parameters:
            m (object): Model with information to estimate grain yield
            YIELD_FACTOR (float): yield factor
            is_VPDStress (bool): VPD stress. Default is False

        Returns: 
            (array): A number or array of Day Temperature
    
    ''' 
    def getGYield(self, m=None, YIELD_FACTOR=None, is_VPDStress=False):
        if (self.inputWPN is None):
            print("Input data not valid")
            return
        #if ('GPP' not in self.inputWPN):
        #    print("Input GPP values are not valid")
        #if ((is_VPDStress is True) and ('GPP_SFvpd' not in self.inputWPN)):
        #    print("Input GPP stressed VPD values are not valid")
        
        if (('GPP' not in self.inputWPN) and ('GPP_SFvpd' not in self.inputWPN)):
            print("Input GPP values are not valid")
            return
            
        if (YIELD_FACTOR is None and m is not None):
            YIELD_FACTOR = m.parameters["YIELD_FACTOR"]
        elif (YIELD_FACTOR is None):
            YIELD_FACTOR = self.params["YIELD_FACTOR"]
                
        try:
            if (is_VPDStress is True):
                if ('GPP_SFvpd' in self.inputWPN):
                    sumGPP_SFvpd = np.sum(self.inputWPN['GPP_SFvpd'].to_numpy())
                    simyield_SFvpd = round(sumGPP_SFvpd * YIELD_FACTOR, 2)
                    self.attributes['SimYield_SFvpd'] = simyield_SFvpd
            else:
                if ('GPP' in self.inputWPN):
                    sumGPP = np.sum(self.inputWPN['GPP'].to_numpy())
                    simyield = round(sumGPP * YIELD_FACTOR, 2)
                    self.attributes['SimYield'] = simyield
            
        except Exception as err:
            print("Error estimating grain yield {}. Error: {}".format(self.uid, err))
            self.errors.append({"uid": self.uid, 
                                "error": "Calculating grain yield. Error: {}".format(err)})
        
        
    # =============================
    # Fit model
    # Depreciated: Dec, 2022
    # =============================
    
    def fit(self, m=None, ft='PRFT', is_VPDStress=False, verbose=False):
        ''' Run a model to fit yield for a selected site using a specific temperature function

            Warning: Deprecated.
                This function was depreciated on Dec, 2022.
            
            Parameters:
                m (object): Model to run
                ft (str): Name of the temperature response function: Default is 'PRFT'
                is_VPDStress (bool): Vapor pressure deficit stress. Default is `False`
                verbose (bool): Display comments during processing
            
            Returns: 
                (array): An array of Sites with intermediate results
        
        '''
        if (m is None):
            print("Model parameters not valid")
            return
        #
        new_attr = {}
        try:
            if (ft=='PRFT'):
                #if ('Tdaymax' not in self.inputWPN):
                _ = self.getTDay(m)
                _ = self.getPRFT(m)
                if (is_VPDStress is True):
                    #_ = self.getSFvpd(m)
                    _ = self.getGPP(m, is_VPDStress=True)
                else:
                    _ = self.getGPP(m) # Default No stress condition
                    
                self.inputWPN['UID'] = self.uid
                # Yield
                if (is_VPDStress is True):
                    _ = self.getGYield(m, is_VPDStress=True)
                else:
                    _ = self.getGYield(m)
                
                # Select attributes to return
                new_attr["UID"] = self.attributes["UID"]
                new_attr["country"] = self.attributes["country"]
                new_attr["location"] = self.attributes["location"]
                new_attr["loc_code"] = self.attributes["loc_code"]
                new_attr["cycle"] = self.attributes["cycle"]
                new_attr["ObsYield"] = self.attributes["ObsYield"]
                
                if (is_VPDStress is True):
                    new_attr["SimYield"] = self.attributes['SimYield_SFvpd']
                else:
                    new_attr["SimYield"] = self.attributes["SimYield"]
                
                new_attr['RUE']=m.parameters['RUE']
                new_attr['Topt']=m.parameters['Topt']
                new_attr['TminFactor']= m.parameters['TMIN_PERC_FACTOR']
                new_attr['TmaxFactor']= 1 - m.parameters['TMIN_PERC_FACTOR']
                
                if (is_VPDStress is True):
                    new_attr['Lvpd']=m.parameters['Lvpd']
                    new_attr['Uvpd']=m.parameters['Uvpd']
                    new_attr['SFvpd_Lthres']=m.parameters['SFvpd_Lthres']
                    new_attr['SFvpd_Uthres']=m.parameters['SFvpd_Uthres']
                    
                # calculate stats
                
            
            elif (ft=='WEFT'):
                pass
            elif (ft=='TPF'):
                pass
            else:
                print("Temperature function not found")
            
        except:
            print("Error fitting the model in site UID:{}".format(self.uid))
            
        return new_attr #self.getAttr() #self.inputWPN #.getAttr()
        
        
    