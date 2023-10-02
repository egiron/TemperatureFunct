# coding=utf-8
from __future__ import absolute_import, division, print_function

import warnings
warnings.simplefilter('ignore')
warnings.filterwarnings('ignore')

import sys, os, gc
import itertools
import pandas as pd
from datetime import datetime
from joblib import Parallel, delayed
from tqdm import tqdm

from numba import cuda
import pyarrow as pa
import pyarrow.parquet as pq

from . import *
from ..data import *
from ..util import *
from ..util import mergeParquetFiles

from . import tday, ipar, prft, gpp, sfvpd, gyield, metrics

sys.path.insert(0, r"../tfunct")
import tfunct


class Model(object):
    __version__ = "tfunct version 1.0.0"
    __author__ = "Azam Lashkari, Urs Christoph Schulthess, Ernesto Giron Echeverry"
    __copyright__ = "Copyright (c) 2022 CIMMYT-Henan Collaborative Innovation Center"
    __license__ = "Public Domain"
    
    
    def __init__(self, config, params=None):
        self.config = config
        self.parameters = dict(
                RUE = 3,
                DRYMATTER = 0.8,
                FACTOR_TON_HA = 0.01,
                YIELD_FACTOR = 0.8 * 0.01,
                TMIN_PERC_FACTOR = 0.25,
                NDVI_lowerThreshold = 0.16,
                Toptmin = 15,
                Topt = 18,
                Toptmax = 25,
                Tmin = 9,
                Tmax = 34,
                Lvpd = 1,
                Uvpd = 4,
                SFvpd_Lthres = 0.2,
                SFvpd_Uthres = 1,
            )
        if (params is not None):
            self.parameters = {**self.parameters, **params}
            
    def load_raw_datasets(self):
        ''' Load raw phenology and AgERA5 datasets 
        
            Return:
                A existing dataset including WeatherFile, PhenoFile, NDVIFile.
                The raw data can be seen in config object. eg. config['PhenoFile']
        
        '''
        if (self.config is None):
            print("Configuration is not valid")
            return
        
        WeatherFile, PhenoFile, NDVIFile = None, None, None
        try:
            # ---------------------
            # Read Phenology
            # ---------------------
            sites_phenology_path = os.path.join(self.config['PHENO_FILE'])
            if (os.path.exists(sites_phenology_path)):
                # Load Phenology file
                PhenoFile = pd.read_csv(sites_phenology_path)
                # convert sowing date to DATE format
                PhenoFile['Sowing_date']=pd.to_datetime(PhenoFile['SowingDateQC'].astype(str), format='%Y-%m-%d')
                # convert days to heading and days to maturity to date (DATE), sowing date = 0
                PhenoFile['Heading_date'] = PhenoFile[['Sowing_date', 'Days_To_Heading']]\
                .apply(lambda row: row['Sowing_date'] + pd.DateOffset(days=row['Days_To_Heading']), axis=1)
                PhenoFile['Maturity_date'] = PhenoFile[['Sowing_date', 'Days_To_Maturity']]\
                .apply(lambda row: row['Sowing_date'] + pd.DateOffset(days=row['Days_To_Maturity']), axis=1)
                # remove cols
                PhenoFile.drop(columns=['SowingDateQC'], inplace=True)
                PhenoFile.reset_index(drop=True, inplace=True)
                PhenoFile['UID'] = PhenoFile.index + 1
            else:
                raise ("Error reading crop phenology file")
                
            # ---------------------
            # Read NDVI
            # ---------------------
            sites_ndvi_path = os.path.join(self.config['NDVI_FILE'])
            if (os.path.exists(sites_ndvi_path)):
                # Load NDVI file
                NDVIFile = pd.read_csv(sites_ndvi_path)
                NDVIFile['phenotype_date']=pd.to_datetime(NDVIFile['phenotype_date'].astype(str), format='%Y-%m-%d')
            else:
                raise ("Error reading NDVI file")
                
            # ---------------------
            # Read Weather
            # ---------------------
            sites_weather_path = os.path.join(self.config['WEATHER_FILE'])
            if (os.path.exists(sites_weather_path)):
                # Load Weather file
                WeatherFile = pd.read_csv(sites_weather_path)
                # add Date column to WeatherFile
                WeatherFile['Date'] = pd.to_datetime(
                    WeatherFile['Year'].astype(str) + '-' + WeatherFile['Month'].astype(str) + '-' +
                    WeatherFile['Day'].astype(str), 
                    format='%Y-%m-%d')
                
                WeatherFile.rename(columns={
                    'Shortwave Radiation [MJ/m2/d]': 'SolRad',
                    'TemperatureMax [C]':'TMAX', 
                    'TemperatureMin [C]':'TMIN',
                    #'Canopy Temperature [C]':'TC',
                    'Precipitation [mm/d]':'PCP', 
                    'Vapor Pressure Deficit max [kPa]':'VPDMAX',
                    #'Vapor Pressure Deficit at Tc [kPa]':'VPDTC', 
                    'Wind Speed 2m [m/s]':'WINDSPEED',
                    'Relative Humidity max [%]':'RhumMax',
                    'Relative Humidity min [%]':'RhumMin', 
                }, inplace=True)
            else:
                raise ("Error reading weather file")
        
        except Exception as err:
            print("Problem reading raw datasets. Error: {}".format(err))
            
        self.config['WeatherFile'] = WeatherFile
        self.config['PhenoFile'] = PhenoFile
        self.config['NDVIFile'] = NDVIFile

    # ---------------------------------------
    # Prepare dataset
    # ---------------------------------------
    def preprocess_raw_datasets(self, data):
        ''' Preprocess raw datasets 

            Parameters:
                data (dictionary): The phenology, ndvi and weather data for each location in example dataset.

            Return:
                A dataset in a specific format for the tfunct package
        '''
        if (self.config is None):
            print("Configuration is not valid")
            return
        
        WeatherFile, PhenoFile, NDVIFile = None, None, None
        try:
            try:
                # ---------------------
                # Process Phenology
                # ---------------------
                PhenoFile = data['Pheno'].copy()
                # convert days to heading and days to maturity to date (DATE), sowing date = 0
                PhenoFile['Sowing_date']= PhenoFile['SowingDateQC']
                PhenoFile['Heading_date'] = PhenoFile[['Sowing_date', 'Days_To_Heading']]\
                .apply(lambda row: row['Sowing_date'] + pd.DateOffset(days=row['Days_To_Heading']), axis=1)
                PhenoFile['Maturity_date'] = PhenoFile[['Sowing_date', 'Days_To_Maturity']]\
                .apply(lambda row: row['Sowing_date'] + pd.DateOffset(days=row['Days_To_Maturity']), axis=1)
                # remove cols
                PhenoFile.drop(columns=['SowingDateQC'], inplace=True)
                PhenoFile.reset_index(drop=True, inplace=True)
                PhenoFile['UID'] = PhenoFile.index + 1
            except Exception as err:
                raise ("Error reading crop phenology file", err)
                
            # ---------------------
            # Process NDVI
            # ---------------------
            NDVIFile = data['NDVI'].copy()
           
            # ---------------------
            # Process Weather
            # ---------------------
            try:
                WeatherFile = data['Weather'].copy()
                # add Date column to WeatherFile
                WeatherFile['Date'] = pd.to_datetime(
                    WeatherFile['Year'].astype(str) + '-' + WeatherFile['Month'].astype(str) + '-' +
                    WeatherFile['Day'].astype(str), format='%Y-%m-%d')
                    
                WeatherFile.rename(columns={
                    'Shortwave Radiation [MJ/m2/d]': 'SolRad',
                    'TemperatureMax [C]':'TMAX', 
                    'TemperatureMin [C]':'TMIN',
                    #'Canopy Temperature [C]':'TC',
                    'Precipitation [mm/d]':'PCP', 
                    'Vapor Pressure Deficit max [kPa]':'VPDMAX',
                    #'Vapor Pressure Deficit at Tc [kPa]':'VPDTC', 
                    'Wind Speed 2m [m/s]':'WINDSPEED',
                    'Relative Humidity max [%]':'RhumMax',
                    'Relative Humidity min [%]':'RhumMin', 
                }, inplace=True)
            except Exception as err:
                raise ("Error reading weather file", err)
        
        except Exception as err:
            print("Problem reading raw datasets. Error: {}".format(err))
            
        self.config['WeatherFile'] = WeatherFile
        self.config['PhenoFile'] = PhenoFile
        self.config['NDVIFile'] = NDVIFile
    
    # ---------------------------------------
    # Prepare dataset
    # ---------------------------------------
    '''
        Prepare data for further analysis
    '''
    def prepareData(self):
        if (self.config['PhenoFile'] is None):
            print("Phenology file not found")
            return
        if (self.config['NDVIFile'] is None):
            print("NDVI file not found")
            return
        if (self.config['WeatherFile'] is None):
            print("Weather file not found")
            return
        
        sites = []
        for i, row in tqdm(self.config['PhenoFile'].iterrows()):
            #s = dict(row.iteritems())  # iteritems was removed in 2.0.0
            s = dict(row.items())
            s['Sowing_date'] = str(s['Sowing_date']).split(' ')[0] #pd.to_datetime(str(s['Sowing_date']), format='%Y-%m-%d')
            s['Heading_date'] = str(s['Heading_date']).split(' ')[0]
            s['Maturity_date'] = str(s['Maturity_date']).split(' ')[0]

            site = tfunct.Site(uid=s['UID'], attributes=s)
            # Get phenology dates between Heading and Maturity 
            site.pheno_dates = pd.date_range(start=s['Heading_date'], end=s['Maturity_date']) #, inclusive="both")
            # Agregar NDVI
            ndvi = self.config['NDVIFile'][((self.config['NDVIFile']['loc_code']==s['loc_code']) 
                                            & (self.config['NDVIFile']['cycle']==s['cycle'])) ][['phenotype_date', 'NDVI']]
            site.raw_ndvi = ndvi.to_dict()

            # Add interpolated NDVI
            # fill missing dates, or populate missing dates
            df = ndvi.sort_values(['phenotype_date']).set_index('phenotype_date')
            df = df.resample('1D').mean().interpolate(method='linear').ffill().bfill() #method='nearest', axis=0 limit_area=None, limit_direction='forward', axis=0
            df.drop_duplicates(inplace=True)
            pheno_ndvi = pd.merge(pd.DataFrame(site.pheno_dates).rename(columns={0:'phenotype_date'}),
                                  df, on=['phenotype_date'], how='left')
            # Get last row: if NDVI at maturity is null, set the NDVI at last row NDVI = 0.16
            pheno_ndvi.reset_index(drop=True, inplace=True)
            if (pheno_ndvi.iloc[-1][['NDVI']].isna().any()):
                pheno_ndvi.loc[len(pheno_ndvi)-1, 'NDVI'] = self.parameters['NDVI_lowerThreshold']
            # fill missing dates, or populate missing dates
            df = pheno_ndvi.sort_values(['phenotype_date']).set_index('phenotype_date')
            df = df.resample('1D').mean().interpolate(method='linear').ffill().bfill() 
            df.reset_index(drop=False, inplace=True)
            df.drop_duplicates(inplace=True)
            # store interp NDVI
            #site.ndvi = df #.to_dict()

            # add Weather data
            _mask = ((self.config['WeatherFile']['Date'] >= s['Heading_date']) 
                     & (self.config['WeatherFile']['Date'] <= s['Maturity_date']) 
                     & (self.config['WeatherFile']['location']==s['location']) )
            weather = self.config['WeatherFile'][_mask].reset_index(drop=True)
            site.weather = weather
            #weather = self.config['WeatherFile'][self.config['WeatherFile']['location']==s['location']]
            inputWPN = pd.merge(df, weather, left_on=['phenotype_date'], right_on=['Date'], how='left')
            # Add some useful columns
            inputWPN["loc_code"] = s['loc_code']
            inputWPN["cycle"] = s['cycle']
            inputWPN["lat"] = s['lat'] # Replace in case 
            inputWPN["lon"] = s['lon']
            inputWPN["DOY"] = inputWPN['phenotype_date'].apply(lambda x: pd.Timestamp(x).dayofyear)
            inputWPN['DOY'] = inputWPN['DOY'].astype(int)

            # Calculate iPAR
            #inputWPN["iPAR"] = ipar._getiPAR(inputWPN["NDVI"].to_numpy())
            # if (cuda.is_available()):
            ipar_ndvi = np.zeros((inputWPN["NDVI"].to_numpy().shape[0]))
            ipar.iPAR_gu(inputWPN["NDVI"].to_numpy(), ipar_ndvi)
            inputWPN["iPAR"] = ipar_ndvi
            #            
            res_Tdaymax = np.zeros((inputWPN['TMIN'].to_numpy().shape[0]))
            tday.tDay_gu(inputWPN['TMIN'].to_numpy(),inputWPN['TMAX'].to_numpy(), 
                            self.parameters['TMIN_PERC_FACTOR'], res_Tdaymax) 
            inputWPN['Tdaymax'] = res_Tdaymax
            # else:
            #     inputWPN["iPAR"] = self.getIPAR(m=self, ndvi=inputWPN["NDVI"].to_numpy()) # In parallel for large datasets
            #     #inputWPN['Tdaymax'] = tday._getTDay(inputWPN['TMIN'].to_numpy(), inputWPN['TMAX'].to_numpy(), 0.25) 
            #     #inputWPN['Tdaymax'] = self.getTDay(self, inputWPN) # In parallel for large datasets
            #     #
            #     #inputWPN['Tdaymax'] = self.getTDay(m=self, w=weather, tminFactor=self.parameters['TMIN_PERC_FACTOR'])
                
            #     inputWPN['Tdaymax'] = tday.estimate_TDay(Tmin=inputWPN['TMIN'].to_numpy(), 
            #                                              Tmax=inputWPN['TMAX'].to_numpy(), 
            #                                              tminFactor=self.parameters['TMIN_PERC_FACTOR'])
                
            # Remove some colums
            inputWPN.drop(columns=['Year', 'Month', 'Day', 'Date','RhumMax', 'RhumMin', 'WINDSPEED'], inplace=True)
            site.inputWPN = inputWPN
            
            # number of days with TMAX > 34 C and number of days with TMIN < 9 C
            site.attributes['ndays_tmn_lt9'] = len(inputWPN[inputWPN['TMIN'] < 9])
            site.attributes['ndays_tmx_gt34'] = len(inputWPN[inputWPN['TMAX'] > 34])
            # Average of NDVI, iPAR, etc from heading to maturity
            site.attributes['avg_Tdaymax'] = round(np.nanmean(inputWPN['Tdaymax']),3)
            site.attributes['avg_NDVI'] = round(np.nanmean(inputWPN['NDVI']),3)
            site.attributes['avg_iPAR'] = round(np.nanmean(inputWPN['iPAR']),3)
            
            sites.append(site)
        
        return sites
    
    # ---------------------------------------
    # Calculate day time temperature - TDay
    # ---------------------------------------
    def getTDay(self, m=None, w=None, tminFactor=None):
        '''Calculate day time temperature

        Parameters:
            m (object): A tfunct model with sites and the necessary information
            w (array): Table of weather data with minimum and maximum temperature records
            tminFactor (float): Minimum Temperature factor. Default is usually 0.25. 
                                It can be reviewed in configuration parameters such as `TMIN_PERC_FACTOR`.

        Returns: 
            (array): A number or array of Day Temperature

        ''' 
        result = []
        if (w is None):
            #w = self.config['WeatherFile']
            #if (w is None):
            print("Table of weather data with Minimum and Maximum Temperatures is not valid")
            return
        if (tminFactor is None and m is not None):
            tminFactor = self.parameters["TMIN_PERC_FACTOR"]
        elif (tminFactor is None):
            tminFactor = self.parameters["TMIN_PERC_FACTOR"]
        try:
            if (('TMIN' in list(w)) and ('TMAX' in list(w)) ) :
                result = tday.estimate_TDay(w['TMIN'].to_numpy(), w['TMAX'].to_numpy(), tminFactor )
            else:
                print("Values for TMIN and TMAX were not found")
        except Exception as err:
            print("Error calculating Day temperature. Error: {}".format(err))

        return result
    
    
    # ---------------------------------------
    # Total light interception - iPAR
    # ---------------------------------------
    def getIPAR(self, m=None, ndvi=None):
        ''' Total light interception - iPAR

            Reference:
                iPAR = NDVI * 1.25 - 0.19 # between heading and maturity (Campos et al. 2018)
                
                iPAR = NDVI * 1.25 - 0.21 Daughtry et al. (1992)


                - Asrar, G., Fuchs, M., Kanemasu, E.T., Hatfield, J.L., 1984. 
                Estimating absorbed photosynthetic radiation and leaf area index from spectral reflectance 
                in wheat. Agron. J. 76, 30-306.

                - Campos 2018 Remote sensing-based crop biomass with water or light-driven crop growth models in 
                wheat commercial fields.

            Parameters:
                m (object): A tfunct model 
                ndvi (array): Array of float values

            Returns: 
               (array): An array of Total light interception values

        '''
        result = []
        if (ndvi is None):
            print("NDVI data is not valid")
            return
        try:
            result = ipar.estimate_IPAR(ndvi)
        except Exception as err:
            print("Error calculating Total light interception - iPAR. Error: {}".format(err))

        return result
    
    
    #
    # -------------------------------------------------
    # Combinations for grain yield - PRFT, WETF, TPF
    # -------------------------------------------------
    ''' Get yield using a defined temperature function and default parameters
        
        Warning: Deprecated
            Stop using this function.

        Parameters:
            tfun (str): Name of the temperature function
            sites (dictionary): List of sites with the basic information to estimate yield.
            is_VPDStress (bool): Vapor Pressure Deficit stress. Default is False.
            expCSV (bool): Export grain yield results in an individual file. Default is True.

        Returns:
            (array): A dataframe with simulated grain yield values

    '''
    def getYield_v0(self, tfun='PRFT', sites=None, is_VPDStress=False, expCSV=True):
        
        arr_tn, arr_tx, arr_ndvi, arr_solrad, arr_vpdmax, arr_ipar, GPP = gpp.prepareVectors_dataset(sites)
        data_input = np.stack((arr_tn, arr_tx, arr_ndvi, arr_solrad, arr_vpdmax, arr_ipar, GPP))

        if (tfun=='PRFT'):
            # if (cuda.is_available()):
            gpp.GPP(arr_tn, arr_tx, arr_solrad, arr_vpdmax, arr_ipar,
                [self.parameters['TMIN_PERC_FACTOR']], [self.parameters['Topt']], [self.parameters['RUE']], 
                [is_VPDStress], [self.parameters['Lvpd']], [self.parameters['Uvpd']], 
                [self.parameters['SFvpd_Lthres']], [self.parameters['SFvpd_Uthres']], GPP)

            for obs in range(0, len(GPP)):
                sites[obs].inputWPN['GPP'] = GPP[obs][~np.isnan(GPP[obs])]

            attr = []
            for _id in range(0, len(sites)):
                sumGPP = np.sum(sites[_id].inputWPN['GPP'].to_numpy())
                sites[_id].attributes['SimYield'] = round(sumGPP * self.parameters['YIELD_FACTOR'], 2)
                attr.append(sites[_id].attributes)

            #print(sites[19].attributes)
            df_GYield = pd.DataFrame(attr)[['UID', 'country', 'location', 'loc_code', 'cycle', 'ObsYield',
                   'SimYield']]

            RUE = str(self.parameters['RUE'])
            Topt = str(self.parameters['Topt'])
            TminFactor = str(self.parameters['TMIN_PERC_FACTOR'])
            TmaxFactor = str(1 - self.parameters['TMIN_PERC_FACTOR'])
            Lvpd = str(self.parameters['Lvpd'])
            Uvpd = str(self.parameters['Uvpd'])
            SFvpd_Lthres = str(self.parameters['SFvpd_Lthres'])
            SFvpd_Uthres = str(self.parameters['SFvpd_Uthres'])

            if (is_VPDStress is True): 
                    Lvpd = str(self.parameters['Lvpd'])
                    Uvpd = str(self.parameters['Uvpd'])
                    SFvpd_Lthres = str(self.parameters['SFvpd_Lthres'])
                    SFvpd_Uthres = str(self.parameters['SFvpd_Uthres'])

            df_GYield['RUE'] = RUE
            df_GYield['Topt'] = Topt
            df_GYield['TminFactor'] = TminFactor
            df_GYield['TmaxFactor'] = TmaxFactor

            if (is_VPDStress is True): 
                df_GYield['Lvpd'] = Lvpd
                df_GYield['Uvpd'] = Uvpd
                df_GYield['SFvpd_Lthres'] = SFvpd_Lthres
                df_GYield['SFvpd_Uthres'] = SFvpd_Uthres

            # Save grain yield model in individual files
            if (expCSV is True):
                TFname = "PRFT_noStress"
                fname = "Yield_RUE" + RUE + "_Topt" + Topt + "_TmnFact" + TminFactor + "_" +TFname
                if (is_VPDStress is True):
                    TFname = "PRFT_SFvpd"
                    fname = "Yield_RUE" + RUE + "_Topt" + Topt + "_L" + Lvpd + "_U" + Uvpd + "_SFvpd" \
                    + SFvpd_Lthres + "-" + SFvpd_Uthres + "_TmnFact" + TminFactor + "_" +TFname

                res_path = os.path.join(self.config['RESULTS_PATH'],"PRFT",TFname)
                if not os.path.isdir(res_path):
                    os.makedirs(res_path)
                df_GYield.to_csv(os.path.join(res_path,'{}.csv'.format(fname)), index=False)
        #
        elif (tfun=='WETF'):
            attr = []
            for _id in range(0, len(sites)):
                attr.append(sites[_id].attributes)

            df_GYield = pd.DataFrame(attr)[['UID', 'country', 'location', 'loc_code', 'cycle', 'ObsYield']]
            df_GYield['ObsYield'] = df_GYield['ObsYield'].astype(float).round(2)

            if (is_VPDStress==False): # No stress conditions
                isVPDStress = [False]
                #array_params_to_run = np.array([[3.0, 9.0, 18.0, 34.0, 0.25]]) #RUE, Tmin, Topt, Tmax, tminFactor
                array_params_to_run = np.array([[self.parameters['RUE'], 
                                                 self.parameters['Tmin'], self.parameters['Topt'], 
                                                 self.parameters['Tmax'], self.parameters['TMIN_PERC_FACTOR']] ])
                rows = array_params_to_run.shape[0]
                cols = arr_tn.shape[0]
                array_results = np.zeros(rows*cols, dtype=np.float64).reshape(rows, cols)
                #gpp.GYield_WETF(data_input, array_params_to_run, isVPDStress, array_results)
                gyield.estimate(data_input, array_params_to_run, is_VPDStress, [2], array_results)
                df_GYield = pd.DataFrame(gpp.createDF(df_GYield.to_numpy(), array_params_to_run, array_results))
                df_GYield.columns = ['UID', 'country', 'location', 'loc_code', 'cycle', 'ObsYield', 
                                         'RUE','Tmin', 'Topt','Tmax', 'TminFactor','SimYield']


            else:
                isVPDStress = [True]
                #RUE, Tmin, Topt, Tmax, tminFactor, Lvpd, Uvpd, SFvpd_Lthres, SFvpd_Uthres
                #array_params_to_run = np.array([[3.0, 9.0, 18.0, 34.0, 0.25, 1.0, 4.0, 0.2, 1.0]])
                array_params_to_run = np.array([[self.parameters['RUE'], 
                                                 self.parameters['Tmin'], self.parameters['Topt'], 
                                                 self.parameters['Tmax'], self.parameters['TMIN_PERC_FACTOR'],
                                                 self.parameters['Lvpd'], self.parameters['Uvpd'], 
                                                 self.parameters['SFvpd_Lthres'], self.parameters['SFvpd_Uthres']]
                                               ])
                rows = array_params_to_run.shape[0]
                cols = arr_tn.shape[0]
                array_results = np.zeros(rows*cols, dtype=np.float64).reshape(rows, cols)
                #gpp.GYield_WETF(data_input, array_params_to_run, isVPDStress, array_results)
                gyield.estimate(data_input, array_params_to_run, is_VPDStress, [2], array_results)
                df_GYield = pd.DataFrame(gpp.createDF(df_GYield.to_numpy(), array_params_to_run, array_results))
                df_GYield.columns = ['UID', 'country', 'location', 'loc_code', 'cycle', 'ObsYield', 
                                      'RUE','Tmin', 'Topt','Tmax', 'TminFactor','Lvpd','Uvpd',
                                      'SFvpd_Lthres','SFvpd_Uthres','SimYield']
            # 
            # Save grain yield model in individual files
            if (expCSV is True):
                RUE = str(self.parameters['RUE'])
                Tmin = str(self.parameters['Tmin'])
                Topt = str(self.parameters['Topt'])
                Tmax = str(self.parameters['Tmax'])
                TminFactor = str(self.parameters['TMIN_PERC_FACTOR'])
                TmaxFactor = str(1 - self.parameters['TMIN_PERC_FACTOR'])
                Lvpd = str(self.parameters['Lvpd'])
                Uvpd = str(self.parameters['Uvpd'])
                SFvpd_Lthres = str(self.parameters['SFvpd_Lthres'])
                SFvpd_Uthres = str(self.parameters['SFvpd_Uthres'])
                TFname = "WETF_noStress"
                fname = "Yield_RUE" + RUE + "_Tmin" + Tmin + "_Topt" + Topt + "_Tmax" + Tmax + "_TmnFact" \
                + TminFactor + "_" +TFname
                if (is_VPDStress is True):
                    TFname = "WETF_SFvpd"
                    fname = "Yield_RUE" + RUE + "_Tmin" + Tmin + "_Topt" + Topt + "_Tmax" + Tmax + "_L" \
                    + Lvpd + "_U" + Uvpd + "_SFvpd" + SFvpd_Lthres + "-" + SFvpd_Uthres + "_TmnFact" \
                    + TminFactor + "_" +TFname

                res_path = os.path.join(self.config['RESULTS_PATH'],"WETF",TFname)
                if not os.path.isdir(res_path):
                    os.makedirs(res_path)
                df_GYield.to_csv(os.path.join(res_path,'{}.csv'.format(fname)), index=False)
                #df_GYield.to_parquet(os.path.join(res_path,'{}.parquet'.format(fname)), index=False)


        elif (tfun=='TPF'):
            pass

        return df_GYield
    
    # --------------------- Updated GYield Sep 30, 2023 -----------------------
    '''
        Set up all the required parameters to run analisis using global variables in the model
        This parameters are organize for a single observation of each parameter or model run.
    '''
    def setupParamters(self, tfun='PRFT', sites=None, is_VPDStress=False):
        arr_tn, arr_tx, arr_ndvi, arr_solrad, arr_vpdmax, arr_ipar, GPP = gpp.prepareVectors_dataset(sites)
        data_input = np.stack((arr_tn, arr_tx, arr_ndvi, arr_solrad, arr_vpdmax, arr_ipar, GPP))
        temfun = [1] if tfun=='PRFT' else [2] if tfun=='WETF' else [3] if tfun=='TPF' else None

        attr = []
        for _id in range(0, len(sites)):
            attr.append(sites[_id].attributes)
        df_GYield = pd.DataFrame(attr)[['UID', 'country', 'location', 'loc_code', 'cycle', 'ObsYield']]
        df_GYield['ObsYield'] = df_GYield['ObsYield'].astype(float).round(2)

        RUE = [self.parameters['RUE']]
        Tmin = [self.parameters['Tmin']]
        Toptmin = [self.parameters['Toptmin']]
        Topt = [self.parameters['Topt']]
        Toptmax = [self.parameters['Toptmax']]
        Tmax = [self.parameters['Tmax']]
        TminFactor = [self.parameters['TMIN_PERC_FACTOR']]
        Lvpd = [self.parameters['Lvpd']]
        Uvpd = [self.parameters['Uvpd']]
        SFvpd_Lthres = [self.parameters['SFvpd_Lthres']]
        SFvpd_Uthres = [self.parameters['SFvpd_Uthres']]

        # No stress conditions
        cols = arr_tn.shape[0]
        array_params_to_run, array_results = self.getCombinations(functype=tfun, cols=cols, RUE=RUE, 
                                                                   Tmin=Tmin, Toptmin=Toptmin, Topt=Topt, 
                                                                   Toptmax=Toptmax, 
                                                                   Tmax=Tmax, TminFactor=TminFactor, 
                                                                   Lvpd=Lvpd, Uvpd=Uvpd, SFvpd_Lthres=SFvpd_Lthres, 
                                                                   SFvpd_Uthres=SFvpd_Uthres,
                                                                   isVPDStress=is_VPDStress)

        rows = array_params_to_run.shape[0]
        array_results = np.zeros(rows*cols, dtype=np.float64).reshape(rows, cols)

        TFname = f"{tfun}_noStress" if is_VPDStress is False else f"{tfun}_SFvpd"
        if (tfun=='PRFT'):
            fname = "Yield_RUE" + str(RUE[0]) + "_Topt" + str(Topt[0]) + "_TmnFact" + str(TminFactor[0]) + "_" +TFname
            if (is_VPDStress is True):
                fname = "Yield_RUE" + str(RUE[0]) + "_Topt" + str(Topt[0]) + "_L" + str(Lvpd[0]) + "_U" + str(Uvpd[0]) + "_SFvpd" \
                + str(SFvpd_Lthres[0]) + "-" + str(SFvpd_Uthres[0]) + "_TmnFact" + str(TminFactor[0]) + "_" +TFname
        elif (tfun=='WETF'):
            fname = "Yield_RUE" + str(RUE[0]) + "_Tmin" + str(Tmin[0]) + "_Topt" + str(Topt[0]) + "_Tmax" + str(Tmax[0]) + "_TmnFact" + str(TminFactor[0]) + "_" +TFname
            if (is_VPDStress is True):
                fname = "Yield_RUE" + str(RUE[0]) + "_Tmin" + str(Tmin[0]) + "_Topt" + str(Topt[0]) + "_Tmax" + str(Tmax[0]) + "_L" + str(Lvpd[0]) + "_U" + str(Uvpd[0]) + "_SFvpd" \
                + str(SFvpd_Lthres[0]) + "-" + str(SFvpd_Uthres[0]) + "_TmnFact" + str(TminFactor[0]) + "_" +TFname
        elif (tfun=='TPF'):
            fname = "Yield_RUE" + str(RUE[0]) +  "_Tmin" + str(Tmin[0]) + "_Toptmin" + str(Toptmin[0]) + "_Toptmax" + str(Toptmax[0]) + "_Tmax" + str(Tmax[0]) + "_TmnFact" + str(TminFactor[0]) + "_" +TFname
            if (is_VPDStress is True):
                fname = "Yield_RUE" + str(RUE[0]) + "_Tmin" + str(Tmin[0]) + "_Toptmin" + str(Toptmin[0]) + "_Toptmax" + str(Toptmax[0]) + "_Tmax" + str(Tmax[0]) + "_L" + str(Lvpd[0]) + "_U" + str(Uvpd[0]) + "_SFvpd" \
                + str(SFvpd_Lthres[0]) + "-" + str(SFvpd_Uthres[0]) + "_TmnFact" + str(TminFactor[0]) + "_" +TFname

        return df_GYield, data_input, array_params_to_run, temfun, array_results, TFname, fname

    def getYield(self, tfun='PRFT', sites=None, is_VPDStress=False, expCSV=True):
        ''' Get yield using a defined temperature function and default parameters

            Parameters:
                tfun (str): Name of the temperature function
                sites (dictionary): List of sites with the basic information to estimate yield.
                is_VPDStress (bool): Vapor Pressure Deficit stress. Default is False.
                expCSV (bool): Export grain yield results in an individual file. Default is True.

            Returns:
                (array): A dataframe with simulated grain yield values

        '''
        df_GYield, data_input, array_params_to_run, temfun, array_results, TFname, fname = \
        self.setupParamters(tfun=tfun, sites=sites, is_VPDStress=is_VPDStress)
        # Estimate grain yield
        gyield.estimate(data_input, array_params_to_run, [is_VPDStress], temfun, array_results)
        df_GYield = pd.DataFrame(gpp.createDF(df_GYield.to_numpy(), array_params_to_run, array_results))
        if (tfun=='PRFT'):
            if (is_VPDStress is False):
                df_GYield.columns = ['UID', 'country', 'location', 'loc_code', 'cycle', 'ObsYield', 
                                     'RUE','TminFactor','Topt','SimYield']
            else:
                df_GYield.columns = ['UID', 'country', 'location', 'loc_code', 'cycle', 'ObsYield', 
                                     'RUE','TminFactor','Topt', 'Lvpd','Uvpd', 'SFvpd_Lthres','SFvpd_Uthres', 
                                     'SimYield']
            
            # gpp.estimate(data_input, array_params_to_run, [is_VPDStress], temfun, array_results)
            #for obs in range(0, len(GPP)):
            #    sites[obs].inputWPN['GPP'] = GPP[obs][~np.isnan(GPP[obs])]
            #attr = []
            #for _id in range(0, len(sites)):
            #    sumGPP = np.sum(sites[_id].inputWPN['GPP'].to_numpy())
            #    sites[_id].attributes['SimYield'] = round(sumGPP * self.parameters['YIELD_FACTOR'], 2)
            #    attr.append(sites[_id].attributes)

        #
        elif (tfun=='WETF'):
            #gyield.estimate(data_input, array_params_to_run, [is_VPDStress], temfun, array_results)
            #df_GYield = pd.DataFrame(gpp.createDF(df_GYield.to_numpy(), array_params_to_run, array_results))
            if (is_VPDStress is False):
                df_GYield.columns = ['UID', 'country', 'location', 'loc_code', 'cycle', 'ObsYield', 
                                     'RUE','Tmin', 'Topt', 'Tmax', 'TminFactor','SimYield']
            else:
                df_GYield.columns = ['UID', 'country', 'location', 'loc_code', 'cycle', 'ObsYield', 
                                     'RUE','TminFactor','Tmin', 'Topt', 'Tmax', 'Lvpd','Uvpd', 
                                     'SFvpd_Lthres','SFvpd_Uthres', 'SimYield']

        elif (tfun=='TPF'):
            #gyield.estimate(data_input, array_params_to_run, [is_VPDStress], temfun, array_results)
            #df_GYield = pd.DataFrame(gpp.createDF(df_GYield.to_numpy(), array_params_to_run, array_results))
            if (is_VPDStress is False):
                df_GYield.columns = ['UID', 'country', 'location', 'loc_code', 'cycle', 'ObsYield', 
                                     'RUE','Tmin', 'Toptmin', 'Toptmax', 'Tmax', 'TminFactor','SimYield']
            else:
                df_GYield.columns = ['UID', 'country', 'location', 'loc_code', 'cycle', 'ObsYield', 
                                     'RUE','TminFactor','Tmin', 'Toptmin', 'Toptmax', 'Tmax', 'Lvpd','Uvpd', 
                                     'SFvpd_Lthres','SFvpd_Uthres', 'SimYield']

        # Save grain yield model in individual files
        if (expCSV is True):
            res_path = os.path.join(self.config['RESULTS_PATH'],tfun,TFname)
            if not os.path.isdir(res_path):
                os.makedirs(res_path)
            df_GYield.to_csv(os.path.join(res_path,'{}.csv'.format(fname)), index=False)

        return df_GYield
    
    
    
    # ---------------------------
    # Combinations for grain yield - PRFT
    ''' Get yield using **PRFT** function

        Warning: Deprecated
            Stop using this function.

        Parameters:
            sites (dictionary): List of sites with the basic information to estimate yield.
            is_VPDStress (bool): Vapor Pressure Deficit stress. Default is False.
            expCSV (bool): Export grain yield results in an individual file. Default is True.

        Returns:
            (array): A dataframe with simulated grain yield values

    '''
    def getYield_PRFT(self, sites=None, is_VPDStress=False, expCSV=True):
        
        arr_tn, arr_tx, arr_ndvi, arr_solrad, arr_vpdmax, arr_ipar, output_array = gpp.prepareVectors_dataset(sites)
        
        # if (cuda.is_available()):
        gpp.GPP(arr_tn, arr_tx, arr_solrad, arr_vpdmax, arr_ipar,
            [self.parameters['TMIN_PERC_FACTOR']], [self.parameters['Topt']], [self.parameters['RUE']], 
            [is_VPDStress], [self.parameters['Lvpd']], [self.parameters['Uvpd']], 
            [self.parameters['SFvpd_Lthres']], [self.parameters['SFvpd_Uthres']], output_array)
    
        for obs in range(0, len(output_array)):
            sites[obs].inputWPN['GPP'] = output_array[obs][~np.isnan(output_array[obs])]

        attr = []
        for _id in range(0, len(sites)):
            sumGPP = np.sum(sites[_id].inputWPN['GPP'].to_numpy())
            sites[_id].attributes['SimYield'] = round(sumGPP * self.parameters['YIELD_FACTOR'], 2)
            attr.append(sites[_id].attributes)

        #print(sites[19].attributes)
        df_GYield = pd.DataFrame(attr)[['UID', 'country', 'location', 'loc_code', 'cycle', 'ObsYield',
               'SimYield']]

        RUE = str(self.parameters['RUE'])
        Topt = str(self.parameters['Topt'])
        TminFactor = str(self.parameters['TMIN_PERC_FACTOR'])
        TmaxFactor = str(1 - self.parameters['TMIN_PERC_FACTOR'])
        Lvpd = str(self.parameters['Lvpd'])
        Uvpd = str(self.parameters['Uvpd'])
        SFvpd_Lthres = str(self.parameters['SFvpd_Lthres'])
        SFvpd_Uthres = str(self.parameters['SFvpd_Uthres'])

        if (is_VPDStress is True): 
                Lvpd = str(self.parameters['Lvpd'])
                Uvpd = str(self.parameters['Uvpd'])
                SFvpd_Lthres = str(self.parameters['SFvpd_Lthres'])
                SFvpd_Uthres = str(self.parameters['SFvpd_Uthres'])

        df_GYield['RUE'] = RUE
        df_GYield['Topt'] = Topt
        df_GYield['TminFactor'] = TminFactor
        df_GYield['TmaxFactor'] = TmaxFactor

        if (is_VPDStress is True): 
            df_GYield['Lvpd'] = Lvpd
            df_GYield['Uvpd'] = Uvpd
            df_GYield['SFvpd_Lthres'] = SFvpd_Lthres
            df_GYield['SFvpd_Uthres'] = SFvpd_Uthres

        # Save grain yield model in individual files
        if (expCSV is True):
            TFname = "PRFT_noStress"
            fname = "Yield_RUE" + RUE + "_Topt" + Topt + "_TmnFact" + TminFactor + "_" +TFname
            if (is_VPDStress is True):
                TFname = "PRFT_SFvpd"
                fname = "Yield_RUE" + RUE + "_Topt" + Topt + "_L" + Lvpd + "_U" + Uvpd + "_SFvpd" \
                + SFvpd_Lthres + "-" + SFvpd_Uthres + "_TmnFact" + TminFactor + "_" +TFname

            res_path = os.path.join(self.config['RESULTS_PATH'],"PRFT",TFname)
            if not os.path.isdir(res_path):
                os.makedirs(res_path)
            df_GYield.to_csv(os.path.join(res_path,'{}.csv'.format(fname)), index=False)


        return df_GYield
    
    
    
    # ---------------------------
    # Combinations for PRFT 
    # ---------------------------
    def runCombinations_PRFT(self, data=None, sites=None, comb=None, is_VPDStress=False, 
                             expCSV=True, individualCSV=False, target=None):
        ''' Run several combinations using different parameters for PRFT model 
        
            Warning: Deprecated
                    Stop using this function.

            Parameters:
                data (object): Array of base data.
                sites (object): Array of sites
                comb (object): Array of combinations
                is_VPDStress (bool): Vapor Pressure Deficit stress. Default is False
                expCSV (bool): Save results in CSV format. Default is True
                individualCSV (bool): Default is False
                target (str): Name of the target device. Available values are `cpu`, `parallel` and `cuda`. Default is `cpu`
        
            Returns:
                (object): A dataframe with estimated grain yield of all combinations.
        
        '''
        outputs = []
        if (data is None):
            arr_tn, arr_tx, arr_ndvi, arr_solrad, arr_vpdmax, arr_ipar, output_array = gpp.prepareVectors_dataset(sites)
        else:
            arr_tn, arr_tx, arr_ndvi, arr_solrad, arr_vpdmax, arr_ipar, output_array = data
        
        if (comb is None):
            comb = [dict(
                RUE = 3,
                TMIN_PERC_FACTOR = 0.25,
                Topt = 18,
                Lvpd = 1,
                Uvpd = 4,
                SFvpd_Lthres = 0.2,
                SFvpd_Uthres = 1,
            )]
        for p in tqdm(comb):
            self.parameters = {**self.parameters, **p}
            if (target=='cuda'):
                rows, cols = arr_tn.shape
                output_array = np.ones(rows*cols, dtype=np.float64).reshape(rows, cols)

                # invoke on CUDA with manually managed memory
                dev_tn = cuda.to_device(arr_tn) # alloc and copy input data
                dev_tx = cuda.to_device(arr_tx) # alloc and copy input data
                dev_ndvi = cuda.to_device(arr_ndvi)
                dev_solrad = cuda.to_device(arr_solrad)
                dev_vpdmax = cuda.to_device(arr_vpdmax)
                #dev_ipar = cuda.to_device(arr_ipar)
                
                dev_tminFactor = cuda.to_device([float(self.parameters['TMIN_PERC_FACTOR'])])
                dev_Topt = cuda.to_device([float(self.parameters['Topt'])])
                dev_RUE = cuda.to_device([float(self.parameters['RUE'])])
                dev_is_VPDStress = cuda.to_device([0.0])
                dev_Lvpd = cuda.to_device([float(self.parameters['Lvpd'])])
                dev_Uvpd = cuda.to_device([float(self.parameters['Uvpd'])])
                dev_SFvpd_Lthres = cuda.to_device([float(self.parameters['SFvpd_Lthres'])])
                dev_SFvpd_Uthres = cuda.to_device([float(self.parameters['SFvpd_Uthres'])])

                dev_res = cuda.to_device(output_array) 

                gpp.GPP_cuda(dev_tn, dev_tx, dev_ndvi, arr_solrad, arr_vpdmax, dev_tminFactor, 
                         dev_Topt, dev_RUE, dev_is_VPDStress, dev_Lvpd, dev_Uvpd,
                         dev_SFvpd_Lthres, dev_SFvpd_Uthres, dev_res)

                dev_res.copy_to_host(output_array) # retrieve the result
            elif (target=='parallel'):
                gpp.GPP_parallel(arr_tn, arr_tx, arr_ndvi, arr_solrad, arr_vpdmax, 
                    [self.parameters['TMIN_PERC_FACTOR']], [self.parameters['Topt']], [self.parameters['RUE']], 
                    [is_VPDStress], [self.parameters['Lvpd']], [self.parameters['Uvpd']], 
                    [self.parameters['SFvpd_Lthres']], [self.parameters['SFvpd_Uthres']], output_array)
            elif (target=='cpu'):
                gpp.GPP_cpu(arr_tn, arr_tx, arr_ndvi, arr_solrad, arr_vpdmax, 
                    [self.parameters['TMIN_PERC_FACTOR']], [self.parameters['Topt']], [self.parameters['RUE']], 
                    [is_VPDStress], [self.parameters['Lvpd']], [self.parameters['Uvpd']], 
                    [self.parameters['SFvpd_Lthres']], [self.parameters['SFvpd_Uthres']], output_array)
            else:
                # Try to send ready ipar data but with no changes in speed up
                gpp.GPP(arr_tn, arr_tx, arr_solrad, arr_vpdmax, arr_ipar,
                    [self.parameters['TMIN_PERC_FACTOR']], [self.parameters['Topt']], [self.parameters['RUE']], 
                    [is_VPDStress], [self.parameters['Lvpd']], [self.parameters['Uvpd']], 
                    [self.parameters['SFvpd_Lthres']], [self.parameters['SFvpd_Uthres']], output_array)

            for obs in range(0, len(output_array)):
                sites[obs].inputWPN['GPP'] = output_array[obs][~np.isnan(output_array[obs])]

            attr = []
            for _id in range(0, len(sites)):
                sumGPP = np.sum(sites[_id].inputWPN['GPP'].to_numpy())
                sites[_id].attributes['SimYield'] = round(sumGPP * self.parameters['YIELD_FACTOR'], 2)
                attr.append(sites[_id].attributes)

            #print(sites[19].attributes)
            df_GYield = pd.DataFrame(attr)[['UID', 'country', 'location', 'loc_code', 'cycle', 'ObsYield', 'SimYield']]

            RUE = str(self.parameters['RUE'])
            Topt = str(self.parameters['Topt'])
            TminFactor = str(self.parameters['TMIN_PERC_FACTOR'])
            TmaxFactor = str(1 - self.parameters['TMIN_PERC_FACTOR'])
            
            if (is_VPDStress is True): 
                Lvpd = str(self.parameters['Lvpd'])
                Uvpd = str(self.parameters['Uvpd'])
                SFvpd_Lthres = str(self.parameters['SFvpd_Lthres'])
                SFvpd_Uthres = str(self.parameters['SFvpd_Uthres'])

            df_GYield['RUE'] = RUE
            df_GYield['Topt'] = Topt
            df_GYield['TminFactor'] = TminFactor
            df_GYield['TmaxFactor'] = TmaxFactor
            
            if (is_VPDStress is True): 
                df_GYield['Lvpd'] = Lvpd
                df_GYield['Uvpd'] = Uvpd
                df_GYield['SFvpd_Lthres'] = SFvpd_Lthres
                df_GYield['SFvpd_Uthres'] = SFvpd_Uthres
            
            # 
            outputs.append(df_GYield)

            # Save grain yield model in individual files
            if (individualCSV is True):
                TFname = "PRFT_noStress"
                fname = "Yield_RUE" + RUE + "_Topt" + Topt + "_TmnFact" + TminFactor + "_" +TFname
                if (is_VPDStress is True):
                    TFname = "PRFT_SFvpd"
                    fname = "Yield_RUE" + RUE + "_Topt" + Topt + "_L" + Lvpd + "_U" + Uvpd + "_SFvpd" \
                    + SFvpd_Lthres + "-" + SFvpd_Uthres + "_TmnFact" + TminFactor + "_" +TFname

                res_path = os.path.join(self.config['RESULTS_PATH'],"PRFT",TFname)
                if not os.path.isdir(res_path):
                    os.makedirs(res_path)
                df_GYield.to_csv(os.path.join(res_path,'{}.csv'.format(fname)), index=False)

        del arr_tn, arr_tx, arr_ndvi, arr_solrad, arr_vpdmax, output_array, attr, df_GYield
        _ = gc.collect()
        
        df = pd.concat([x for x in outputs ])
        if (expCSV is True):
            TFname = "PRFT_noStress"
            fname = "combinations_Yield_PRFT_noStress"
            if (is_VPDStress is True):
                TFname = "PRFT_SFvpd"
                fname = "combinations_Yield_PRFT_SFvpd"
            res_path = os.path.join(self.config['RESULTS_PATH'],"PRFT",TFname)
            if not os.path.isdir(res_path):
                os.makedirs(res_path)
            df.to_csv(os.path.join(res_path,'{}.csv'.format(fname)), index=False)
            
        return df
    
    # =============================
    # Fit model (Optimized)
    # Depreciated: Dec, 2022
    # =============================
    def fit(self, data=None, sites=None, comb=None, ft='PRFT', is_VPDStress=False, expCSV=True, 
            individualCSV=False, n_jobs=4, verbose=False):
        '''
            Run a model to fit yield

            Warning: Deprecated on Dec, 2022.
                    Stop using this function.
            
            Parameters:
                data (object): Array of base data.
                sites (object): Array of sites
                comb (object): Array of combinations
                ft (str): Name of the temperature response function. Default is 'PRFT'. 
                is_VPDStress (bool): Vapor Pressure Deficit stress. Default is False
                expCSV (bool): Save results in CSV format. Default is True
                individualCSV (bool): Default is False
                n_jobs (int): Number of CPU cores to use in paralell processing. Default is 4.
                verbose (bool): Display comments during processing. Default is False
        
            Returns:
                (array): An array of sites with intermediate results
        
        '''
        if (sites is None):
            print("Input parameters are not valid")
            return
        #
        df_GYield = None
        try:
            if (ft=='PRFT'):
                
                # Run in parallel
                if (comb is not None):
                    prt = int(len(comb)/n_jobs)
                    splitted_runs = [comb[:prt], comb[prt:prt*2], comb[prt*2:prt*3], comb[prt*3:]]
                    with Parallel(n_jobs=n_jobs, verbose=5) as parallel:
                        delayed_funcs = [delayed(lambda p: self.runCombinations_PRFT(data=data, sites=sites,
                                                comb=p, is_VPDStress=is_VPDStress, expCSV=expCSV, individualCSV=individualCSV)
                                                )(run) for run in splitted_runs]
                        output = parallel(delayed_funcs)

                    df_GYield = pd.concat([x for x in output ])
                else:
                    df_GYield = self.runCombinations_PRFT(data=data, sites=sites, comb=comb, 
                                                      is_VPDStress=is_VPDStress, expCSV=expCSV, individualCSV=individualCSV)
                    
                
                # calculate stats
                
                #
                
            
            elif (ft=='WETF'):
                pass
            elif (ft=='TPF'):
                pass
            else:
                print("Temperature function not found")
            
        except:
            print("Error fitting the model")
            
        return df_GYield
    
    
    # ---------------------------------------------
    # Export WETF results arrow file in parquet format
    # ---------------------------------------------
    ''' Create arrays (columns) as an input for the exported arrow file in parquet format '''
    def export_WETF_Results(self, df_gy, batch, array_results, fname="combinations_Yield_WETF_SFvpd", TFname="WETF_SFvpd"):
        col_UID = []
        col_country = []
        col_location = []
        col_loc_code = []
        col_cycle = []
        col_ObsYield = []
        col_RUE = []
        col_Tmin = []
        col_Topt = []
        col_Tmax = []
        col_TminFactor = []
        col_Lvpd = []
        col_Uvpd = []
        col_SFvpd_Lthres = []
        col_SFvpd_Uthres = []
        col_SimYield = []
        for i in range(0, len(batch)):
            for j in range(len(df_gy)):
                #reg = list(df_gy[j]) + list(array_params_to_run[i]) + list([array_results[i][j]])
                #UID, country, location, loc_code, cycle, ObsYield, RUE, Tmin, Topt, Tmax, TminFactor, \
                # Lvpd, Uvpd, SFvpd_Lthres, SFvpd_Uthres, SimYield = reg
                UID, country, location, loc_code, cycle, ObsYield = df_gy[j]
                RUE, Tmin, Topt, Tmax, TminFactor, Lvpd, Uvpd, SFvpd_Lthres, SFvpd_Uthres = batch[i]
                SimYield = array_results[i][j]

                col_UID.append(UID)
                col_country.append(country)
                col_location.append(location)
                col_loc_code.append(loc_code)
                col_cycle.append(cycle)
                col_ObsYield.append(ObsYield)
                col_RUE.append(RUE)
                col_Tmin.append(Tmin)
                col_Topt.append(Topt)
                col_Tmax.append(Tmax)
                col_TminFactor.append(TminFactor)
                col_Lvpd.append(Lvpd)
                col_Uvpd.append(Uvpd)
                col_SFvpd_Lthres.append(SFvpd_Lthres)
                col_SFvpd_Uthres.append(SFvpd_Uthres)
                col_SimYield.append(SimYield)
        # create export file
        ndarray_table = pa.table(
            {
                "UID": col_UID,
                "country": col_country,
                "location": col_location,
                "loc_code": col_loc_code,
                "cycle": col_cycle,
                "ObsYield": col_ObsYield,
                "RUE": col_RUE,
                "Tmin": col_Tmin,
                "Topt": col_Topt,
                "Tmax": col_Tmax,
                "TminFactor": col_TminFactor,
                "Lvpd": col_Lvpd,
                "Uvpd": col_Uvpd,
                "SFvpd_Lthres": col_SFvpd_Lthres,
                "SFvpd_Uthres": col_SFvpd_Uthres,
                "SimYield": col_SimYield
            }
        )
        res_path = os.path.join(self.config['RESULTS_PATH'],"WETF",TFname)
        if not os.path.isdir(res_path):
            os.makedirs(res_path)
        pq.write_table(ndarray_table, os.path.join(res_path,'{}.parquet'.format(fname)) )

        del col_UID, col_country, col_location, col_loc_code, col_cycle, col_ObsYield
        del col_RUE, col_Tmin, col_Topt, col_Tmax, col_TminFactor, col_Lvpd, col_Uvpd
        del col_SFvpd_Lthres, col_SFvpd_Uthres, col_SimYield
        _ = gc.collect()
        
        
    # ---------------------------------------------
    # Export WETF results arrow file in parquet format
    # ---------------------------------------------
    def exportResults_in_batch_WETF_SFvpd(self, df_gy, comb, res, batch_size=4, merge=False, mergedfname='merged', 
                                          TFname="WETF_SFvpd_parts", removeParts=False):
        assert len(comb) == len(res)
        prt = int(len(comb)/batch_size)
        for i in tqdm(range(1, batch_size+1)):
            #print(i, len(comb[prt*(i-1):prt*i]) )
            self.export_WETF_Results(df_gy, comb[prt*(i-1):prt*i], res[prt*(i-1):prt*i], 
                              fname="combinations_Yield_WETF_SFvpd_part_{}".format(i), TFname=TFname)
        # Merge results
        if (merge is True):
            res_path = os.path.join(self.config['RESULTS_PATH'],"WETF",TFname)
            out_path = os.path.join(self.config['RESULTS_PATH'],"WETF","WETF_SFvpd")
            mergeParquetFiles(res_path, out_path, fname=mergedfname, removeParts=removeParts)
    
    # ---------------------------------------------
    # Export TPF results arrow file in parquet format
    # ---------------------------------------------
    ''' Create arrays (columns) as an input for the exported arrow file in parquet format '''
    def export_TPF_Results(self, df_gy, batch, array_results, fname="combinations_Yield_TPF_",
                           TFname="TPF_", is_VPDStress=False):
        col_UID = []
        col_country = []
        col_location = []
        col_loc_code = []
        col_cycle = []
        col_ObsYield = []
        col_RUE = []
        col_Tmin = []
        col_Toptmin = []
        col_Toptmax = []
        col_Tmax = []
        col_TminFactor = []
        col_Lvpd = []
        col_Uvpd = []
        col_SFvpd_Lthres = []
        col_SFvpd_Uthres = []
        col_SimYield = []
        for i in range(0, len(batch)):
            for j in range(len(df_gy)):
                UID, country, location, loc_code, cycle, ObsYield = df_gy[j]
                if (is_VPDStress is True):
                    RUE, Tmin, Toptmin, Toptmax, Tmax, TminFactor, Lvpd, Uvpd, SFvpd_Lthres, SFvpd_Uthres = batch[i]
                else:
                    RUE, Tmin, Toptmin, Toptmax, Tmax, TminFactor = batch[i]
                    
                SimYield = array_results[i][j]

                col_UID.append(UID)
                col_country.append(country)
                col_location.append(location)
                col_loc_code.append(loc_code)
                col_cycle.append(cycle)
                col_ObsYield.append(ObsYield)
                col_RUE.append(RUE)
                col_Tmin.append(Tmin)
                col_Toptmin.append(Toptmin)
                col_Toptmax.append(Toptmax)
                col_Tmax.append(Tmax)
                col_TminFactor.append(TminFactor)
                
                if (is_VPDStress is True):
                    col_Lvpd.append(Lvpd)
                    col_Uvpd.append(Uvpd)
                    col_SFvpd_Lthres.append(SFvpd_Lthres)
                    col_SFvpd_Uthres.append(SFvpd_Uthres)
                #
                col_SimYield.append(SimYield)
        # create export file
        ndarray_table = None
        #fname = "combinations_Yield_TPF_noStress"
        #TFname = "TPF_noStress"
        if (is_VPDStress is True):
            ndarray_table = pa.table(
                {
                    "UID": col_UID,
                    "country": col_country,
                    "location": col_location,
                    "loc_code": col_loc_code,
                    "cycle": col_cycle,
                    "ObsYield": col_ObsYield,
                    "RUE": col_RUE,
                    "Tmin": col_Tmin,
                    "Toptmin": col_Toptmin,
                    "Toptmax": col_Toptmax,
                    "Tmax": col_Tmax,
                    "TminFactor": col_TminFactor,
                    "Lvpd": col_Lvpd,
                    "Uvpd": col_Uvpd,
                    "SFvpd_Lthres": col_SFvpd_Lthres,
                    "SFvpd_Uthres": col_SFvpd_Uthres,
                    "SimYield": col_SimYield
                }
            )
        else:
            ndarray_table = pa.table(
                {
                    "UID": col_UID,
                    "country": col_country,
                    "location": col_location,
                    "loc_code": col_loc_code,
                    "cycle": col_cycle,
                    "ObsYield": col_ObsYield,
                    "RUE": col_RUE,
                    "Tmin": col_Tmin,
                    "Toptmin": col_Toptmin,
                    "Toptmax": col_Toptmax,
                    "Tmax": col_Tmax,
                    "TminFactor": col_TminFactor,
                    "SimYield": col_SimYield
                }
            )
            
        res_path = os.path.join(self.config['RESULTS_PATH'],"TPF",TFname)
        if not os.path.isdir(res_path):
            os.makedirs(res_path)
        pq.write_table(ndarray_table, os.path.join(res_path,'{}.parquet'.format(fname)) )

        del col_UID, col_country, col_location, col_loc_code, col_cycle, col_ObsYield
        del col_RUE, col_Tmin, col_Toptmin, col_Toptmax, col_Tmax, col_TminFactor, col_Lvpd, col_Uvpd
        del col_SFvpd_Lthres, col_SFvpd_Uthres, col_SimYield
        _ = gc.collect()
    
    ''' Create arrays (columns) as an input for the exported arrow file in parquet format '''
    def export_TPF_Results_and_Metrics(self, df_gy, batch, array_results, fname="combinations_Yield_TPF_",
                           TFname="TPF_", is_VPDStress=False, calc_metrics=True):
        
        col_UID = []
        col_country = []
        col_location = []
        col_loc_code = []
        col_cycle = []
        col_ObsYield = []
        col_RUE = []
        col_Tmin = []
        col_Toptmin = []
        col_Toptmax = []
        col_Tmax = []
        col_TminFactor = []
        col_Lvpd = []
        col_Uvpd = []
        col_SFvpd_Lthres = []
        col_SFvpd_Uthres = []
        col_SimYield = []
        # Columns to create Metrics file in parquet format
        colm_RUE = []
        colm_Tmin = []
        colm_Toptmin = []
        colm_Toptmax = []
        colm_Tmax = []
        colm_TminFactor = []
        colm_Lvpd = []
        colm_Uvpd = []
        colm_SFvpd_Lthres = []
        colm_SFvpd_Uthres = []
        colm_mae = [] 
        colm_mse = [] 
        colm_rmse = [] 
        colm_rmsre = [] 
        colm_MAPE = [] 
        colm_pvalue = [] 
        colm_R2 = [] 
        colm_EF = [] 
        colm_intercept = [] 
        colm_slope = [] 
        colm_Cb = [] 
        colm_CCC = [] 
        colm_accuracy = []
        val_stats_ObsYield = [x for x in df_gy[:,5]] # Select only grain yield column
        for i in range(0, len(batch)):
            val_stats_SimYield = []
            if (is_VPDStress is True):
                RUE, Tmin, Toptmin, Toptmax, Tmax, TminFactor, Lvpd, Uvpd, SFvpd_Lthres, SFvpd_Uthres = batch[i]
            else:
                RUE, Tmin, Toptmin, Toptmax, Tmax, TminFactor = batch[i]
            for j in range(len(df_gy)):

                UID, country, location, loc_code, cycle, ObsYield = df_gy[j]

                SimYield = array_results[i][j]
                val_stats_SimYield.append(array_results[i][j]) # for metrics

                col_UID.append(UID)
                col_country.append(country)
                col_location.append(location)
                col_loc_code.append(loc_code)
                col_cycle.append(cycle)
                col_ObsYield.append(ObsYield)
                col_RUE.append(RUE)
                col_Tmin.append(Tmin)
                col_Toptmin.append(Toptmin)
                col_Toptmax.append(Toptmax)
                col_Tmax.append(Tmax)
                col_TminFactor.append(TminFactor)
                col_SimYield.append(SimYield)

                if (is_VPDStress is True):
                    col_Lvpd.append(Lvpd)
                    col_Uvpd.append(Uvpd)
                    col_SFvpd_Lthres.append(SFvpd_Lthres)
                    col_SFvpd_Uthres.append(SFvpd_Uthres)
                #
            # Metrics
            if (calc_metrics is True):
                #m = metrics.calculateMetrics(np.array(val_stats_ObsYield), np.array(val_stats_SimYield))
                #st.append(list(batch[i]) + list(m))
                m_mae, m_mse, m_rmse, m_rmsre, m_MAPE, m_pvalue, m_R2, m_EF, m_intercept, m_slope, m_Cb, m_CCC, m_accuracy = \
                metrics.calculateMetrics(np.array(val_stats_ObsYield), np.array(val_stats_SimYield))
                colm_RUE.append(RUE)
                colm_Tmin.append(Tmin)
                colm_Toptmin.append(Toptmin)
                colm_Toptmax.append(Toptmax)
                colm_Tmax.append(Tmax)
                colm_TminFactor.append(TminFactor)
                if (is_VPDStress is True):
                    colm_Lvpd.append(Lvpd)
                    colm_Uvpd.append(Uvpd)
                    colm_SFvpd_Lthres.append(SFvpd_Lthres)
                    colm_SFvpd_Uthres.append(SFvpd_Uthres)
                colm_mae.append(m_mae)
                colm_mse.append(m_mse) 
                colm_rmse.append(m_rmse)
                colm_rmsre.append(m_rmsre)
                colm_MAPE.append(m_MAPE)
                colm_pvalue.append(m_pvalue)
                colm_R2.append(m_R2)
                colm_EF.append(m_EF) 
                colm_intercept.append(m_intercept)
                colm_slope.append(m_slope)
                colm_Cb.append(m_Cb)
                colm_CCC.append(m_CCC)
                colm_accuracy.append(m_accuracy)

        # create export files
        ndarray_table = None
        ndarray_metrics = None
        if (is_VPDStress is True):
            ndarray_table = pa.table(
                {
                    "UID": col_UID,
                    "country": col_country,
                    "location": col_location,
                    "loc_code": col_loc_code,
                    "cycle": col_cycle,
                    "ObsYield": col_ObsYield,
                    "RUE": col_RUE,
                    "Tmin": col_Tmin,
                    "Toptmin": col_Toptmin,
                    "Toptmax": col_Toptmax,
                    "Tmax": col_Tmax,
                    "TminFactor": col_TminFactor,
                    "Lvpd": col_Lvpd,
                    "Uvpd": col_Uvpd,
                    "SFvpd_Lthres": col_SFvpd_Lthres,
                    "SFvpd_Uthres": col_SFvpd_Uthres,
                    "SimYield": col_SimYield
                }
            )
            # metrics
            if (calc_metrics is True):
                ndarray_metrics = pa.table(
                    {
                        "RUE": colm_RUE,
                        "Tmin": colm_Tmin,
                        "Toptmin": colm_Toptmin,
                        "Toptmax": colm_Toptmax,
                        "Tmax": colm_Tmax,
                        "TminFactor": colm_TminFactor,
                        "Lvpd": colm_Lvpd,
                        "Uvpd": colm_Uvpd,
                        "SFvpd_Lthres": colm_SFvpd_Lthres,
                        "SFvpd_Uthres": colm_SFvpd_Uthres,
                        "MAE": colm_mae,
                        "MSE": colm_mse,
                        "RMSE": colm_rmse,
                        "RMSRE": colm_rmsre,
                        "MAPE": colm_MAPE,
                        "pvalue": colm_pvalue,
                        "R2": colm_R2,
                        "EF": colm_EF,
                        "intercept": colm_intercept,
                        "slope": colm_slope,
                        "Cb": colm_Cb,
                        "CCC": colm_CCC,
                        "Accuracy": colm_accuracy 
                    }
                )
        else:
            ndarray_table = pa.table(
                {
                    "UID": col_UID,
                    "country": col_country,
                    "location": col_location,
                    "loc_code": col_loc_code,
                    "cycle": col_cycle,
                    "ObsYield": col_ObsYield,
                    "RUE": col_RUE,
                    "Tmin": col_Tmin,
                    "Toptmin": col_Toptmin,
                    "Toptmax": col_Toptmax,
                    "Tmax": col_Tmax,
                    "TminFactor": col_TminFactor,
                    "SimYield": col_SimYield
                }
            )
            # metrics
            if (calc_metrics is True):
                ndarray_metrics = pa.table(
                    {
                        "RUE": colm_RUE,
                        "Tmin": colm_Tmin,
                        "Toptmin": colm_Toptmin,
                        "Toptmax": colm_Toptmax,
                        "Tmax": colm_Tmax,
                        "TminFactor": colm_TminFactor,
                        "MAE": colm_mae,
                        "MSE": colm_mse,
                        "RMSE": colm_rmse,
                        "RMSRE": colm_rmsre,
                        "MAPE": colm_MAPE,
                        "pvalue": colm_pvalue,
                        "R2": colm_R2,
                        "EF": colm_EF,
                        "intercept": colm_intercept,
                        "slope": colm_slope,
                        "Cb": colm_Cb,
                        "CCC": colm_CCC,
                        "Accuracy": colm_accuracy 
                    }
                )

        res_path = os.path.join(self.config['RESULTS_PATH'],"TPF",TFname)
        if not os.path.isdir(res_path):
            os.makedirs(res_path)
        #print("Saving results...")
        pq.write_table(ndarray_table, os.path.join(res_path,'{}.parquet'.format(fname)) )
        # saving metrics
        if (calc_metrics is True):
            #print("Saving metrics...")
            res_path = os.path.join(self.config['RESULTS_PATH'],"TPF","{}_metrics".format(TFname))
            if not os.path.isdir(res_path):
                os.makedirs(res_path)
            pq.write_table(ndarray_metrics, os.path.join(res_path,'metrics_{}.parquet'.format(fname)) )
        #
        del col_UID, col_country, col_location, col_loc_code, col_cycle, col_ObsYield
        del col_RUE, col_Tmin, col_Toptmin, col_Toptmax, col_Tmax, col_TminFactor, col_Lvpd, col_Uvpd
        del col_SFvpd_Lthres, col_SFvpd_Uthres, col_SimYield
        del colm_RUE, colm_Tmin, colm_Toptmin, colm_Toptmax, colm_Tmax, colm_TminFactor
        del colm_Lvpd, colm_Uvpd, colm_SFvpd_Lthres, colm_SFvpd_Uthres
        del colm_mae, colm_mse, colm_rmse, colm_rmsre, colm_MAPE, colm_pvalue, colm_R2  
        del colm_EF, colm_intercept, colm_slope, colm_Cb, colm_CCC, colm_accuracy
        del val_stats_SimYield, val_stats_ObsYield
        _ = gc.collect()
    
        
    def exportResults_in_batch_TPF_noStress(self, df_gy, comb, res, batch_size=4, merge=False, mergedfname='merged', 
                                          TFname="TPF_noStress_parts", removeParts=False):
        assert len(comb) == len(res)
        prt = int(len(comb)/batch_size)
        for i in tqdm(range(1, batch_size+1)):
            self.export_TPF_Results(df_gy, comb[prt*(i-1):prt*i], res[prt*(i-1):prt*i], 
                              fname="combinations_Yield_TPF_noStress_part_{}".format(i), TFname=TFname, is_VPDStress=False)
        # Merge results
        if (merge is True):
            res_path = os.path.join(self.config['RESULTS_PATH'],"TPF",TFname)
            out_path = os.path.join(self.config['RESULTS_PATH'],"TPF","TPF_noStress")
            mergeParquetFiles(res_path, out_path, fname=mergedfname, removeParts=removeParts)
            
    def exportResults_in_batch_TPF_SFvpd(self, df_gy, comb, res, batch_size=4, merge=False, mergedfname='merged', 
                                          TFname="TPF_SFvpd_parts", removeParts=False):
        assert len(comb) == len(res)
        prt = int(len(comb)/batch_size)
        for i in tqdm(range(1, batch_size+1)):
            self.export_TPF_Results(df_gy, comb[prt*(i-1):prt*i], res[prt*(i-1):prt*i], 
                              fname="combinations_Yield_TPF_SFvpd_part_{}".format(i), TFname=TFname, is_VPDStress=True)
        # Merge results
        if (merge is True):
            res_path = os.path.join(self.config['RESULTS_PATH'],"TPF",TFname)
            out_path = os.path.join(self.config['RESULTS_PATH'],"TPF","TPF_SFvpd")
            mergeParquetFiles(res_path, out_path, fname=mergedfname, removeParts=removeParts)
            
    def exportResults_in_batch_TPF(self, df_gy, comb, res, batch_size=4, is_VPDStress=False, 
                                   merge=False, mergedfname='merged', removeParts=False, calc_metrics=True):
        TFname="TPF_noStress_parts"
        fname="combinations_Yield_TPF_noStress"
        if (is_VPDStress is True): 
            TFname="TPF_SFvpd_parts"
            fname="combinations_Yield_TPF_SFvpd"
        assert len(comb) == len(res)
        prt = int(len(comb)/batch_size)
        for i in tqdm(range(1, batch_size+1)):
            self.export_TPF_Results_and_Metrics(df_gy, comb[prt*(i-1):prt*i], res[prt*(i-1):prt*i], 
                              fname="{}_part_{}".format(fname, i), TFname=TFname, is_VPDStress=is_VPDStress,
                                                calc_metrics=calc_metrics)
        # Merge results
        if (merge is True):
            res_path = os.path.join(self.config['RESULTS_PATH'],"TPF",TFname)
            out_path = os.path.join(self.config['RESULTS_PATH'],"TPF","TPF_noStress")
            if (is_VPDStress is True): 
                out_path = os.path.join(self.config['RESULTS_PATH'],"TPF","TPF_SFvpd")
            mergeParquetFiles(res_path, out_path, fname=mergedfname, removeParts=removeParts)
        # metrics allways merged
        if (calc_metrics is True):
            res_path = os.path.join(self.config['RESULTS_PATH'],"TPF","{}_metrics".format(TFname))
            out_path = os.path.join(self.config['RESULTS_PATH'],"TPF","TPF_noStress")
            if (is_VPDStress is True): 
                out_path = os.path.join(self.config['RESULTS_PATH'],"TPF","TPF_SFvpd")
            #mergeParquetFiles(res_path, out_path, fname="{}_metrics".format(mergedfname), removeParts=removeParts)
            mergeParquetFiles(res_path, out_path, fname="{}".format(fname.replace('combinations','metrics')),
                              removeParts=removeParts)
            
        
    # ---------------------------------------------
    # Merge or join parquet files
    # ---------------------------------------------
    def mergedParquetFiles(self, in_path, out_path, fname='merge', removeParts=False):
        #paths = Path(in_path).glob("*.parquet")
        #print(list(paths))
        #coalesce_parquets(paths, outpath=os.path.join(out_path, "{}.parquet".format(fname)))
        #print(pq.ParquetFile(os.path.join(out_path, "{}.parquet".format(fname) )).metadata)
        mergeParquetFiles(in_path, out_path, fname=fname, removeParts=removeParts)
            
    # ------------------------------------        
    # Old versions
    # Export directly to parquet
    ''' Export directly to parquet 
    
        Warning: Deprecated
            Stop using this function.
        
    '''
    def export_to_parquet(self, sim_res, ft='WETF', is_VPDStress=False):
        ndarray_table = None
        if (ft=='WETF'):
            TFname = "WETF_noStress"
            fname = "combinations_Yield_WETF_noStress"
            res_path = os.path.join(self.config['RESULTS_PATH'],"WETF",TFname)
            if not os.path.isdir(res_path):
                os.makedirs(res_path)

            if (is_VPDStress is False):
                ndarray_table = pa.table(
                    {
                        "UID": np.array([x[0] for x in sim_res], dtype='int32'),
                        "country": np.array([x[1] for x in sim_res], dtype='object'),
                        "location": np.array([x[2] for x in sim_res], dtype='object'),
                        "loc_code": np.array([x[3] for x in sim_res], dtype='object'),
                        "cycle": np.array([x[4] for x in sim_res], dtype='int32'),
                        "ObsYield": np.array([x[5] for x in sim_res], dtype='float32'),
                        "RUE": np.array([x[6] for x in sim_res], dtype='float32'),
                        "Tmin": np.array([x[7] for x in sim_res], dtype='float32'),
                        "Topt": np.array([x[8] for x in sim_res], dtype='float32'),
                        "Tmax": np.array([x[9] for x in sim_res], dtype='float32'),
                        "TminFactor": np.array([x[10] for x in sim_res], dtype='float32'),
                        "SimYield": np.array([x[11] for x in sim_res], dtype='float32')
                    }
                )
            else:
                TFname = "WETF_SFvpd"
                fname = "combinations_Yield_WETF_SFvpd"
                res_path = os.path.join(self.config['RESULTS_PATH'],"WETF",TFname)
                if not os.path.isdir(res_path):
                    os.makedirs(res_path)
                ndarray_table = pa.table(
                    {
                        "UID": np.array([x[0] for x in sim_res], dtype='int32'),
                        "country": np.array([x[1] for x in sim_res], dtype='object'),
                        "location": np.array([x[2] for x in sim_res], dtype='object'),
                        "loc_code": np.array([x[3] for x in sim_res], dtype='object'),
                        "cycle": np.array([x[4] for x in sim_res], dtype='int32'),
                        "ObsYield": np.array([x[5] for x in sim_res], dtype='float32'),
                        "RUE": np.array([x[6] for x in sim_res], dtype='float32'),
                        "Tmin": np.array([x[7] for x in sim_res], dtype='float32'),
                        "Topt": np.array([x[8] for x in sim_res], dtype='float32'),
                        "Tmax": np.array([x[9] for x in sim_res], dtype='float32'),
                        "TminFactor": np.array([x[10] for x in sim_res], dtype='float32'),
                        "Lvpd": np.array([x[11] for x in sim_res], dtype='float32'),
                        "Uvpd": np.array([x[12] for x in sim_res], dtype='float32'),
                        "SFvpd_Lthres": np.array([x[13] for x in sim_res], dtype='float32'),
                        "SFvpd_Uthres": np.array([x[14] for x in sim_res], dtype='float32'),
                        "SimYield": np.array([x[15] for x in sim_res], dtype='float32')
                    }
                )
            #print(ndarray_table)
        pq.write_table(ndarray_table, os.path.join(res_path,'{}.parquet'.format(fname)) )

    #          
    # Export directly to parquet
    ''' Export directly to parquet '''
    def export_to_parquet_WETF(self, sim_res, ft='WETF', is_VPDStress=False):
        
        ndarray_table = None
        if (ft=='WETF'):
            TFname = "WETF_noStress"
            fname = "combinations_Yield_WETF_noStress"
            res_path = os.path.join(self.config['RESULTS_PATH'],"WETF",TFname)
            if not os.path.isdir(res_path):
                os.makedirs(res_path)

            if (is_VPDStress is False):
                ndarray_table = pa.table(
                    {
                        "UID": [x[0] for x in sim_res],
                        "country": [x[1] for x in sim_res],
                        "location": [x[2] for x in sim_res],
                        "loc_code": [x[3] for x in sim_res],
                        "cycle": [x[4] for x in sim_res],
                        "ObsYield": [x[5] for x in sim_res],
                        "RUE": [x[6] for x in sim_res],
                        "Tmin": [x[7] for x in sim_res],
                        "Topt": [x[8] for x in sim_res],
                        "Tmax": [x[9] for x in sim_res],
                        "TminFactor": [x[10] for x in sim_res],
                        "SimYield": [x[11] for x in sim_res]
                    }
                )
            else:
                TFname = "WETF_SFvpd"
                fname = "combinations_Yield_WETF_SFvpd"
                res_path = os.path.join(self.config['RESULTS_PATH'],"WETF",TFname)
                if not os.path.isdir(res_path):
                    os.makedirs(res_path)
                ndarray_table = pa.table(
                    {
                        "UID": [x[0] for x in sim_res],
                        "country": [x[1] for x in sim_res],
                        "location": [x[2] for x in sim_res],
                        "loc_code": [x[3] for x in sim_res],
                        "cycle": [x[4] for x in sim_res],
                        "ObsYield": [x[5] for x in sim_res],
                        "RUE": [x[6] for x in sim_res],
                        "Tmin": [x[7] for x in sim_res],
                        "Topt": [x[8] for x in sim_res],
                        "Tmax": [x[9] for x in sim_res],
                        "TminFactor": [x[10] for x in sim_res],
                        "Lvpd": [x[11] for x in sim_res],
                        "Uvpd": [x[12] for x in sim_res],
                        "SFvpd_Lthres": [x[13] for x in sim_res],
                        "SFvpd_Uthres": [x[14] for x in sim_res],
                        "SimYield": [x[15] for x in sim_res]
                    }
                )
        pq.write_table(ndarray_table, os.path.join(res_path,'{}.parquet'.format(fname)) )

    #
    # -----------------
    # Metrics
    # -----------------
    ''' Create statistics or metrics for the results '''
    def export_metrics_v0(self, df_gy, batch, array_results, ft="PRFT", is_VPDStress=False):
        
        st = []
        for i in range(0, len(batch)):
            val_stats_SimYield = []
            val_stats_ObsYield = []
            for j in range(len(df_gy)):
                UID, country, location, loc_code, cycle, ObsYield = df_gy[j]
                #if (is_VPDStress is True):
                #    if (ft=="PRFT"):
                #        RUE, tminFactor, Topt, Lvpd, Uvpd, SFvpd_Lthres, SFvpd_Uthres = batch[i]
                #    elif (ft=="WETF"):
                #        RUE, Tmin, Topt, Tmax, tminFactor, Lvpd, Uvpd, SFvpd_Lthres, SFvpd_Uthres = batch[i]
                #    elif (ft=="TPF"):
                #        RUE, Tmin, Toptmin, Toptmax, Tmax, tminFactor, Lvpd, Uvpd, SFvpd_Lthres, SFvpd_Uthres = batch[i]
                #else:
                #    if (ft=="PRFT"):
                #        RUE, tminFactor, Topt = batch[i]
                #    elif (ft=="WETF"):
                #        RUE, Tmin, Topt, Tmax, tminFactor = batch[i]
                #    elif (ft=="TPF"):
                #        RUE, Tmin, Toptmin, Toptmax, Tmax, tminFactor = batch[i]

                val_stats_SimYield.append(array_results[i][j])
                val_stats_ObsYield.append(ObsYield)
                #
            # Metrics
            #print("val_stats_ObsYield ->", val_stats_ObsYield)
            #print("val_stats_SimYield ->", val_stats_SimYield)
            #m = metrics.calculateMetrics_v2(np.array(val_stats_ObsYield), np.array(val_stats_SimYield))
            m = calculateMetrics_v2(np.array(val_stats_ObsYield), np.array(val_stats_SimYield))
            st.append(list(batch[i]) + list(m))

        return st
    
    ''' Create statistics or metrics for the results '''
    def export_metrics(self, df_gy, batch, array_results):
        
        st = []
        for i in tqdm(range(0, len(batch))):
            val_stats_SimYield = []
            val_stats_ObsYield = [x for x in df_gy[:,5]] #[] #df_gy[:,5]
            for j in range(len(df_gy)):
                val_stats_SimYield.append(array_results[i][j])
            # Metrics
            m = metrics.calculateMetrics(np.array(val_stats_ObsYield), np.array(val_stats_SimYield))
            st.append(list(batch[i]) + list(m))

        return st

    # -----------------
    # Combinations
    # -----------------
    def setup_dataInput_forCombinations(self, sites):
        ''' Set up data for each sites to estimate grain yield from different combinations
    
            Parameters:
                sites (dictionary): 
                
            Return:
                Array of arrays with weather, iPAR and GPP data
        '''
        attr = []
        for _id in range(0, len(sites)):
            attr.append(sites[_id].attributes)
        
        df_GYield = pd.DataFrame(attr)[['UID', 'country', 'location', 'loc_code', 'cycle', 'ObsYield']]
        df_GYield['ObsYield'] = df_GYield['ObsYield'].astype(float).round(2)
        
        # input datasets
        arr_tn, arr_tx, arr_ndvi, arr_solrad, arr_vpdmax, arr_ipar, GPP = gpp.prepareVectors_dataset(sites)
        data_input = np.stack((arr_tn, arr_tx, arr_ndvi, arr_solrad, arr_vpdmax, arr_ipar, GPP))
        cols = arr_tn.shape[0]
        
        del attr, arr_tn, arr_tx, arr_ndvi, arr_solrad, arr_vpdmax, arr_ipar, GPP
        _ = gc.collect()
        return df_GYield, data_input, cols
    #
    #
    def getCombinations(self, functype='PRFT', cols=None, RUE=None, Tmin=None, Toptmin=None,
                        Topt=None, Toptmax=None, Tmax=None, TminFactor=None, 
                        Lvpd=None, Uvpd=None, SFvpd_Lthres=None, SFvpd_Uthres=None,
                        isVPDStress=False):
        '''
            Generate combinations to establish the optimum temperature response for grain-filling period 
            and estimate grain yield.
    
            Parameters:
                functype (str): Type of temperature response function (eg. PRFT, WETF, TPF)
                cols (int): Number of daily records in growing period
                RUE (array): Array of RUE values
                Tmin (array): Array of minimum temperature values
                Toptmin (array): Array of Optimun minimum temperature values
                Topt (array): Array of Optimun temperature values
                Toptmax (array): Array of Optimun maximum temperature values
                Tmax (array): Array of maximum temperature values
                TminFactor (array): Array of minimum temperature factor
                Lvpd (array): Array of lower VPD values
                Uvpd (array): Array of lower VPD values
                SFvpd_Lthres (array): Array of lower VPD stress factor values
                SFvpd_Uthres (array): Array of lower VPD stress factor values
                isVPDStress (bool): True/False value for using VPD stress condition
    
            Return:
                Array of combinations and an array to save results 
            
        '''
        df_cmbs = None
        array_params_to_run = None
        array_results = None
        if ((functype=='PRFT') and (isVPDStress is False)):
            # 1) Combinations for PRFT (no stress)
            if (RUE is None):
                RUE = [3.0] #[2.8, 2.9, 3.0, 3.1, 3.2]
            if (Topt is None):
                Topt = [x for x in range(15, 26)]
            if (TminFactor is None):
                TminFactor = [0.0, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5]
            combination_args_PRFT_noStress = [RUE, Topt, TminFactor]
            
            df_cmbs = pd.DataFrame(data=list(itertools.product(*combination_args_PRFT_noStress)),
                                   columns=['RUE', 'Topt', 'TminFactor'])
            df_cmbs.drop_duplicates(inplace=True)
            
            array_params_to_run = []
            for i, cmb in df_cmbs.iterrows():
                array_params_to_run.append([
                    float(cmb['RUE']),
                    float(cmb['TminFactor']),
                    float(cmb['Topt'])
                ])
            
            array_params_to_run = np.array(array_params_to_run)
    
            is_VPDStress = [isVPDStress]
            rows = array_params_to_run.shape[0]
            #cols = arr_tn.shape[0]
            array_results = np.zeros(rows*cols, dtype=np.float64).reshape(rows, cols)
    
        elif ((functype=='PRFT') and (isVPDStress is True)):
            # 2) Combinations for Topt in PRFT (VPD stress)
            if (RUE is None):
                RUE = [3.0]
            if (Topt is None):
                Topt = [x for x in range(15, 26)]
            if (TminFactor is None):
                TminFactor = [0.0, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5]
            if (Lvpd is None):
                Lvpd = [0.5, 1, 1.5, 2, 2.5, 3, 3.5]
            if (Uvpd is None):
                Uvpd = [1, 1.5, 2, 2.5, 3, 3.5, 4]
            if (SFvpd_Lthres is None):
                SFvpd_Lthres = [0.2, 0.4, 0.6, 0.8] 
            if (SFvpd_Uthres is None):
                SFvpd_Uthres = [1]
            combination_args_PRFT_vpdStress = [RUE, Topt, TminFactor, Lvpd, Uvpd, SFvpd_Lthres, SFvpd_Uthres]
            
            df_cmbs = pd.DataFrame(data=list(itertools.product(*combination_args_PRFT_vpdStress)),
                            columns=['RUE', 'Topt', 'TminFactor', 'Lvpd', 'Uvpd', 'SFvpd_Lthres', 'SFvpd_Uthres'])
            df_cmbs.drop_duplicates(inplace=True)
            # Remove combinations with Lvpd >= Uvpd
            df_cmbs = df_cmbs[~(df_cmbs['Lvpd'] >= df_cmbs['Uvpd'])]
            
            array_params_to_run = []
            for i, cmb in df_cmbs.iterrows():
                array_params_to_run.append([
                    float(cmb['RUE']),
                    float(cmb['TminFactor']),
                    float(cmb['Topt']),
                    float(cmb['Lvpd']),
                    float(cmb['Uvpd']),
                    float(cmb['SFvpd_Lthres']),
                    float(cmb['SFvpd_Uthres'])
                ])
            
            array_params_to_run = np.array(array_params_to_run)

            # Setup array for simulation results 
            is_VPDStress = [isVPDStress] #[True]
            rows = array_params_to_run.shape[0]
            #cols = arr_tn.shape[0]  
            array_results = np.zeros(rows*cols, dtype=np.float64).reshape(rows, cols)
        
        elif ((functype=='WETF') and (isVPDStress is False)):
            # 1) Combinations for WETF (no stress)
            if (RUE is None):
                RUE = [3.0]
            if (Tmin is None):
                Tmin = [x for x in range(0, 11)]
            if (Topt is None):
                Topt = [x for x in range(10, 26)]
            if (Tmax is None):
                Tmax = [x for x in range(30, 46)]
            if (TminFactor is None):
                TminFactor = [0.25]
            combination_args_WETF_noStress = [RUE, Tmin, Topt, Tmax, TminFactor]

            df_cmbs = pd.DataFrame(data=list(itertools.product(*combination_args_WETF_noStress)),
                            columns=['RUE', 'Tmin', 'Topt', 'Tmax', 'TminFactor'])
            df_cmbs.drop_duplicates(inplace=True)
            df_cmbs = df_cmbs[~(df_cmbs['Tmin'] >= df_cmbs['Topt'])]
            df_cmbs = df_cmbs[~(df_cmbs['Topt'] >= df_cmbs['Tmax'])]
            array_params_to_run = df_cmbs.values
            #array_params_to_run = []
            #for i, cmb in df_cmbs.iterrows():
            #    array_params_to_run.append([
            #        float(cmb['RUE']),
            #        float(cmb['TminFactor']),
            #        float(cmb['Tmin']),
            #        float(cmb['Topt']),
            #        float(cmb['Tmax'])
            #    ])
            #array_params_to_run = np.array(array_params_to_run)
            #is_VPDStress = [isVPDStress]
            rows = array_params_to_run.shape[0]
            #cols = arr_tn.shape[0]
            array_results = np.zeros(rows*cols, dtype=np.float64).reshape(rows, cols)
            
        elif ((functype=='WETF') and (isVPDStress is True)):
            # 2) Combinations for Topt in WETF (VPD stress)
            if (RUE is None):
                RUE = [3.0]
            if (Tmin is None):
                Tmin = [x for x in range(0, 11)]
            if (Topt is None):
                Topt = [x for x in range(10, 26)]
            if (Tmax is None):
                Tmax = [x for x in range(30, 46)]
            if (TminFactor is None):
                TminFactor = [0.25] #[0.0, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5]
            if (Lvpd is None):
                Lvpd = [0.5, 1, 1.5, 2, 2.5, 3, 3.5]
            if (Uvpd is None):
                Uvpd = [1, 1.5, 2, 2.5, 3, 3.5, 4]
            if (SFvpd_Lthres is None):
                SFvpd_Lthres = [0.2, 0.4, 0.6, 0.8] 
            if (SFvpd_Uthres is None):
                SFvpd_Uthres = [1]
            combination_args_WETF_vpdStress = [RUE,  Tmin, Topt, Tmax, TminFactor,
                                               Lvpd, Uvpd, SFvpd_Lthres, SFvpd_Uthres]

            df_cmbs = pd.DataFrame(data=list(itertools.product(*combination_args_WETF_vpdStress)),
                                   columns=['RUE', 'Tmin','Topt','Tmax', 'TminFactor', 
                                            'Lvpd', 'Uvpd', 'SFvpd_Lthres', 'SFvpd_Uthres'])
            df_cmbs.drop_duplicates(inplace=True)
            # Remove combinations with Lvpd >= Uvpd
            df_cmbs = df_cmbs[~(df_cmbs['Lvpd'] >= df_cmbs['Uvpd'])]
            df_cmbs = df_cmbs[~(df_cmbs['Tmin'] >= df_cmbs['Topt'])]
            df_cmbs = df_cmbs[~(df_cmbs['Topt'] >= df_cmbs['Tmax'])]

            array_params_to_run = df_cmbs.values
            # For testing
            #array_params_to_run = np.array([[3.0, 9.0, 18.0, 34.0, 0.25, 1.0, 4.0, 0.2, 1.0]]) 
            #RUE, Tmin, Topt, Tmax, tminFactor, Lvpd, Uvpd, SFvpd_Lthres, SFvpd_Uthres

            # For run simulations
            #is_VPDStress = [True]
            rows = array_params_to_run.shape[0]
            #cols = arr_tn.shape[0]  
            array_results = np.zeros(rows*cols, dtype=np.float64).reshape(rows, cols)
            
        elif ((functype=='TPF') and (isVPDStress is False)):
            # 1) Combinations for TPF (no stress)
            if (RUE is None):
                RUE = [3.0]
            if (Tmin is None):
                Tmin = [x for x in range(0, 5)]
            if (Toptmin is None):
                Toptmin = [x for x in range(14, 20)]
            if (Toptmax is None):
                Toptmax = [x for x in range(15, 20)]
            if (Tmax is None):
                Tmax = [x for x in range(30, 40)]
            if (TminFactor is None):
                TminFactor = [0.25]
            combination_args_TPF_noStress = [RUE, Tmin, Toptmin, Toptmax, Tmax, TminFactor]

            df_cmbs = pd.DataFrame(data=list(itertools.product(*combination_args_TPF_noStress)),
                            columns=['RUE', 'Tmin', 'Toptmin', 'Toptmax', 'Tmax', 'TminFactor'])
            df_cmbs.drop_duplicates(inplace=True)
            df_cmbs = df_cmbs[~(df_cmbs['Tmin'] >= df_cmbs['Toptmin'])]
            df_cmbs = df_cmbs[~(df_cmbs['Toptmin'] >= df_cmbs['Toptmax'])]
            df_cmbs = df_cmbs[~(df_cmbs['Toptmax'] >= df_cmbs['Tmax'])]

            array_params_to_run = df_cmbs.values
            rows = array_params_to_run.shape[0]
            array_results = np.zeros(rows*cols, dtype=np.float64).reshape(rows, cols)
            
        elif ((functype=='TPF') and (isVPDStress is True)):
            # 2) Combinations for Topt in TPF (VPD stress)
            if (RUE is None):
                RUE = [3.0]
            if (Tmin is None):
                Tmin = [x for x in range(0, 5)]
            if (Toptmin is None):
                Toptmin = [x for x in range(14, 20)]
            if (Toptmax is None):
                Toptmax = [x for x in range(15, 20)]
            if (Tmax is None):
                Tmax = [x for x in range(30, 40)]
            if (TminFactor is None):
                TminFactor = [0.25]
            if (Lvpd is None):
                Lvpd = [1.5, 2, 2.5, 3]
            if (Uvpd is None):
                Uvpd = [2, 2.5, 3, 3.5]
            if (SFvpd_Lthres is None):
                SFvpd_Lthres = [0.6, 0.8] 
            if (SFvpd_Uthres is None):
                SFvpd_Uthres = [1]

            combination_args_TPF_vpdStress = [RUE, Tmin, Toptmin, Toptmax, Tmax, TminFactor, 
                                              Lvpd, Uvpd, SFvpd_Lthres, SFvpd_Uthres]

            df_cmbs = pd.DataFrame(data=list(itertools.product(*combination_args_TPF_vpdStress)),
                                   columns=['RUE','Tmin','Toptmin','Toptmax','Tmax', 'TminFactor', 
                                            'Lvpd', 'Uvpd', 'SFvpd_Lthres', 'SFvpd_Uthres'])
            df_cmbs.drop_duplicates(inplace=True)
            # Remove combinations with Lvpd >= Uvpd
            df_cmbs = df_cmbs[~(df_cmbs['Lvpd'] >= df_cmbs['Uvpd'])]
            df_cmbs = df_cmbs[~(df_cmbs['Tmin'] >= df_cmbs['Toptmin'])]
            df_cmbs = df_cmbs[~(df_cmbs['Toptmin'] >= df_cmbs['Toptmax'])]
            df_cmbs = df_cmbs[~(df_cmbs['Toptmax'] >= df_cmbs['Tmax'])]

            array_params_to_run = df_cmbs.values
            rows = array_params_to_run.shape[0]
            array_results = np.zeros(rows*cols, dtype=np.float64).reshape(rows, cols)
    
        #
        del df_cmbs
        _ = gc.collect()
        return array_params_to_run, array_results
    #
    #
    def getGYield_forCombinations(self, functype, df_GYield, data_input, array_params_to_run, isVPDStress, 
                                  array_results, saveFile=True, returnDF=True, fmt='parquet'):
        ''' 
            Estimate grain yield using parameters for each combinations
    
            Parameters:
                df_GYield (array): Dataframe with observed grain yield
                data_input (array): Array of inputs containing weather, iPAR, PRFT and GPP data for each site.
                array_params_to_run (array): Array of combinations (RUE, TOpt and TminFactor) to simulate
                isVPDStress (bool): True/False value for using VPD stress condition
                array_results (array): Array used to save results
                saveFile (bool): True if save file in results folder
                returnDF (bool): True if return a pandas dataframe. Use False when the combinations are too large.
                fmt (str): Format to save file. Comma separate value (csv) or Apache parquet (parquet). 
                           Default is `parquet`
    
            Return:
                A dataframe or table with results
                
        '''
        is_VPDStress = [isVPDStress]
        temfun = [1] if functype=='PRFT' else [2] if functype=='WETF' else [3] if functype=='TPF' else None
        df = None
        TFname = "{}_{}".format(functype, 'SFvpd' if isVPDStress else 'noStress')
        fname = "combinations_Yield_{}_{}".format(functype, 'SFvpd' if isVPDStress else 'noStress')
        # Estimate grain yield for each combinations
        gyield.estimate(data_input, array_params_to_run, is_VPDStress, temfun, array_results)
        #if (functype=='PRFT'):
        #    #gpp.GYield(data_input, array_params_to_run, is_VPDStress, array_results)
        #    gyield.estimate(data_input, array_params_to_run, is_VPDStress, [1], array_results)
        #elif (functype=='WETF'):
        #    #gpp.GYield_WETF(data_input, array_params_to_run, is_VPDStress, array_results)
        #    gyield.estimate(data_input, array_params_to_run, is_VPDStress, [2], array_results)
        #elif (functype=='TPF'):
        #    #gpp.GYield_TPF(data_input, array_params_to_run, is_VPDStress, array_results)
        #    gyield.estimate(data_input, array_params_to_run, is_VPDStress, [3], array_results)
            
        #
        # Join results and create final dataframe
        if (returnDF is True):
            df = pd.DataFrame(gpp.createDF(df_GYield.to_numpy(), array_params_to_run, array_results))
            if (isVPDStress is False):
                if (functype=='PRFT'):
                    df.columns = ['UID', 'country', 'location', 'loc_code', 'cycle', 'ObsYield', 
                                  'RUE','TminFactor','Topt','SimYield']
                elif (functype=='WETF'):
                    df.columns = ['UID', 'country', 'location', 'loc_code', 'cycle', 'ObsYield', 
                                  'RUE','Tmin', 'Topt', 'Tmax', 'TminFactor', 'SimYield']
                elif (functype=='TPF'):
                    df.columns = ['UID', 'country', 'location', 'loc_code', 'cycle', 'ObsYield', 
                                  'RUE','Tmin', 'Toptmin', 'Toptmax', 'Tmax','TminFactor', 'SimYield']
            else:
                if (functype=='PRFT'):
                    df.columns = ['UID', 'country', 'location', 'loc_code', 'cycle', 'ObsYield', 
                                  'RUE','TminFactor','Topt','Lvpd','Uvpd','SFvpd_Lthres','SFvpd_Uthres','SimYield']
                elif (functype=='WETF'):
                    df.columns = ['UID', 'country', 'location', 'loc_code', 'cycle', 'ObsYield', 
                                  'RUE', 'Tmin', 'Topt', 'Tmax', 'TminFactor',
                                  'Lvpd','Uvpd','SFvpd_Lthres','SFvpd_Uthres','SimYield']
                elif (functype=='TPF'):
                    df.columns = ['UID', 'country', 'location', 'loc_code', 'cycle', 'ObsYield', 
                                  'RUE','Tmin', 'Toptmin', 'Toptmax', 'Tmax', 'TminFactor',
                                  'Lvpd','Uvpd','SFvpd_Lthres','SFvpd_Uthres','SimYield']
            
            if (saveFile is True):
                res_path = os.path.join(self.config['RESULTS_PATH'],functype,TFname)
                if not os.path.isdir(res_path):
                    os.makedirs(res_path)
                if (fmt=='csv'):
                    df.to_csv(os.path.join(res_path,'{}.csv'.format(fname)), index=False)
                elif (fmt=='parquet'):
                    df.to_parquet(os.path.join(res_path,'{}.parquet'.format(fname)), index=False)
                    
            #
            return df
        else:
            # No return
            if (saveFile is True):
                # Save results by default the format is Parquet due to the size of the files
                if ((functype=='WETF') and (isVPDStress is False)):
                    # Save results
                    sim_res = gpp.createDF(df_GYield.to_numpy(), array_params_to_run, array_results)
                    self.export_to_parquet_WETF(sim_res, ft='WETF', is_VPDStress=False)
                    del sim_res
                elif ((functype=='WETF') and (isVPDStress is True)):
                    self.exportResults_in_batch_WETF_SFvpd(df_GYield.to_numpy(), array_params_to_run, 
                                                       array_results, batch_size=4, 
                                                       merge=True, mergedfname='combinations_Yield_WETF_SFvpd', 
                                                       TFname="WETF_SFvpd_parts", removeParts=True)
                elif ((functype=='TPF') and (isVPDStress is False)):
                    # Save results
                    print("Saving TPF combinations and metrics to parque file...")
                    self.exportResults_in_batch_TPF(df_GYield.to_numpy(), array_params_to_run, array_results, 
                                                     batch_size=4, is_VPDStress=False, merge=True, 
                                                     mergedfname='combinations_Yield_TPF_noStress', 
                                                     removeParts=True, calc_metrics=True)

                if ((functype=='TPF') and (isVPDStress is True)):
                    self.exportResults_in_batch_TPF(df_GYield.to_numpy(), array_params_to_run, array_results, 
                                                     batch_size=16, is_VPDStress=True, merge=True, 
                                                     mergedfname='combinations_Yield_TPF_SFvpd', 
                                                     removeParts=True, calc_metrics=True)
        #
        _ = gc.collect()
    #
    #
    def getCombinations_Metrics(self, functype, isVPDStress, df_GYield, array_params_to_run, 
                                array_results, saveFile=True, fmt='parquet'):
        '''
            Get evaluation metrics for each simulation
    
            Parameters:
                functype (str): Type of temperature response function (eg. PRFT, WETF, TPF)
                isVPDStress (bool): True/False value for using VPD stress condition
                df_GYield (array): Dataframe with observed grain yield
                array_params_to_run (array): Array of combinations (RUE, TOpt and TminFactor) to simulate
                array_results (array): Array used to save results
                saveFile (bool): True if save file in results folder
                fmt (str): Format to save file. Comma separate value (csv) or Apache parquet (parquet). 
                           Default is `parquet`
            
            Return:
                A dataframe with several metrics 

            Returns:
                MAE (float): Mean Absolute Error
                MSE (float): Mean Squared Error
                RMSE (float): Root Mean Squared Error
                RMSRE (float): Root Mean Squared Relative Error
                MAPE (float): Mean Absolute Percentage Error
                pvalue (float): p-value 
                R2 (float): R Squared metric
                EF (float): Nash-Sutcliffe metric
                intercept (float):  Intercept of the regression model
                slope (float): Slope of the regression model
                Cb (float): A bias correction factor
                CCC (float): Concordance correlation coefficient
                Accuracy (float): Accuracy in percentage

        '''
        is_VPDStress = [isVPDStress]
        m = self.export_metrics(df_GYield.to_numpy(), array_params_to_run, array_results)
        TFname = "{}_{}".format(functype, 'SFvpd' if isVPDStress else 'noStress')
        fname = "metrics_Yield_{}_{}".format(functype, 'SFvpd' if isVPDStress else 'noStress')
        if (isVPDStress is True):
            if (functype=='PRFT'):
                col_names = ['RUE','TminFactor','Topt','Lvpd','Uvpd','SFvpd_Lthres','SFvpd_Uthres', 
                             "MAE","MSE","RMSE","RMSRE","MAPE","pvalue","R2","EF",
                             "intercept","slope","Cb","CCC","Accuracy" ]
            elif (functype=='WETF'):
                col_names = ['RUE','Tmin', 'Topt', 'Tmax', 'TminFactor',
                             'Lvpd','Uvpd','SFvpd_Lthres','SFvpd_Uthres', 
                             "MAE","MSE","RMSE","RMSRE","MAPE","pvalue","R2","EF",
                             "intercept","slope","Cb","CCC","Accuracy" ]
            elif (functype=='TPF'):
                col_names = ['RUE', 'Tmin', 'Toptmin', 'Toptmax', 'Tmax','TminFactor',
                             'Lvpd','Uvpd','SFvpd_Lthres','SFvpd_Uthres', 
                             "MAE","MSE","RMSE","RMSRE","MAPE","pvalue","R2","EF",
                             "intercept","slope","Cb","CCC","Accuracy" ]
        else:
            if (functype=='PRFT'):
                col_names = ["RUE", "TminFactor", "Topt", "MAE","MSE","RMSE","RMSRE","MAPE","pvalue",
                             "R2","EF", "intercept","slope","Cb","CCC","Accuracy" ]
            elif (functype=='WETF'):
                col_names = ['RUE','Tmin', 'Topt', 'Tmax', 'TminFactor',
                             "MAE","MSE","RMSE","RMSRE","MAPE","pvalue","R2","EF",
                             "intercept","slope","Cb","CCC","Accuracy" ]
            elif (functype=='TPF'):
                col_names = ['RUE','Tmin', 'Toptmin', 'Toptmax', 'Tmax','TminFactor',
                             "MAE","MSE","RMSE","RMSRE","MAPE","pvalue","R2","EF",
                             "intercept","slope","Cb","CCC","Accuracy" ]
            
        _metrics = pd.DataFrame(m, columns=col_names)
    
        if (saveFile is True):
            res_path = os.path.join(self.config['RESULTS_PATH'],functype,TFname)
            if not os.path.isdir(res_path):
                os.makedirs(res_path)
            if (fmt=='csv'):
                _metrics.to_csv(os.path.join(res_path,'{}.csv'.format(fname)), index=False)
            elif (fmt=='parquet'):
                _metrics.to_parquet(os.path.join(res_path,'{}.parquet'.format(fname)), index=False)
        
        return _metrics
    
    