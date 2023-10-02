# coding=utf-8
from __future__ import absolute_import, division, print_function

import warnings
warnings.simplefilter('ignore')
warnings.filterwarnings('ignore')

#from . import *
import os, sys
import pandas as pd
import datetime as dt

def load_dataset():
    '''
        Load example (_Phenology, NDVI, and Weather_) datasets from International Wheat Improvement Network (IWIN) sites and 
        the other locations across the globe to estimate yield under non-stressed and Vapor pressure Deficit (VPD) 
        stressed conditions as a function of temperature.

        Return:
            A data dictionary with all raw data and information needed to carry out the demo.


        Examples:
            ``` raw
                >>> from tfunct.data import load_dataset
                >>> # Load example dataset (Phenology, NDVI and Weather data for each site)
                >>> data = load_dataset()
                >>> print(data.keys()) # dict_keys(['Pheno', 'NDVI', 'Weather'])
                >>> # Display Pheno data
                >>> data['Pheno']
                >>> # Display NDVI data
                >>> data['NDVI']
                >>> # Display Weather data
                >>> data['Weather']
            ```        
    '''
    file_path = os.path.realpath(__file__)
    example_data_path = os.path.join(file_path.replace('__init__.py',''), 'example')
    #print(example_data_path)
    #current_path = os.getcwd() #print(os.path.abspath(""))
    #if (csv is True):
    #dateparse = lambda x: dt.datetime.strptime(str(x), '%Y-%m-%d') if (str(x)!='' and str(x)!='nan') else None #'%Y-%m-%d %H:%M:%S')
    #    Weather = pd.read_csv(os.path.join(example_data_path, "Meteo_data_Final_20221219.csv"), 
    #                          index_col=False)
    #    NDVI = pd.read_csv(os.path.join(example_data_path, "NDVI_Ave_Final_20221219.csv"), 
    #                       index_col=False, parse_dates=['phenotype_date'], date_parser=dateparse)
    #    Pheno = pd.read_csv(os.path.join(example_data_path, "Pheno_Ave_Final_20221219.csv"), 
    #                        index_col=False, parse_dates=['SowingDateQC'], date_parser=dateparse)
    Pheno = pd.read_parquet(os.path.join(example_data_path, "pheno.parquet"))
    NDVI = pd.read_parquet(os.path.join(example_data_path, "ndvi.parquet"))
    Weather = pd.read_parquet(os.path.join(example_data_path, "weather.parquet"))
    
    data = {
        "Pheno":Pheno,
        "NDVI":NDVI,
        "Weather":Weather
    }
    return data