# coding=utf-8
from __future__ import absolute_import, division, print_function, annotations

import warnings
warnings.simplefilter('ignore')
warnings.filterwarnings('ignore')

#from . import *

import os, gc
import shutil
import numpy as np
import pandas as pd
from datetime import date, datetime
import pandas as pd

import seaborn as sns
from matplotlib import pyplot as plt

from pathlib import Path
from typing import Callable, Iterable, TypeVar

import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.compute as pc

#from . import figures

#from sklearn import metrics  # Conflict created with metrics.py in model module
from sklearn.metrics import silhouette_score, mean_squared_error, r2_score


'''Colors class for highlight comments in CLI or notebook outputs

    Two subclasses fg for foreground and bg for background.
    Use as colors.subclass.colorname.
    i.e. colors.fg.red or colors.bg.green

    Also, the generic bold, disable, underline, reverse, strikethrough,
    and invisible work with the main class
    i.e. colors.bold

    Reset all colors with colors.reset
'''
class HT:
    reset='\033[0m'
    bold='\033[01m'
    disable='\033[02m'
    underline='\033[04m'
    reverse='\033[07m'
    strikethrough='\033[09m'
    invisible='\033[08m'
    class fg:
        black='\033[30m'
        red='\033[31m'
        green='\033[32m'
        orange='\033[33m'
        blue='\033[34m'
        purple='\033[35m'
        cyan='\033[36m'
        lightgrey='\033[37m'
        darkgrey='\033[90m'
        lightred='\033[91m'
        lightgreen='\033[92m'
        yellow='\033[93m'
        lightblue='\033[94m'
        pink='\033[95m'
        lightcyan='\033[96m'
    class bg:
        black='\033[40m'
        red='\033[41m'
        green='\033[42m'
        orange='\033[43m'
        blue='\033[44m'
        purple='\033[45m'
        cyan='\033[46m'
        lightgrey='\033[47m'


def is_number(s):
    try:
        float(s) # for int, long and float
    except ValueError:
        try:
            complex(s) # for complex
        except ValueError:
            return False
    return True

def formatInt(f):
    try:
        f = int(f)
    except:
        f = np.nan
    return f

#define function to swap columns
def swap_columns(df, col1, col2):
    col_list = list(df.columns)
    x, y = col_list.index(col1), col_list.index(col2)
    col_list[y], col_list[x] = col_list[x], col_list[y]
    df = df[col_list]
    return df

''' Get Day of the Year '''
def getDOY(d):
    day_of_year = datetime.strptime(d, "%Y-%m-%d").timetuple().tm_yday
    return day_of_year

def getDate(str_d):
    try:
        return pd.Timestamp(str_d)
    except (ValueError, TypeError):
        return np.nan    # leave unchanged
    
def getBackEmerDate(YearofEme, DOYofEme):
    try:
        ed = pd.Timestamp('{}-01-01'.format(int(YearofEme)))
        edoy = pd.DateOffset(days=int(DOYofEme))
        return ed + edoy
    except (ValueError, TypeError):
        return np.nan    # leave unchanged
    
def getHeadingDate(sowingdate, daysHead):
    try:
        return sowingdate + pd.DateOffset(days=int(daysHead))
    except (ValueError, TypeError):
        return np.nan

def getMaturityDate(sowingdate, daysMat):
    try:
        return sowingdate + pd.DateOffset(days=int(daysMat))
    except (ValueError, TypeError):
        return np.nan    
    
    
def getPhenologyDateAfterSowing(sowingdate, daysaftersowing):
    try:
        return sowingdate + pd.DateOffset(days=int(daysaftersowing))
    except (ValueError, TypeError):
        return np.nan    
    
def getObsDaysHM(matu, head):
    try:
        return int((matu - head).days)
    except (ValueError, TypeError):
        return np.nan    
    
# ---------------------------------------------
# Find nearest value or index 
# to a user define value from array
# ---------------------------------------------
def find_nearest_value(array, value):
    '''
        Find nearest value to a user define value from array
        
        Parameters:
            array (array): Array of values
            value (int): value to find into the array
        
        Returns: 
            (int): a number with the nearest value found
    '''
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def find_nearest_index(array, values):
    '''
        Find nearest index to a user define value from array
        
        Parameters:
            array (array): Array of values
            values (int): value to find into the array
        
        Returns: 
            (int): a number with the nearest index found
    '''
    values = np.atleast_1d(values)
    indices = np.abs(np.int64(np.subtract.outer(array, values))).argmin(0)
    out = array[indices]
    return out

def find_nearest(array, value):
    '''
        Find nearest index and value to a user define value from array
        
        Parameters:
            array (array): Array of values
            value (int): value to find into the array
        
        Returns: 
            (int): a number with the nearest value found
    '''
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]

'''
    Get nearest record in array

    Parameters:
        arr (array): Array of values
        value (int): Value to find into the array
        verbose (bool): Display comments if True

'''
def getNearestRow(arr=None, value=0, verbose=False):
    if (arr is None or len(arr)<=0):
        print("Please check out your input values...")
        return
    #n_value = find_nearest_value(arr, value)
    #n_index = find_nearest_index(arr, value)[0]
    n_index, n_value = find_nearest(arr, value)
    if (verbose is True):
        #print("Temperature adjusted heading to maturity days:")
        print("Nearest value: {}".format(n_value))
        print("Nearest Index: {}".format(n_index))
    return n_index, n_value

# ---------------------------
# Parquet file utils
# ---------------------------
"""coalesce_parquets.py
    gist of how to coalesce small row groups into larger row groups.
    Solves the problem described in https://issues.apache.org/jira/browse/PARQUET-1115
"""
def stream_to_parquet(path: Path, tables: Iterable[pa.Table]) -> None:
    try:
        first = next(tables)
    except StopIteration:
        return
    schema = first.schema
    with pq.ParquetWriter(path, schema) as writer:
        writer.write_table(first)
        for table in tables:
            table = table.cast(schema)  # enforce schema
            writer.write_table(table)


def stream_from_parquet(path: Path) -> Iterable[pa.Table]:
    reader = pq.ParquetFile(path)
    for batch in reader.iter_batches():
        yield pa.Table.from_batches([batch])


def stream_from_parquets(paths: Iterable[Path]) -> Iterable[pa.Table]:
    for path in paths:
        yield from stream_from_parquet(path)


"""Coalesce items into chunks. Tries to maximize chunk size and not exceed max_size.
    If an item is larger than max_size, we will always exceed max_size, so make a
    best effort and place it in its own chunk.
    You can supply a custom sizer function to determine the size of an item.
    Default is len.
    >>> list(coalesce([1, 2, 11, 4, 4, 1, 2], 10, lambda x: x))
    [[1, 2], [11], [4, 4, 1], [2]]
"""
T = TypeVar("T")
def coalesce( items: Iterable[T], max_size: int, sizer: Callable[[T], int] = len ) -> Iterable[list[T]]:
    batch = []
    current_size = 0
    for item in items:
        this_size = sizer(item)
        if current_size + this_size > max_size:
            yield batch
            batch = []
            current_size = 0
        batch.append(item)
        current_size += this_size
    if batch:
        yield batch


def coalesce_parquets(paths: Iterable[Path], outpath: Path, max_size: int = 2**20) -> None:
    tables = stream_from_parquets(paths)
    # Instead of coalescing using number of rows as your metric, you could
    # use pa.Table.nbytes or something.
    # table_groups = coalesce(tables, max_size, sizer=lambda t: t.nbytes)
    table_groups = coalesce(tables, max_size)
    coalesced_tables = (pa.concat_tables(group) for group in table_groups)
    stream_to_parquet(outpath, coalesced_tables)

def mergeParquetFiles(in_path, out_path, fname='merge', removeParts=False):
    paths = Path(in_path).glob("*.parquet")
    #print(list(paths))
    if not os.path.isdir(out_path):
        os.makedirs(out_path, exist_ok=True)
    
    coalesce_parquets(paths, outpath=os.path.join(out_path, "{}.parquet".format(fname)))
    print(pq.ParquetFile(os.path.join(out_path, "{}.parquet".format(fname) )).metadata)
    if (removeParts is True):
        try:
            shutil.rmtree(in_path)
        except OSError as e:
            print("Error: %s - %s." % (e.filename, e.strerror))

# ------------------------------
# some metrics
# ------------------------------
def CCC(y_true, y_pred):
    '''Lin's Concordance correlation coefficient

        Computes Lin's (1989, 2000) concordance correlation coefficient for 
        agreement on a continuous measure obtained by two methods. The 
        concordance correlation coefficient combines measures of both precision 
        and accuracy to determine how far the observed data deviate from the 
        line of perfect concordance (that is, the line at 45 degrees on a square 
        scatter plot). 

        Parameters:
            y_true (array): Array of observed values
            y_pred (array): Array of predicted values
        
        Returns:
            (float): Concordance correlation coefficient
    
    '''
    
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

def Cb(x,y):
    '''
        A bias correction factor that measures how far the best-fit line deviates 
        from a line at 45 degrees (a measure of accuracy). 
        
        No deviation from the 45 degree line occurs when Cb = 1. See Lin (1989 page 258).

        Parameters:
            x (array): Array of observed values
            y (array): Array of predicted values
        
        Returns:
            (float): Bias correction factor
    
    '''
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

def getScores(df, fld1=None, fld2=None):
    ''' Get stats for model results 

        Parameters:
            df (array): A pandas dataframe with Observed and Simulated values 
            fld1 (str): Name of the columns or field with observed values
            fld2 (str): Name of the columns or field with predicted values

        Returns:
            r2score (float): R squared metric
            mape (float): Mean absolute percentage error
            rmse (float): Root mean squared error
            n_rmse (float): Normalized RMSE
            d_index (float): d-index metric
            ef (float): Nash-Sutcliffe metric
            ccc (float): Concordance correlation coefficient
            cb (float): A bias correction factor
            accuracy (float): Accuracy in percentage
        
    '''
    if (df is None):
        print("Input data not valid")
        return
    if (fld1 is None or fld2 is None):
        print("Variable are not valid")
        return
    df_notnull = df[[fld1, fld2]].dropna()
    y_test = df_notnull[fld1].astype('double') #float16
    y_predicted = df_notnull[fld2].astype('double') #float16
    accuracy = round(getAccuracy(y_test, y_predicted),2)
    r2score = round(r2_score(y_test.values, y_predicted.values), 2)
    #r2 = round((np.corrcoef(y_test.values,y_predicted.values)[0, 1])**2, 2)
    rmse = round(mean_squared_error(y_test.values, y_predicted.values, squared=False),2)
    #rmse2 = np.sqrt(np.mean((y_test.values - y_predicted.values)**2))
    #n_rmse = round((rmse / y_test.values.mean()) * 100,2)
    n_rmse = round((rmse / y_test.values.mean()), 3)
    #rmsre = 100 * np.sqrt(np.mean(((y_test.values - y_test.values) / y_test.values)**2))
    mape = round(np.mean(np.abs((y_test.values - y_predicted.values)/y_test.values))*100, 2)
    d1 = ((y_test.values - y_predicted.values).astype('double') ** 2).sum()
    d2 = ((np.abs(y_predicted.values - y_test.values.mean()) + np.abs(y_test.values - y_test.values.mean())).astype('double') ** 2).sum()
    d_index = round(1 - (d1 / d2) ,3)
    # Nash–Sutcliffe model efficiency (EF)
    ef = round(1 - ( np.sum((y_test.values - y_predicted.values)**2) / np.sum((y_test.values - np.mean(y_test.values))**2) ), 2)
    # Concordance correlation coefficient
    ccc = round(CCC(y_test.values, y_predicted.values),2)
    # A bias correction factor
    cb = round(Cb(y_test.values, y_predicted.values),2)
    return r2score, mape, rmse, n_rmse, d_index, ef, ccc, cb, accuracy

'''
    Calculate accuracy

    Parameters:
        y_true (array): Array of observed values
        y_predicted (array): Array of predicted values
    
    Returns:
        (float): Accuracy in percentage
'''
def getAccuracy(y_true, y_predicted):
    mape = np.mean(np.abs((y_true - y_predicted)/y_true))*100
    if (mape<=100):
        accuracy = np.round((100 - mape), 2)
    else:
        mape = np.mean(np.abs((y_predicted - y_true)/ y_predicted))*100
        accuracy = np.round((100 - mape), 2)
    return accuracy

# ---------------------------

# ---------------------------
# Filter combinations
# ---------------------------
def filterTopCombinations(df_m, df_cmb, fnct='PRFT', VPDstress=False, top=3, Cb=1.0, ccc=0.86, rmsre=100):
    '''
        Filter combinations for selecting Top 3 of the best simulations
        
        Warning: Deprecated
            Stop using this function, instead use `filterSimulations`.
            
        Parameters:
            df_m (array): A dataframe with metrics from each simulations
            df_cmb (array): A dataframe with combinations results
            fnct (str): Temperature response function. default 'PRFT'
            VPDstress (bool): Stressed VPD. default False
            top (int): Number of selected records. default 3
            Cb (float): A threshold for Cb metric. default 1.0
            ccc (float): A threshold for CCC metric. default 0.86
            rmsre (float): A threshold for RMSRE metric. default 100

        Returns:
           (object): A dataframe with all filtered combinations
            
    '''
    # Filter the metrics
    df_filtered = df_m[(
        (df_m['RUE']==3.0) & (df_m['Cb']==Cb) & (df_m['CCC'] > ccc) & (df_m['RMSRE'] < rmsre)
    )].sort_values(['Cb','CCC','Accuracy'], ascending=False)
    #.sort_values(['Cb','Accuracy','CCC'], ascending=False) #.sort_values(['RMSRE'], ascending=True)
    
    if (fnct=='PRFT' and VPDstress is False):
        cmb_filtered = pd.DataFrame()
        for idx in df_filtered[:top].index:
            Topt = df_filtered['Topt'][idx]
            TminFactor = df_filtered['TminFactor'][idx]
            df_tmp = df_cmb[( (df_cmb['RUE']==3) & (df_cmb['TminFactor']==TminFactor) 
                             & (df_cmb['Topt']==Topt) )].reset_index(drop=True)
            cmb_filtered = pd.concat([cmb_filtered, df_tmp])
    elif (fnct=='PRFT' and VPDstress is True):
        cmb_filtered = pd.DataFrame()
        for idx in df_filtered[:top].index:
            Topt = df_filtered['Topt'][idx]
            TminFactor = df_filtered['TminFactor'][idx]
            Lvpd = df_filtered['Lvpd'][idx]
            Uvpd = df_filtered['Uvpd'][idx]
            SFvpd_Lthres = df_filtered['SFvpd_Lthres'][idx]
            SFvpd_Uthres = df_filtered['SFvpd_Uthres'][idx]

            df_tmp = df_cmb[( (df_cmb['RUE']==3) & (df_cmb['TminFactor']==TminFactor) & (df_cmb['Topt']==Topt)
                             & (df_cmb['Lvpd']==Lvpd) & (df_cmb['Uvpd']==Uvpd) & (df_cmb['SFvpd_Lthres']==SFvpd_Lthres) 
                             & (df_cmb['SFvpd_Uthres']==SFvpd_Uthres)
            )].reset_index(drop=True)
            cmb_filtered = pd.concat([cmb_filtered, df_tmp])
            
    elif (fnct=='WETF' and VPDstress is False):
        cmb_filtered = pd.DataFrame()
        for idx in df_filtered[:top].index:
            Tmin = df_filtered['Tmin'][idx]
            Topt = df_filtered['Topt'][idx]
            Tmax = df_filtered['Tmax'][idx]
            TminFactor = df_filtered['TminFactor'][idx]

            df_tmp = df_cmb[( (df_cmb['RUE']==3) & (df_cmb['TminFactor']==TminFactor) & (df_cmb['Tmin']==Tmin)
                    & (df_cmb['Topt']==Topt) & (df_cmb['Tmax']==Tmax)
            )].reset_index(drop=True)
            cmb_filtered = pd.concat([cmb_filtered, df_tmp])
    elif (fnct=='WETF' and VPDstress is True):
        cmb_filtered = pd.DataFrame()
        for idx in df_filtered[:top].index:
            Tmin = df_filtered['Tmin'][idx]
            Topt = df_filtered['Topt'][idx]
            Tmax = df_filtered['Tmax'][idx]
            TminFactor = df_filtered['TminFactor'][idx]
            Lvpd = df_filtered['Lvpd'][idx]
            Uvpd = df_filtered['Uvpd'][idx]
            SFvpd_Lthres = df_filtered['SFvpd_Lthres'][idx]
            SFvpd_Uthres = df_filtered['SFvpd_Uthres'][idx]

            df_tmp = df_cmb[( (df_cmb['RUE']==3) & (df_cmb['TminFactor']==TminFactor) & (df_cmb['Tmin']==Tmin)
                    & (df_cmb['Topt']==Topt) & (df_cmb['Tmax']==Tmax)
                    & (df_cmb['Lvpd']==Lvpd) & (df_cmb['Uvpd']==Uvpd) & (df_cmb['SFvpd_Lthres']==SFvpd_Lthres) 
                    & (df_cmb['SFvpd_Uthres']==SFvpd_Uthres)
            )].reset_index(drop=True)
            cmb_filtered = pd.concat([cmb_filtered, df_tmp])
            
    elif (fnct=='TPF' and VPDstress is False):
        cmb_filtered = pd.DataFrame()
        for idx in df_filtered[:top].index:
            Tmin = df_filtered['Tmin'][idx]
            Toptmin = df_filtered['Toptmin'][idx]
            Toptmax = df_filtered['Toptmax'][idx]
            Tmax = df_filtered['Tmax'][idx]
            TminFactor = df_filtered['TminFactor'][idx]
            
            df_tmp = df_cmb[( (df_cmb['RUE']==3) & (df_cmb['TminFactor']==TminFactor) & (df_cmb['Tmin']==Tmin)
                    & (df_cmb['Toptmin']==Toptmin) & (df_cmb['Toptmax']==Toptmax) & (df_cmb['Tmax']==Tmax)
            )].reset_index(drop=True)
            cmb_filtered = pd.concat([cmb_filtered, df_tmp])
    elif (fnct=='TPF' and VPDstress is True):
        cmb_filtered = pd.DataFrame()
        for idx in df_filtered[:top].index:
            Tmin = df_filtered['Tmin'][idx]
            Toptmin = df_filtered['Toptmin'][idx]
            Toptmax = df_filtered['Toptmax'][idx]
            Tmax = df_filtered['Tmax'][idx]
            TminFactor = df_filtered['TminFactor'][idx]
            Lvpd = df_filtered['Lvpd'][idx]
            Uvpd = df_filtered['Uvpd'][idx]
            SFvpd_Lthres = df_filtered['SFvpd_Lthres'][idx]
            SFvpd_Uthres = df_filtered['SFvpd_Uthres'][idx]

            df_tmp = df_cmb[( (df_cmb['RUE']==3) & (df_cmb['TminFactor']==TminFactor) & (df_cmb['Tmin']==Tmin)
                    & (df_cmb['Toptmin']==Toptmin) & (df_cmb['Toptmax']==Toptmax) & (df_cmb['Tmax']==Tmax)
                    & (df_cmb['Lvpd']==Lvpd) & (df_cmb['Uvpd']==Uvpd) & (df_cmb['SFvpd_Lthres']==SFvpd_Lthres) 
                    & (df_cmb['SFvpd_Uthres']==SFvpd_Uthres)
            )].reset_index(drop=True)
            cmb_filtered = pd.concat([cmb_filtered, df_tmp])
    return cmb_filtered.reset_index(drop=True)

# --------------------------------------
# Improved filter combinations function 
# --------------------------------------
def filterSimulations(functype='PRFT', VPDstress=False, cmb=None, met=None, 
                      cmb_filters=None, met_filters=None, pdFormat=False, dispFig=True, 
                      saveFig=False, figname='Fig_topComb_avgYield', figfmt='jpg',
                      saveResults=True, outputPath='./', fmt='parquet'):
    '''
        Filter a table with several combinations results from the three temperature functions.
        This function is similar to `util.filterTopCombinations` function but 
        optimized to work only with Apache parquet files.
        
        Parameters:
            functype (str): Name of the temperature response function. Default 'PRFT'
            VPDstress (bool): Stressed VPD. default False
            cmb (str): Full path of the combinations result file
            met (str): Full path of the metrics result file
            cmb_filters (dict): Dictionary with the parameters to use as a constraints in combonation filters.
                                it must match with the respective temperature function parameters. 
                                Default filters: dict( RUE = ('RUE', '=', 3.0 ), TminFactor = ('TminFactor', '=', 0.25 ) )
            met_filters (dict): Dictionary with the parameters to use as a constraints in metrics filters.
                                Default values: metfilters = dict( Cb = ('>=', 0.9), CCC = ('>=', 0.8) )
            pdFormat (bool): Export filtered dataset in pandas format. Default `False` 
            dispFig (bool): Display figure of grain yield correlation. Default is `True`
            saveFig (bool): Save figure. Default is `False`
            figname (str): Name of the figure file to be saved
            figfmt (str): Format of the figure file. `JPEG` or `PDF` formats are the available options. Default is `pdf`.
            saveResults (bool): Save filtered data. Default `True`
            outputPath (str): Output folder to save the data
            fmt (str): File format to save in csv or parquet. Default is `parquet`
            
        Returns:
           (object): A dataframe or parquet file with the filtered dataset
            
    '''
    if (functype not in ["PRFT", "WETF", "TPF"]):
        print("Not valid temperature function")
        return
    #filters = [
    #    ('RUE', '=', 3.0), ('TminFactor', '=', 0.25), 
    #    ('Tmin', '>=', 0.0), ('Tmin', '<=', 10.0), 
    #    ('Toptmin', '>=', 10.0), ('Toptmin', '<=', 25.0), 
    #    ('Toptmax', '>=', 20.0), ('Toptmax', '<=', 30.0), 
    #    ('Tmax', '>=', 30.0), ('Tmax', '<=', 45.0), 
    #    ('Lvpd', '>=', 0.5), ('Lvpd', '<=', 3.5), 
    #    ('Uvpd', '>=', 1.0), ('Uvpd', '<=', 4.0), 
    #    ('SFvpd_Lthres', '>=', 0.4), ('SFvpd_Lthres', '<=', 0.8), 
    #    ('SFvpd_Uthres', '=', 1.0)
    #]
    parameters = dict(
                RUE = ('RUE', '=', 3.0 ),
                TminFactor = ('TminFactor', '=', 0.25 )
            )
    if (cmb_filters is not None):
        parameters = {**parameters, **cmb_filters}
    
    cmb_filters = list(parameters.values())
    #
    metfilters = dict(
        Cb = ('>=', 0.9),
        CCC = ('>=', 0.8)
    )
    if (met_filters is not None):
        metfilters = {**metfilters, **met_filters}
    
    mf = None
    for v in metfilters:
        if (mf is None):
            mf = ('pc.field("{}") {} {:.2f}'.format(v, metfilters[v][0], metfilters[v][1]))
        else:
            mf = f'({mf})' +' & '+ '(pc.field("{}") {} {:.2f})'.format(v, metfilters[v][0], metfilters[v][1])
    #print(eval(mf))
    # -------------------------
    # Load combinations
    if (isinstance(cmb, str) ):
        print("Loading combinations...")
        #cmb = pd.read_parquet(cmb)
        cmb = pq.read_table( source=cmb, use_threads=True, filters=cmb_filters)
    elif (isinstance(cmb, object) and cmb is not None):
        #print("Parquet file loaded yet!")
        pass
    elif (cmb is None):
        print("Combination file not defined")
        return
    # Load metrics
    if (isinstance(met, str) ):
        print("Loading metrics...")
        #met = pd.read_parquet(met)
        met = pq.read_table( source=met, use_threads=True, filters=cmb_filters)
        print("{} combinations found".format(met.shape[0]))
    elif (isinstance(met, object) and met is not None):
        #print("Parquet file loaded yet!")
        pass
    elif (met is None):
        print("Metrics file not defined")
        return
    
    filtered_cmb = None
    #filtered_met = None
    
    # Join
    #print(cmb.shape, met.shape)
    if (VPDstress is True):
        if (functype=='PRFT'):
            fp_cm = cmb.join(met, keys=['RUE', 'TminFactor','Topt',  'Lvpd', 'Uvpd', 
                                        'SFvpd_Lthres', 'SFvpd_Uthres'])
        elif (functype=='WETF'):
            fp_cm = cmb.join(met, keys=['RUE', 'TminFactor','Tmin', 'Topt', 'Tmax',  
                                        'Lvpd', 'Uvpd', 'SFvpd_Lthres', 'SFvpd_Uthres'])
        elif (functype=='TPF'):
            fp_cm = cmb.join(met, keys=['RUE', 'TminFactor','Tmin', 'Toptmin', 'Toptmax', 
                                        'Tmax',  'Lvpd', 'Uvpd', 'SFvpd_Lthres', 'SFvpd_Uthres'])
    else:
        if (functype=='PRFT'):
            fp_cm = cmb.join(met, keys=['RUE', 'TminFactor','Topt'])
        elif (functype=='WETF'):
            fp_cm = cmb.join(met, keys=['RUE', 'TminFactor','Tmin', 'Topt', 'Tmax'])
        elif (functype=='TPF'):
            fp_cm = cmb.join(met, keys=['RUE', 'TminFactor','Tmin', 'Toptmin', 'Toptmax', 'Tmax'])
    #
    #print(fp_cm.shape)
    # Apply filter
    #filtered_cmb = fp_cm.filter((pc.field("Cb") >= 0.98) & (pc.field("CCC") >= 0.8))
    filtered_cmb = fp_cm.filter(eval(mf))
    filtered_met = met.filter(eval(mf))
    print("{} combinations found after applied filter".format(filtered_met.shape[0]))
    
    # Save file
    if (saveResults is True):
        hoy = datetime.now().strftime('%Y%m%d')
        PATH_TO_SAVE_FILTER = os.path.join(outputPath, f"{functype}_filtered")
        if not os.path.exists(PATH_TO_SAVE_FILTER):
            os.makedirs(PATH_TO_SAVE_FILTER, exist_ok=True)
        fname = "filtered_metrics_combinations_Yield_{}_{}_{}.{}".format(functype, 'SFvpd' if VPDstress else 'noStress', hoy, fmt)
        if (fmt=='csv'):
            # This could take long time if the combinations file is big
            filtered_cmb.to_pandas().to_csv(os.path.join(PATH_TO_SAVE_FILTER, fname), index=False)
        elif (fmt=='parquet'):
            pq.write_table(filtered_cmb, os.path.join(PATH_TO_SAVE_FILTER, fname))
        print("Filtered data set saved at {}".format(os.path.join(PATH_TO_SAVE_FILTER, fname)) )
    
    avgYield = None
    if (dispFig is True):# saveFig
        avgYield = filtered_cmb.group_by(["UID", "location", "loc_code"]).aggregate([("ObsYield", "mean"), ("SimYield", "mean")]).to_pandas()
        dirname = os.path.join(outputPath, f"{functype}_filtered", "Figures")
        title = '{}\n{}'.format(functype, 'VPD streess condition' if VPDstress else 'No streess condition')
        #fname = "Fig_topComb_avgYield_{}_{}_{}".format(functype, 'SFvpd' if VPDstress else 'noStress', hoy)
        figures.chart_compareResults(df_result=avgYield, fld1="ObsYield_mean", fld2="SimYield_mean",
                                     alpha=.75, s=45, xy_lim=2, hue='location', loc_leg=2, ncol=2, ha='left', va='top',
                                     title=title, #xlabel='', ylabel='', 
                                     dirname=dirname, fname=figname,  dispScore=True, dispLegend=True, 
                                     saveFig=saveFig, showFig=dispFig, fmt=figfmt)

    if (pdFormat is True):
        filtered_cmb = filtered_cmb.to_pandas(split_blocks=True, self_destruct=True)
        #del filtered_cmb # not necessary, but a good practice
    
    del avgYield
    
    _ = gc.collect()
    return filtered_cmb
    
    
# ---------------------------


