# coding=utf-8
#******************************************************************************
#
# Create figures to support paper for wheat grain yield using 3 temperature functions
# 
# version: 1.0
# Copyright: (c) Aug 15, 2023 - CIMMYT
# Authors: Azam Lashkari (azam.lashkari@jic.ac.uk)
#          Urs christoph schulthess (U.Schulthess@cgiar.org)
#          Ernesto Giron (e.giron.e@gmail.com)
#
# This source is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# This code is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# A copy of the GNU General Public License is available on the World Wide Web
# at <http://www.gnu.org/copyleft/gpl.html>. You can also obtain it by writing
# to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
# MA 02111-1307, USA.
#
#******************************************************************************

from __future__ import absolute_import, division, print_function, annotations

import warnings
warnings.simplefilter('ignore')
warnings.filterwarnings('ignore')

import os, gc
import itertools
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.gridspec as gridspec
from matplotlib.collections import PathCollection
from matplotlib.legend_handler import HandlerPathCollection, HandlerLine2D

import seaborn as sns
from matplotlib import pyplot as plt
import datetime as dt
from datetime import date, datetime
from tqdm import tqdm

#from .util import getScores
from . import *
from ..model import tday, gpp, tpf


# ------------------------------------------------
# One figure to compare Observed and Simulated Yield
# ------------------------------------------------
def chart_compareResults(df_result=None, fld1=None, fld2=None, alpha=.75, s=15, xy_lim=2, hue=None, 
                         loc_leg=2, ncol=2, ha='left', va='top',
                         title='Observed vs Simulated grain yield', xlabel='Observed Yield (tha$^{-1}$)', 
                         ylabel='Simulated Yield (tha$^{-1}$)', dirname='Figures_tFunt', fname='Fig_model_', 
                         dispScore=True, dispLegend=True, saveFig=False, showFig=True, fmt='pdf'):
    '''
        Display a scatter plot to compare two variables in the results
        
        Parameters:
            df_result (array): A pandas DataFrame with the results and variables to compare
            fld1 (str): Variable or column name to compare
            fld2 (str): Variable or column name to compare
            alpha (float): Transparency of the points in chart
            s (float): Size of the points in chart
            xy_lim (int): Used to extend the x-axis limit. Default 2 units
            hue (str): Variable to classify or discriminate the results in colors
            title (str): Title of the figure
            xlabel (str): Label of the x-axis
            ylabel (str): Label of the y-axis
            dirname (str): Folder name to save results
            fname (str): File name to save the figure
            dispScore (bool): Display the accurracy and others stats of the model
            dispLegend (bool): Display the legend of the chart
            saveFig (bool): Save file in JPG or PDF format
            fmt (str): Format of the output
        
        Returns:
            (object): A figure in JPG or PDF format with the filename specified into the folder name 
    
    '''
    if (df_result is None):
        print("Input data not valid")
        return
    if (fld1 is None or fld2 is None):
        print("Variable are not valid")
        return
    
    df = df_result.copy()
    df.dropna(subset=[fld1, fld2], inplace=True)
    
    r2score, mape, rmse, n_rmse, d_index, ef, ccc, cb, accuracy = getScores(df, fld1=fld1, fld2=fld2)
    fig, (ax1) = plt.subplots(figsize=(8,6), facecolor='white')
    fig.subplots_adjust(right=0.65)
    g1 = sns.scatterplot(x=fld1, y=fld2, data=df, alpha=alpha, color="#000000", hue=hue, s=s, lw=0, ax=ax1);
    ax1.axline((0, 0), slope=1, color='#444', ls="-", linewidth=0.75, zorder=0, label="line 1:1") #c=".5",
    maxlim = int(max(df[fld1].max(), df[fld2].max())) + xy_lim
    g1.set(xlim=(0, maxlim), ylim=(0, maxlim))
    g1.grid(visible=True, which='major', color='#d3d3d3', linewidth=0.5)
    g1.grid(visible=True, which='minor', color='#d3d3d3', linewidth=0.5)
    ax1.set_axisbelow(True)
    ax1.set_title('{}'.format(title), fontsize=15)
    ax1.set_xlabel(xlabel, fontsize=12)
    ax1.set_ylabel(ylabel, fontsize=12)
    ax1.tick_params(labelsize=12)
    if (dispScore==True):
        ax1.text(0.05,0.96,'Observations: {}\nRMSE: {:.1f}'.format(len(df), rmse) + ' tha$^{-1}$' + '\nn-RMSE: {:.3f}\nd-index: {:.2f}\nCb: {:.2f}\nCCC: {:.2f}\nR$^2$: {:.2f}\nAccuracy: {:.2f}%'.format(n_rmse, d_index, cb, ccc, r2score, accuracy), 
                 fontsize=9, ha=ha, va=va, transform=ax1.transAxes)

    if (dispLegend==True):
        plt.legend(bbox_to_anchor=(1.05, 1), loc=loc_leg, ncol=ncol, borderaxespad=0)
    else:
        #plt.legend.remove()
        ax1.get_legend().remove()
    
    # Save in PDF
    if (saveFig is True and fmt=='pdf'):
        hoy = datetime.now().strftime('%Y%m%d')
        figures_path = '{}_{}'.format(dirname, hoy)
        if not os.path.isdir(figures_path):
            os.makedirs(figures_path)
        fig.savefig(os.path.join(figures_path,"{}_{}_{}.pdf".format(fname, title.replace(' ','').replace('\n','')
                                                                   .replace('?',''), hoy)), 
                    bbox_inches='tight', orientation='portrait',  pad_inches=0.5, dpi=300)

    if (saveFig==True and (fmt=='jpg' or fmt=='png')):
        hoy = datetime.now().strftime('%Y%m%d')
        #figures_path = os.path.join(config['RESULTS_PATH'] , '{}_{}'.format(dirname, hoy) )
        figures_path = '{}_{}'.format(dirname, hoy)
        if not os.path.isdir(figures_path):
            os.makedirs(figures_path)
        fig.savefig(os.path.join(figures_path,"{}_{}_{}.{}".format(fname, title.replace(' ','').replace('\n','')
                                                                   .replace('?',''),  hoy, fmt)), 
                    bbox_inches='tight', facecolor=fig.get_facecolor(), edgecolor='none', transparent=False, dpi=300)
    
    if (showFig is True):
        if (saveFig is False):
            plt.tight_layout()
        fig.show()
    else:
        del fig
        plt.close();



# ------------------------------------------------
# Figures to compare Observed and Simulated Yield
# ------------------------------------------------
'''
    Create a correlation figure between Observed and Simulated grain yield 
    
'''
def plot_corrTempFunct(cmb_noStress=None, cmb_noStress_filtered=None, cmb_SFvpd=None, cmb_SFvpd_filtered=None,
                       functype='PRFT', fld1='ObsYield',fld2='SimYield',hue='location', 
                       ncol=6, s=80, alpha=0.45, xy_lim=1, 
                       fonts_axes=12, fonts_titles=14, dispScore=False, errorbar=False, 
                       saveFig=True, showFig=True, path_to_save_results='./', dirname='Figures', 
                       fname='Fig_1', fmt='pdf'):
    
    
    
    fig, axis = plt.subplots(1, 2, figsize=(10,5), facecolor='white', constrained_layout=True, sharex=True, sharey=True)
    #fig.suptitle('{}'.format(title), fontsize=18, y=1.05)

    def addLinearReg(df, ax, fld1, fld2):
        # Add linear regression for GY
        df_cleaned = df.dropna(subset=[fld1, fld2]).reset_index(drop=True)
        x = df_cleaned[fld1].to_numpy()
        y = df_cleaned[fld2].to_numpy()
        # determine best fit line
        par = np.polyfit(x, y, 1, full=True)
        pend=par[0][0]
        intercept=par[0][1]
        y_predicted = [pend*i + intercept  for i in x]
        l1 = sns.lineplot(x=x,y=y_predicted, color='blue', ax=ax, ls='--', lw=0.85, label='Linear Regression')
        ax.text(0.95, 0.1,r'$y$ = {:.2f}$X$ + {:.2f}'.format(pend, intercept)+'\n', 
                     fontsize=12.5, ha='right', va='top', transform=ax.transAxes)
        ax.get_legend().remove()

    def createChart(df, ax, title, f='a', errorbar=True, dispScore=False):

        ax.set_title('{}'.format(title), fontweight='bold', fontsize=16)
        ax.set_xlabel('Observed Yield (tha$^{-1}$)', fontsize=fonts_axes)
        ax.set_ylabel('Simulated Yield (tha$^{-1}$)', fontsize=fonts_axes)
        #if (f in ['a']):
        #    ax.set_xlabel('', fontsize=fonts_axes)
        if (f in ['b','e', 'f']):
            ax.set_ylabel('', fontsize=fonts_axes)
        if (f in ['c']): #'b', 
            ax.set_xlabel('', fontsize=fonts_axes)
            ax.set_ylabel('', fontsize=fonts_axes)
        ax.tick_params(labelsize=12)
        g1 = sns.scatterplot(x=fld1, y=fld2, data=df, alpha=alpha, s=s, color="#000000", marker="$\circ$", #ec="face", #markers=mks, #
                             #ci='sd', #err_style='bars',
                             style=hue, hue=hue, lw=0.6, ax=ax); #size='cycle', style="cycle",
        #g2 = sns.pointplot(x=fld1, y=fld2, data=df, errorbar="ci", capsize=.3, ax=ax) #"ci", errorbar="sd", estimator="median"
        if (errorbar is True):
            ax.errorbar(x=fld1, y=fld2, data=df, xerr="ObsYield_std", yerr='SimYield_std', linestyle='', #fmt=None,
                        elinewidth=0.5, label=None, capsize=2, capthick=0.5)

        ax.axline((0, 0), slope=1, color='#444', ls="-", linewidth=0.75, zorder=0, label="line 1:1")
        #maxlim = int(max(df[fld1].max(), df[fld2].max())) + xy_lim
        #g1.set(xlim=(0, maxlim), ylim=(0, maxlim))
        g1.set(xlim=(0, 10), ylim=(0, 10))
        g1.grid(visible=True, which='major', color='#d3d3d3', linewidth=0.5)
        g1.grid(visible=True, which='minor', color='#d3d3d3', linewidth=0.5)
        ax.set_axisbelow(True)
        # Add hztal and vertical lines
        #ax.axvline(df[fld1].mean(), ls='--', c='red', linewidth=1, label="Mean Observed Yield")
        #ax.axhline(df[fld2].mean(), ls='--', c='red', linewidth=1, label="Mean Simulated Yield")
        # Add texts
        minGY = df[fld1].min()
        maxGY = df[fld1].max()
        GY_min = df.loc[[df[fld2].argmin()]]
        GY_max = df.loc[[df[fld2].argmax()]]
        # Put text in the same place for all chart, using relative coordinates of the axes
        ax.text(0.01, 0.99, r"$\bf ("+f+")$", fontsize=14, ha='left', va='top', transform=ax.transAxes)
        if (dispScore==True):
            r2score, mape, rmse, n_rmse, d_index, ef, ccc, cb, accuracy = getScores(df, fld1=fld1, fld2=fld2)
            ax.text(.03, .90,'Observations: {}\nRMSE: {:.1f}'.format(len(df), rmse)+' tha$^{-1}$' + '\nNRMSE: {:.3f}\nCb: {:.2f}\nCCC: {:.2f}\nAccuracy: {:.1f}%\nR$^2$: {:.2f}'.format(n_rmse, cb, ccc, accuracy, r2score), 
                     fontsize=10, ha='left', va='top', transform=ax.transAxes)

        ax.get_legend().remove()


    if (functype=='PRFT'):
        # ------------------------------
        # Chart 1 - PRFT no VPD Stress
        # ------------------------------
        if (errorbar is True):
            #df = cmb_noStress.groupby(['location'], 
            #                          as_index=False).agg({'ObsYield':['mean','std'],'SimYield':['mean','std']})
            #df.columns=['location', 'ObsYield','ObsYield_std','SimYield','SimYield_std']
            df = cmb_noStress.groupby(["UID", "location", "loc_code"], 
                                      as_index=False).agg({'ObsYield':['mean','std'],'SimYield':['mean','std']})
            df.columns=['UID','location', 'loc_code' ,'ObsYield','ObsYield_std','SimYield','SimYield_std']

        else:
            df=cmb_noStress_filtered.copy()
        ax1 = axis[0]
        createChart(df, ax1, title='PRFT\nNo VPD stress', f='a', errorbar=errorbar, dispScore=dispScore)
        addLinearReg(df, ax1, fld1, fld2)
        # ------------------------------
        # Chart 2 - PRFT VPD Stress
        # ------------------------------
        if (errorbar is True):
            #df = cmb_SFvpd.groupby(['location'], 
            #                       as_index=False).agg({'ObsYield':['mean','std'],'SimYield':['mean','std']})
            #df.columns=['location', 'ObsYield','ObsYield_std','SimYield','SimYield_std']
            df = cmb_SFvpd.groupby(["UID", "location", "loc_code"],
                                   as_index=False).agg({'ObsYield':['mean','std'],'SimYield':['mean','std']})
            df.columns=['UID','location', 'loc_code' ,'ObsYield','ObsYield_std','SimYield','SimYield_std']
        else:
            df=cmb_SFvpd_filtered.copy()
        ax2 = axis[1]
        createChart(df, ax2, title='PRFT\nWith VPD stress', f='b', errorbar=errorbar, dispScore=dispScore)
        addLinearReg(df, ax2, fld1, fld2)

    elif (functype=='WETF'):
        # ------------------------------
        # Chart 1 - WETF no VPD Stress
        # ------------------------------
        if (errorbar is True):
            df = cmb_noStress.groupby(["UID", "location", "loc_code"], 
                                      as_index=False).agg({'ObsYield':['mean','std'],'SimYield':['mean','std']})
            df.columns=['UID','location', 'loc_code' ,'ObsYield','ObsYield_std','SimYield','SimYield_std']
        else:
            df=cmb_noStress_filtered.copy()
        ax1 = axis[0]
        createChart(df, ax1, title='WETF\nNo VPD stress', f='a', errorbar=errorbar, dispScore=dispScore)
        addLinearReg(df, ax1, fld1, fld2)
        # ------------------------------
        # Chart 2 - WETF VPD Stress
        # ------------------------------
        if (errorbar is True):
            df = cmb_SFvpd.groupby(["UID", "location", "loc_code"],
                                   as_index=False).agg({'ObsYield':['mean','std'],'SimYield':['mean','std']})
            df.columns=['UID','location', 'loc_code' ,'ObsYield','ObsYield_std','SimYield','SimYield_std']
        else:
            df=cmb_SFvpd_filtered.copy()
        ax2 = axis[1]
        createChart(df, ax2, title='WETF\nWith VPD stress', f='b', errorbar=errorbar, dispScore=dispScore)
        addLinearReg(df, ax2, fld1, fld2)

    elif (functype=='TPF'):
        # ------------------------------
        # Chart 1 - TPF no VPD Stress
        # ------------------------------
        if (errorbar is True):
            df = cmb_noStress.groupby(["UID", "location", "loc_code"], 
                                      as_index=False).agg({'ObsYield':['mean','std'],'SimYield':['mean','std']})
            df.columns=['UID','location', 'loc_code' ,'ObsYield','ObsYield_std','SimYield','SimYield_std']
        else:
            df=cmb_noStress_filtered.copy()
        ax1 = axis[0]
        createChart(df, ax1, title='TPF\nNo VPD stress', f='a', errorbar=errorbar, dispScore=dispScore)
        addLinearReg(df, ax1, fld1, fld2)
        # ------------------------------
        # Chart 1 - TPF VPD Stress
        # ------------------------------
        if (errorbar is True):
            df = cmb_SFvpd.groupby(["UID", "location", "loc_code"],
                                   as_index=False).agg({'ObsYield':['mean','std'],'SimYield':['mean','std']})
            df.columns=['UID','location', 'loc_code' ,'ObsYield','ObsYield_std','SimYield','SimYield_std']
        else:
            df=cmb_SFvpd_filtered.copy()
        ax2 = axis[1]
        createChart(df, ax2, title='TPF\nWith VPD stress', f='b', errorbar=errorbar, dispScore=dispScore)
        addLinearReg(df, ax2, fld1, fld2)

    # ------------------------------
    # Add legend()
    def getLegend_HandlesLabels(ax, handout, lablout):
        handles, labels = ax.get_legend_handles_labels()
        for h,l in zip(handles,labels):
            if l not in lablout:
                lablout.append(l)
                handout.append(h)
        return handout, lablout

    handout=[]
    lablout=[]
    handout, lablout = getLegend_HandlesLabels(ax1, handout, lablout)
    handout, lablout = getLegend_HandlesLabels(ax2, handout, lablout)
    
    fig.legend(handout,lablout, bbox_to_anchor=(0.5, -0.05), loc="center", ncol=ncol, 
               borderaxespad=0,fontsize=10) #, fancybox=True, shadow=True)
    
    fig.tight_layout()
    # Save in PDF
    hoy = dt.datetime.now().strftime('%Y%m%d')
    figures_path = os.path.join(path_to_save_results, '{}_{}'.format(dirname, hoy))
    if not os.path.isdir(figures_path):
        os.makedirs(figures_path)
    if (saveFig is True and fmt=='pdf'):
        fig.savefig(os.path.join(figures_path,"{}_{}.{}".format(fname, hoy, fmt)), 
                    bbox_inches='tight', orientation='portrait',  
                    edgecolor='none', transparent=False, pad_inches=0.5, dpi=300)

    if (saveFig==True and (fmt=='jpg' or fmt=='png')):
        fig.savefig(os.path.join(figures_path,"{}_{}.{}".format(fname, hoy, fmt)), 
                    bbox_inches='tight', facecolor=fig.get_facecolor(), edgecolor='none', transparent=False, dpi=300)
    
    if (showFig is True):
        fig.show()
    else:
        del fig
        plt.close();
        
# ------------------------------------------------
# Comparison between Temperature functions results
# ------------------------------------------------
def plot_TempFunct(fld1='ObsYield',fld2='SimYield',hue='location', ncol=6, s=80, alpha=0.45, xy_lim=1, 
                   fonts_axes=12, fonts_titles=14, dispScore=False, errorbar=False, 
                   saveFig=True, showFig=True, path_to_save_results='./', dirname='Figures', fname='Fig_1', fmt='pdf'):
    
    '''
        
    '''
    
    fig, axis = plt.subplots(2, 3, figsize=(15,10), facecolor='white', constrained_layout=True, sharex=True, sharey=True)
    #fig.suptitle('{}'.format(title), fontsize=18, y=1.05)

    def addLinearReg(df, ax, fld1, fld2):
        # Add linear regression for GY
        df_cleaned = df.dropna(subset=[fld1, fld2]).reset_index(drop=True)
        x = df_cleaned[fld1].to_numpy()
        y = df_cleaned[fld2].to_numpy()
        # determine best fit line
        par = np.polyfit(x, y, 1, full=True)
        pend=par[0][0]
        intercept=par[0][1]
        #print("y = {:.7f}x + {:.7f}".format(pend, intercept))
        y_predicted = [pend*i + intercept  for i in x]
        l1 = sns.lineplot(x=x,y=y_predicted, color='blue', ax=ax, ls='--', lw=0.85, label='Linear Regression')
        ax.text(0.95, 0.1,r'$y$ = {:.2f}$X$ + {:.2f}'.format(pend, intercept)+'\n', #"n: " + r"$\bf" + str(len(df)) + "$" +
                #+'\nRMSE:{:.1f}'.format(rmse)+' tha$^{-1}$' + '\nn-RMSE:{:.1f}%\nd-index:{:.2f}\nR$^2$: {:.2f}\nAccuracy: {:.1f}%'.format(n_rmse, d_index, r2score, accuracy), 
                     fontsize=12.5, ha='right', va='top', transform=ax.transAxes)
        ax.get_legend().remove()

    def createChart(df, ax, title, f='a', errorbar=True, dispScore=False):

        ax.set_title('{}'.format(title), fontweight='bold', fontsize=18)
        ax.set_xlabel('Observed Yield (tha$^{-1}$)', fontsize=fonts_axes)
        ax.set_ylabel('Simulated Yield (tha$^{-1}$)', fontsize=fonts_axes)
        if (f in ['a']):
            ax.set_xlabel('', fontsize=fonts_axes)
        if (f in ['e', 'f']):
            ax.set_ylabel('', fontsize=fonts_axes)
        if (f in ['b', 'c']):
            ax.set_xlabel('', fontsize=fonts_axes)
            ax.set_ylabel('', fontsize=fonts_axes)
        ax.tick_params(labelsize=12)
        g1 = sns.scatterplot(x=fld1, y=fld2, data=df, alpha=alpha, s=s, color="#000000", marker="$\circ$", ec="face", #markers=mks, #
                             #ci='sd', #err_style='bars',
                             style=hue, hue=hue, lw=0.6, ax=ax); #size='cycle', style="cycle",
        #g2 = sns.pointplot(x=fld1, y=fld2, data=df, errorbar="ci", capsize=.3, ax=ax) #"ci", errorbar="sd", estimator="median"
        if (errorbar is True):
            ax.errorbar(x=fld1, y=fld2, data=df, xerr="ObsYield_std", yerr='SimYield_std', linestyle='', #fmt=None,
                        elinewidth=0.5, label=None, capsize=2, capthick=0.5)
            #g2 = sns.boxplot(x=fld1, y=fld2, data=df, hue=None, ax=ax, 
            #                 flierprops={"marker": "x", "markersize":5, "markeredgecolor":"lightgray" },
            #                 boxprops={"facecolor": (.4, .6, .8, .5)}) #medianprops={"color": GScolors[phases_labels[gs]]}

        ax.axline((0, 0), slope=1, color='#444', ls="-", linewidth=0.75, zorder=0, label="line 1:1")
        #maxlim = int(max(df[fld1].max(), df[fld2].max())) + xy_lim
        #g1.set(xlim=(0, maxlim), ylim=(0, maxlim))
        g1.set(xlim=(0, 10), ylim=(0, 10))
        g1.grid(visible=True, which='major', color='#d3d3d3', linewidth=0.5)
        g1.grid(visible=True, which='minor', color='#d3d3d3', linewidth=0.5)
        ax.set_axisbelow(True)
        # Add hztal and vertical lines
        #ax.axvline(df[fld1].mean(), ls='--', c='red', linewidth=1, label="Mean Observed Yield")
        #ax.axhline(df[fld2].mean(), ls='--', c='red', linewidth=1, label="Mean Simulated Yield")
        # Add texts
        minGY = df[fld1].min()
        maxGY = df[fld1].max()
        GY_min = df.loc[[df[fld2].argmin()]]
        GY_max = df.loc[[df[fld2].argmax()]]
        # Put text in the same place for all chart, using relative coordinates of the axes
        ax.text(0.01, 0.99, r"$\bf ("+f+")$", fontsize=18, ha='left', va='top', transform=ax.transAxes)
        #ax.text(xt_tl, yt_tl-.1, "Observations: " + r"$\bf" + str(len(df))  + "$" + "\nLocations: " + r"$\bf"+str(len(df['location'].unique()))+"$",
        #         fontsize=10, ha=ha, va=va, transform=ax.transAxes)
        if (dispScore==True):
            r2score, mape, rmse, n_rmse, d_index, ef, ccc, accuracy = getScores(df, fld1=fld1, fld2=fld2)
            #print('r2score:{}, mape:{}, rmse:{}, n_rmse:{}, d_index:{}, ef:{}, ccc:{}, accuracy:{}'.format(r2score, mape, rmse, n_rmse, d_index, ef, ccc, accuracy))
            #ax.text(.03, .90,'R$^2$: {:.2f}\nRMSE: {:.1f}'.format(r2score, rmse)+' tha$^{-1}$' + '\nNRMSE: {:.3f}\nd-index: {:.2f}\nCCC: {:.2f}\nAccuracy: {:.1f}%'.format(n_rmse, d_index,  ccc, accuracy), 
            #         fontsize=12, ha='left', va='top', transform=ax.transAxes)
            #ax.text(.03, .90,'n: {:.0f}\nR$^2$: {:.2f}\nRMSE: {:.1f}'.format(len(df), r2score, rmse)+' tha$^{-1}$' + '\nNRMSE: {:.3f}\nCCC: {:.2f}\nAccuracy: {:.1f}%'.format(n_rmse, ccc, accuracy), 
            #         fontsize=12, ha='left', va='top', transform=ax.transAxes)
            ax.text(.03, .90,'R$^2$: {:.2f}\nRMSE: {:.1f}'.format(r2score, rmse)+' tha$^{-1}$' + '\nNRMSE: {:.3f}\nCCC: {:.2f}\nAccuracy: {:.1f}%'.format(n_rmse, ccc, accuracy), 
                     fontsize=12, ha='left', va='top', transform=ax.transAxes)

        ax.get_legend().remove()

    # ------------------------------
    # Chart 1 - PRFT no VPD Stress
    # ------------------------------
    if (errorbar is True):
        df = cmb_PRFT_noStress.groupby(['location'], as_index=False).agg({'ObsYield':['mean','std'],'SimYield':['mean','std']})
        df.columns=['location', 'ObsYield','ObsYield_std','SimYield','SimYield_std']
    else:
        df=cmb_PRFT_noStress_filtered.copy()
    ax1 = axis[0, 0]
    createChart(df, ax1, title='PRFT', f='a', errorbar=errorbar, dispScore=dispScore)
    addLinearReg(df, ax1)

    # ------------------------------
    # Chart 2 - WETF no VPD Stress
    # ------------------------------
    if (errorbar is True):
        df = cmb_WETF_noStress.groupby(['location'], as_index=False).agg({'ObsYield':['mean','std'],'SimYield':['mean','std']})
        df.columns=['location', 'ObsYield','ObsYield_std','SimYield','SimYield_std']
    else:
        df=cmb_WETF_noStress_filtered.copy()
    ax2 = axis[0, 1]
    createChart(df, ax2, title='WETF', f='b', errorbar=errorbar, dispScore=dispScore)
    addLinearReg(df, ax2)

    # ------------------------------
    # Chart 3 - TPF no VPD Stress
    # ------------------------------
    if (errorbar is True):
        df = cmb_TPF_noStress.groupby(['location'], as_index=False).agg({'ObsYield':['mean','std'],'SimYield':['mean','std']})
        df.columns=['location', 'ObsYield','ObsYield_std','SimYield','SimYield_std']
    else:
        df=cmb_TPF_noStress_filtered.copy()
    ax3 = axis[0, 2]
    createChart(df, ax3, title='TPF', f='c', errorbar=errorbar, dispScore=dispScore)
    addLinearReg(df, ax3)

    # ------------------------------
    # Chart 4 - PRFT VPD Stress
    # ------------------------------
    if (errorbar is True):
        df = cmb_PRFT_SFvpd.groupby(['location'], as_index=False).agg({'ObsYield':['mean','std'],'SimYield':['mean','std']})
        df.columns=['location', 'ObsYield','ObsYield_std','SimYield','SimYield_std']
    else:
        df=cmb_PRFT_SFvpd_filtered.copy()
    ax4 = axis[1, 0]
    createChart(df, ax4, title='', f='d', errorbar=errorbar, dispScore=dispScore)
    addLinearReg(df, ax4)

    # ------------------------------
    # Chart 5 - WETF VPD Stress
    # ------------------------------
    if (errorbar is True):
        df = cmb_WETF_SFvpd.groupby(['location'], as_index=False).agg({'ObsYield':['mean','std'],'SimYield':['mean','std']})
        df.columns=['location', 'ObsYield','ObsYield_std','SimYield','SimYield_std']
    else:
        df=cmb_WETF_SFvpd_filtered.copy()
    ax5 = axis[1, 1]
    createChart(df, ax5, title='', f='e', errorbar=errorbar, dispScore=dispScore)
    addLinearReg(df, ax5)

    # ------------------------------
    # Chart 6 - TPF VPD Stress
    # ------------------------------
    if (errorbar is True):
        df = cmb_TPF_SFvpd.groupby(['location'], as_index=False).agg({'ObsYield':['mean','std'],'SimYield':['mean','std']})
        df.columns=['location', 'ObsYield','ObsYield_std','SimYield','SimYield_std']
    else:
        df=cmb_TPF_SFvpd_filtered.copy()
    ax6 = axis[1, 2]
    createChart(df, ax6, title='', f='f', errorbar=errorbar, dispScore=dispScore)
    addLinearReg(df, ax6)

    # ------------------------------
    # Add legend()
    def getLegend_HandlesLabels(ax, handout, lablout):
        handles, labels = ax.get_legend_handles_labels()
        for h,l in zip(handles,labels):
            if l not in lablout:
                lablout.append(l)
                handout.append(h)
        return handout, lablout

    handout=[]
    lablout=[]
    handout, lablout = getLegend_HandlesLabels(ax1, handout, lablout)
    handout, lablout = getLegend_HandlesLabels(ax2, handout, lablout)
    handout, lablout = getLegend_HandlesLabels(ax3, handout, lablout)
    handout, lablout = getLegend_HandlesLabels(ax4, handout, lablout)
    
    fig.legend(handout,lablout, bbox_to_anchor=(0.5, -0.05), loc="center", ncol=ncol, 
               borderaxespad=0,fontsize=13) #, fancybox=True, shadow=True)
    
    fig.text(-0.02, 0.28, 'With VPD stress', ha='center', va='center', rotation='vertical', 
             fontweight='bold', fontsize=18)
    fig.text(-0.02, 0.75, 'No VPD stress', ha='center', va='center', rotation='vertical', 
             fontweight='bold', fontsize=18)
    
    fig.tight_layout()
    # Save in PDF
    hoy = dt.datetime.now().strftime('%Y%m%d')
    figures_path = os.path.join(path_to_save_results, '{}_{}'.format(dirname, hoy))
    if not os.path.isdir(figures_path):
        os.makedirs(figures_path)
    if (saveFig is True and fmt=='pdf'):
        fig.savefig(os.path.join(figures_path,"{}_{}.{}".format(fname, hoy, fmt)), 
                    bbox_inches='tight', orientation='portrait',  
                    edgecolor='none', transparent=False, pad_inches=0.5, dpi=300)

    if (saveFig==True and (fmt=='jpg' or fmt=='png')):
        fig.savefig(os.path.join(figures_path,"{}_{}.{}".format(fname, hoy, fmt)), 
                    bbox_inches='tight', facecolor=fig.get_facecolor(), edgecolor='none', transparent=False, dpi=300)
    
    if (showFig is True):
        fig.show()
    else:
        del fig
        plt.close();
#

# ======================================
# Figures to compare several simulations
# --------------------------------------

'''
    Create values for Tday used in detailed figures
    Estimate TDay and Temperature response type I function (PRFT) values from heading to maturity
'''
def createFigure_Type_I_PRFT(sites, cmb=None, roundVal=3, maxTDay=50, saveTable=False, 
                             path_to_save_results='./', fmt='parquet', verbose=False):
    #
    # Setup table with all values for displaying all locations in one figure 
    df = cmb.copy() #[cmb_PRFT_noStress['UID']==uid].reset_index(drop=True)
    nfig = 0
    df_TDays_vs_TempResponse = pd.DataFrame()
    
    def generating_TDayValues_PRFT(site=None, TminFactor=0.25, TOpt=18, isVPDstress=False, 
                                   roundVal=1, maxTDay=40, verbose=False ):
        if (TminFactor is None or TOpt is None):
            print("Parameters not found. Please check the inputs")
            return

        if (verbose):
            print("Getting TDay values for site")
        # Get TDays values
        tDay_site0 = tday.estimate_TDay(site.inputWPN['TMIN'].to_numpy(), 
                                site.inputWPN['TMAX'].to_numpy(), TminFactor )
        #print("tDay_site0", tDay_site0)
        tDay_site = [x/10 for x in range(0, maxTDay*10)] # Create curve values
        tDay_site = sorted(np.unique(np.concatenate([tDay_site, tDay_site0]))) #.round(roundVal)
        prft_site = 1 - 0.0025*(np.asarray(tDay_site)-TOpt)**2 # Get PRFT
        arrPRFT = []
        for x in prft_site: # Clean values lower than zero
            if x > 0.0:
                arrPRFT.append(x)
            else:
                arrPRFT.append(0.0)

        return np.array(tDay_site).round(roundVal), np.array(arrPRFT).round(roundVal)
    
    # Get average values
    def get_mean_TDays_vs_TempResponse(df, TminFactor=0.25, Topt=15, verbose=False):
        df_f = pd.DataFrame()
        for uid in df['UID'].unique():
            df_1 = df[( (df['UID']==uid) & (df['TminFactor']==TminFactor) & (df['Topt']==Topt)
            )][['UID', 'TDay', 'TempResponse']].reset_index(drop=True)
            df_f = pd.concat([df_f, df_1], axis=1)

        TDay_avg = df_f['TDay'].mean(axis=1) #.values
        TempResponse_avg = df_f['TempResponse'].mean(axis=1) #.values
        df_avg = pd.concat([TDay_avg, TempResponse_avg], axis=1)
        df_avg.columns = ['TDay', 'TempResponse']
        df_avg['TminFactor'] = TminFactor
        df_avg['Topt'] = Topt
        df_f = None
        del df_f
        _ = gc.collect()
        return df_avg
    # ------- 
    # Process each site
    if (verbose):
        print("Processing all sites...")
    for uid in tqdm(df['UID'].unique()):
        df2 = df[df['UID']==uid].reset_index(drop=True)
        country = df2['country'].unique()[0]
        loc = df2['location'].unique()[0]
        loc_code = df2['loc_code'].unique()[0]
        cycle = df2['cycle'].unique()[0]
        rue = df2['RUE'].unique()[0]
        obsYield = df2['ObsYield'].mean()
        simYield = df2['SimYield'].mean()
        #print("processing {} - {}".format(uid, loc))
        #
        sel_cmb_prft = df2.groupby(['UID', 'TminFactor', 'Topt'], as_index=False).agg({'ObsYield':'mean'})
        for idx in sel_cmb_prft.index:
            #_uid=sel_cmb_prft['UID'][idx]
            TminFactor=sel_cmb_prft['TminFactor'][idx]
            TOpt=sel_cmb_prft['Topt'][idx]
            site=sites[uid-1]
            # Display charts for each combinations in one figure
            xVals, yVals = generating_TDayValues_PRFT(site=site, TminFactor=TminFactor, TOpt=TOpt, 
                                                      roundVal=roundVal, maxTDay=maxTDay, verbose=verbose)
            #print(xVals, yVals)
            # Create new DF
            #df2['TDay'] = xVals
            #df2['TempResponse'] = yVals
            df3 = pd.DataFrame({'TDay': xVals, 'TempResponse': yVals})
            df3['UID'] = uid
            df3['country'] = country
            df3['location'] = loc
            df3['loc_code'] = loc_code
            df3['cycle'] = cycle
            df3['RUE'] = rue
            df3['TminFactor'] = TminFactor
            df3['Topt'] = TOpt
            df3['ObsYield'] = obsYield
            df3['SimYield'] = simYield
            #df.loc[(df['UID']==uid),'TDay'] = xVals
            #df.loc[(df['UID']==uid),'TempResponse'] = yVals
            df_TDays_vs_TempResponse = pd.concat([df_TDays_vs_TempResponse, df3]) #, ignore_index=True)
    #
    del df3
    df_TDays_vs_TempResponse['TDay'] = df_TDays_vs_TempResponse['TDay'].astype(float).round(2)
    df_TDays_vs_TempResponse['TempResponse'] = df_TDays_vs_TempResponse['TempResponse'].astype(float).round(3)
    
    # ---------------------------
    #if (verbose):
    print("Calculating average value for all simulations...")
    df_TDays_vs_TempResponse_mean_allSites = pd.DataFrame()
    sel_cmb_prft = df_TDays_vs_TempResponse.groupby(['TminFactor', 'Topt'], as_index=False).agg('first')
    for idx in tqdm(sel_cmb_prft.index):
        TminFactor=sel_cmb_prft['TminFactor'][idx]
        Topt=sel_cmb_prft['Topt'][idx]
        df_avg_tmfctor_topt = get_mean_TDays_vs_TempResponse(df=df_TDays_vs_TempResponse, 
                                                             TminFactor=TminFactor, Topt=Topt)
        df_TDays_vs_TempResponse_mean_allSites = pd.concat([df_TDays_vs_TempResponse_mean_allSites, df_avg_tmfctor_topt], axis=0)

    if (saveTable):
        hoy = dt.datetime.now().strftime('%Y%m%d')
        #path_to_save_results = os.path.join(path_to_save_results, '{}_{}'.format('Data_for_Figures', hoy))
        if not os.path.isdir(path_to_save_results):
            os.makedirs(path_to_save_results)
        if (fmt=='csv'):
            df_TDays_vs_TempResponse.to_csv(os.path.join(path_to_save_results, f'PRFT_dataForCharting_TDays_vs_TempResponse_{hoy}.csv'), index=False)
            df_TDays_vs_TempResponse_mean_allSites.to_csv(os.path.join(path_to_save_results , f'PRFT_dataForCharting_TDays_vs_TempResponse_mean_allSites_{hoy}.csv'), index=False)
        elif (fmt=='parquet'):
            df_TDays_vs_TempResponse.to_parquet(os.path.join(path_to_save_results, f'PRFT_dataForCharting_TDays_vs_TempResponse_{hoy}.parquet'), index=False)
            df_TDays_vs_TempResponse_mean_allSites.to_parquet(os.path.join(path_to_save_results , f'PRFT_dataForCharting_TDays_vs_TempResponse_mean_allSites_{hoy}.parquet'), index=False)
    
    return df_TDays_vs_TempResponse, df_TDays_vs_TempResponse_mean_allSites


def diplay_Figure_Type_I_PRFT_InOneFig(cmb=None, fnct='PRFT', df_tdays=None, df_tdays_mean=None,
                                       saveFig=True, showFig=True, fmt='jpg', leg_ncol=3,
                                       path_to_save_results='./', fname=None):

    def createFigTDay_vs_TempResponseBySiteYrs_InOneFig(df_dailyVals=None, uid=1, fnct='PRFT', ax=None):
        '''
            Temperature response type I function for heading to maturity
            Display TDay (°C) vs Relative temperature response for all sites
        '''
        # Get combinations parameters
        df = df_dailyVals[df_dailyVals['UID']==uid].reset_index(drop=True).sort_values(['TDay'])
        #df = df.drop_duplicates(subset=['TDay', 'TempResponse']) # this produce strange plots
        # Unique combinations by site
        sel_cmb = df.groupby(['UID', 'TminFactor', 'Topt'], as_index=False).agg({'first'})
        # Create figure
        #print("Processing curves for each site...")
        for idx in sel_cmb.index:
            uid=sel_cmb['UID'][idx]
            TminFactor=sel_cmb['TminFactor'][idx]
            TOpt=sel_cmb['Topt'][idx]
            df_t = df[((df['TminFactor']==TminFactor) & (df['Topt']==TOpt))]
            x=df_t['TDay'].to_numpy()
            y=df_t['TempResponse'].to_numpy()
            lbl = None if (len(df['UID'].unique())==1) else 'Combinations for all sites'
            ax.plot(np.arange(len(y)), y, marker='o', ls='--', linewidth=0.25, markersize=0.01, 
                    color='gainsboro', zorder=0, label=lbl)
                    #label='TminFactor: {:.2f} - TOpt: {:.0f}°C'.format(TminFactor, TOpt))
            if (idx == (len(sel_cmb)-1)):
                # Just to get the corrected x-axis labels 
                sns.pointplot(x=x, y=y,  color=None, linestyles='', ax=ax, markers='')

    #
    def createFigAverageBySiteYrs_InOneFig(df_dailyVals=None, fnct='PRFT', ax=None):
        '''
            Temperature response type I function for heading to maturity
            Display Average TDay (°C) vs Average Relative temperature response for all sites
        '''
        # Get combinations parameters
        df = df_dailyVals.copy()
        # Unique combinations by site
        sel_cmb = df.groupby(['TminFactor', 'Topt'], as_index=False).agg({'first'})
        # Create figure
        palette = sns.color_palette("Set1", len(sel_cmb))
        line = itertools.cycle(['-', '--', '-.', ':', 'solid', 'dashed', 'dashdot', 'dotted'])
        lines = [next(line) for i in range(len(sel_cmb))]
        marker = itertools.cycle(['o', '^', '*', '8', 's', 'p', 'd', 'v'])
        markers = [next(marker) for i in range(len(sel_cmb))]
        #print("Processing average curves for all site...")
        for idx in tqdm(sel_cmb.index):
            TminFactor=sel_cmb['TminFactor'][idx]
            TOpt=sel_cmb['Topt'][idx]
            df_t = df[((df['TminFactor']==TminFactor) & (df['Topt']==TOpt))].reset_index(drop=True).sort_values(['TDay'])
            x=df_t['TDay'].to_numpy().round(3)
            y=df_t['TempResponse'].to_numpy()
            ax.plot(np.arange(len(y)), y, marker=markers[idx], ls=lines[idx], linewidth=0.25, markersize=0.01, 
                    color=palette[idx], zorder=0, 
                    label='TminFactor: {:.2f} - TOpt: {:.0f}°C'.format(TminFactor, TOpt))
            if (idx == (len(sel_cmb)-1)):
                # Just to get the corrected x-axis labels 
                sns.pointplot(x=x, y=y,  color=None, linestyles='', ax=ax, markers='')
        # 
    # -----------------------

    # Create Figure
    df = cmb.copy()
    
    # ----------
    fig, ax = plt.subplots(1, 1, figsize=(8,6), facecolor='white')
    handout=[]
    lablout=[]

    # All combinations
    # Average combinations
    if (len(df['UID'].unique())>1):
        for nfig in tqdm(range(len(df['UID'].unique()))):
            df2 = df[df['UID']==nfig+1].reset_index(drop=True)
            createFigTDay_vs_TempResponseBySiteYrs_InOneFig(df_tdays, uid=nfig+1, fnct=fnct, ax=ax)
    #
    # Average combinations
    createFigAverageBySiteYrs_InOneFig(df_tdays_mean, fnct=fnct, ax=ax)
    #
    plt.setp(ax.get_lines(),linewidth=0.75)  # set lw for all lines of g axes
    plt.setp(ax.collections, sizes=[0.01])
    
    # Agregar Optimun temperature guide lines 
    v_lines = []
    v_lines_values = ','.join([str(x) for x in df['Topt'].unique()]).split(',')
    #print(v_lines_values)
    for label in ax.get_xticklabels(): #minor=True
        for v in v_lines_values: #['15.0', '16.0', '17.0']:
            if ( round(float(v),1) == round(float(label.get_text()),1) and '{:,.0f}'.format(float(v)) not in v_lines):
                if (len(v_lines) < len(v_lines_values)):
                    v_lines.append(label.get_position()[0])
                    
    #print(v_lines[:len(v_lines_values)])
    ax.vlines(v_lines, 0, 1, ls='-.', color='lightgray', linewidth=0.75) #, label='Optimum Temperature (15, 16, 17°C)')

    ax.set_xticks(ax.get_xticks()[::50]) # Number of TDays
    ax.tick_params(axis='x', labelsize=10, color='lightgray', rotation=90)
    ax.tick_params(axis='y', labelsize=10, color='lightgray')
    ax.grid(color='gainsboro', linestyle='-', linewidth=0.5, zorder=-1)
    ax.set_axisbelow(True)
    
    # Temperature response of type I (one cardinal temperature) function for heading to maturity across all sites.
    # Validate if only one site
    fname = "Figure_PRFT_SiteYrs_Comparison" if fname is None else fname
    if (len(df['UID'].unique())==1):
        country = str(df.loc[0, 'country'])
        loc = str(df.loc[0, 'location'])
        cycle = str(df.loc[0, 'cycle'])
        title1='Temperature response function Type I ({})\n{} - {} - {}'.format(fnct, country, loc, cycle)
        fname = "Figure_PRFT_Site_{}_{}_{}_Comparison".format(loc, country, cycle)
    else:
        title1='Average temperature response function Type I ({})'.format(fnct)
    # Title for the complete figure
    fig.suptitle('{}'.format(title1), fontsize='x-large', fontweight='bold' ) #, fontsize=20, y=1.005)
    ax.set_xlabel('TDay (°C)', fontweight='bold', fontsize=12)
    ax.set_ylabel('Relative response', fontweight='bold', fontsize=12)
    
    def getLegend_HandlesLabels(ax, handout, lablout):
        handles, labels = ax.get_legend_handles_labels()
        for h,l in zip(handles,labels):
            if l not in lablout:
                lablout.append(l)
                handout.append(h)
        return handout, lablout
    
    handout, lablout = getLegend_HandlesLabels(ax, handout, lablout)
    fig.legend(handout,lablout, bbox_to_anchor=(0.5, -0.25), loc="lower center", ncol=leg_ncol, 
               numpoints=1, markerscale=2,
               borderaxespad=0,bbox_transform=plt.gcf().transFigure, fontsize=10 ) #, fancybox=True, shadow=True)

    fig.tight_layout()
    hoy = dt.datetime.now().strftime('%Y%m%d')
    #figures_path = os.path.join(path_to_save_results, 'Figures') #fnct) #'{}_{}'.format(dirname, hoy))
    figures_path = os.path.join(path_to_save_results, '{}_{}'.format('Figures', hoy))
    if not os.path.isdir(figures_path):
        os.makedirs(figures_path)
    if (saveFig==True and (fmt=='jpg' or fmt=='png')):
        fig.savefig(os.path.join(figures_path, "{}_{}.{}".format(fname, hoy,fmt)), 
                    bbox_inches='tight', transparent=False, dpi=300)
    if (saveFig is True and fmt=='pdf'):
        fig.savefig(os.path.join(figures_path, "{}_{}.pdf".format(fname, hoy)), 
                    bbox_inches='tight', orientation='portrait', edgecolor='none', transparent=False, 
                    pad_inches=0.5, dpi=300)

    if (showFig is True):
        fig.show()
    else:
        del fig
        plt.close();
#

def display_FigTDay_vs_TempResponseBySiteYrs(cmb=None, df_tdays=None, roundVal=2, maxTDay=50, fnct='PRFT', 
                                             saveFig=True, showFig=True, fmt='jpg', cols=5, leg_ncol=4,
                                             path_to_save_results='./'):
    df = cmb.copy()
                                       
    def createFigTDay_vs_TempResponseBySiteYrs(df_cmb=None, df_dailyVals=None, uid=1, fnct='PRFT', ax=None):
        '''
            Temperature response type I function for heading to maturity
            Display TDay (°C) vs Relative temperature response by site
        '''

        # Get site data
        df0 = df_cmb[df_cmb['UID']==uid].reset_index(drop=True)
        country = df0['country'].unique()[0]
        loc = df0['location'].unique()[0]
        loc_code = df0['loc_code'].unique()[0]
        cycle = df0['cycle'].unique()[0]
        obsYield = df0['ObsYield'].mean()
        simYield = df0['SimYield'].mean()

        # Get combinations parameters
        df = df_dailyVals[df_dailyVals['UID']==uid].reset_index(drop=True).sort_values(['TDay'])
        #df = df.drop_duplicates(subset=['TDay', 'TempResponse']) # this produce strange plots
        # Unique combinations by site
        sel_cmb = df.groupby(['UID', 'TminFactor', 'Topt'], as_index=False).agg({'ObsYield':'mean'})
        # Create figure
        #fig, ax = plt.subplots(figsize=(8,6))
        palette = sns.color_palette("Set1", len(sel_cmb))
        line = itertools.cycle(['-', '--', '-.', ':', 'solid', 'dashed', 'dashdot', 'dotted'])
        lines = [next(line) for i in range(len(sel_cmb))]
        marker = itertools.cycle(['o', '^', '*', '8', 's', 'p', 'd', 'v'])
        markers = [next(marker) for i in range(len(sel_cmb))]
        for idx in sel_cmb.index:
            uid=sel_cmb['UID'][idx]
            TminFactor=sel_cmb['TminFactor'][idx]
            TOpt=sel_cmb['Topt'][idx]
            df_t = df[((df['TminFactor']==TminFactor) & (df['Topt']==TOpt))]
            x=df_t['TDay'].to_numpy()
            y=df_t['TempResponse'].to_numpy()
            ax.plot(np.arange(len(y)), y, marker=markers[idx], ls=lines[idx], linewidth=0.25, markersize=0.01, 
                    color=palette[idx], zorder=0, 
                    label='TminFactor: {:.2f} - TOpt: {:.0f}°C'.format(TminFactor, TOpt))
            if (idx == (len(sel_cmb)-1)):
                # Just to get the corrected x-axis labels 
                sns.pointplot(x=x, y=y,  color=None, linestyles='', ax=ax, markers='')
        # 
        plt.setp(ax.get_lines(),linewidth=0.75)  # set lw for all lines of g axes
        plt.setp(ax.collections, sizes=[0.01])
        v_lines = []
        v_lines_values = ','.join([str(x) for x in df0['Topt'].unique()]).split(',')
        for label in ax.get_xticklabels(): #minor=True
            for v in v_lines_values: #['15.0', '16.0', '17.0']:
                if ( round(float(v),1) == round(float(label.get_text()),1) and '{:,.0f}'.format(float(v)) not in v_lines):
                    if (len(v_lines) < len(v_lines_values)):
                        v_lines.append(label.get_position()[0])

        ax.vlines(v_lines, 0, 1, ls='-.', color='lightgray', linewidth=0.75) #, label='Optimum Temperature (15, 16, 17°C)')
        # ------------
        ax.set_xticks(ax.get_xticks()[::50])
        ax.tick_params(axis='x', labelsize=6, color='lightgray', rotation=90)
        ax.tick_params(axis='y', labelsize=6, color='lightgray') #, rotation=90)

        ax.set_title('{} - {} - {} - {}'.format(uid, loc, loc_code, cycle)) #, fontweight='bold', fontsize=16)
        #ax.set_xlabel('TDay (°C)', fontweight='bold', fontsize=14)
        #ax.set_ylabel('Relative temperature response (°C)', fontweight='bold', fontsize=14)
        ax.set_xlabel('')  
        ax.set_ylabel('')
        #ax.spines['left'].set_bounds(0, 1)
        #ax.spines['bottom'].set_bounds(0, 40)
        #ax.spines['right'].set_visible(False)
        #ax.spines['top'].set_visible(False)
        ax.grid(color='gainsboro', linestyle='-', linewidth=0.5, zorder=-1)
        ax.set_axisbelow(True)


    def getLegend_HandlesLabels(ax, handout, lablout):
        handles, labels = ax.get_legend_handles_labels()
        for h,l in zip(handles,labels):
            if l not in lablout:
                lablout.append(l)
                handout.append(h)
        return handout, lablout
    
    
    # ------------------------------
    # Create Figure
    fig = plt.figure(figsize=(14, 20), facecolor='white', tight_layout=True) #, constrained_layout=True) 
    #cols = 5
    rows = int(np.ceil(len(df['UID'].unique()) / cols))
    gs = gridspec.GridSpec(nrows=rows, ncols=cols) #, wspace=0, hspace=0)

    handout=[]
    lablout=[]
    nfig = 0
    for nr in tqdm(range(rows)):
        for nc in range(cols):
            #if (nfig>0):
            #    ax = fig.add_subplot(gs[nr, nc], sharey=ax, sharex=ax)
            #else:
            ax = fig.add_subplot(gs[nr, nc])

            #if (nfig<(rows*cols)-(cols+1)):
            #    plt.setp(ax.get_xticklabels(), visible=False)
            #    #plt.setp(ax.get_yticklabels(), visible=False)
            # Create figure
            #print("Processing site-year: ", nfig+1)
            createFigTDay_vs_TempResponseBySiteYrs(df, df_tdays, uid=nfig+1, fnct=fnct, ax=ax)
            # Legend
            handout, lablout = getLegend_HandlesLabels(ax, handout, lablout)
            nfig = nfig + 1

    # Temperature response of type I (one cardinal temperature) function for heading to maturity across all sites.
    fig.text(0.5, -0.01, 'TDay (°C)', ha='center', va='center', fontweight='bold', fontsize=14)
    fig.text(-0.02, 0.5, 'Relative temperature response', ha='center', va='center', 
                     rotation='vertical', fontweight='bold', fontsize=14)
    #title1='Temperature response PRFT function for 50 site-years around the globe' #.format()
    # Title for the complete figure
    #fig.suptitle('{}'.format(title1), fontsize='x-large', fontweight='bold' ) #, fontsize=20, y=1.005)

    def updatescatter(handle, orig):
        handle.update_from(orig)
        handle.set_sizes([4])

    def updateline(handle, orig):
        handle.update_from(orig)
        handle.set_markersize(1)

    fig.legend(handout,lablout, bbox_to_anchor=(0.5, -0.08), loc="lower center", ncol=leg_ncol, 
               numpoints=1, markerscale=2,
               borderaxespad=0,bbox_transform=plt.gcf().transFigure, fontsize=12, 
               handler_map={PathCollection : HandlerPathCollection(update_func=updatescatter),
                            plt.Line2D : HandlerLine2D(update_func = updateline)}
              ) #, fancybox=True, shadow=True)

    fig.tight_layout()
    hoy = dt.datetime.now().strftime('%Y%m%d')
    #figures_path = os.path.join(path_to_save_results, fnct) #'{}_{}'.format(dirname, hoy))
    figures_path = os.path.join(path_to_save_results, '{}_{}'.format('Figures', hoy))
    if not os.path.isdir(figures_path):
        os.makedirs(figures_path)
    if (saveFig==True and (fmt=='jpg' or fmt=='png')):
        fig.savefig(os.path.join(figures_path, "Figure_PRFT_BySiteYrs_{}.{}".format(hoy,fmt)), bbox_inches='tight',
                    transparent=False, dpi=300)
    if (saveFig is True and fmt=='pdf'):
        fig.savefig(os.path.join(figures_path, "Figure_PRFT_BySiteYrs_{}.pdf".format(hoy)), 
                    bbox_inches='tight', orientation='portrait', edgecolor='none', transparent=False, 
                    pad_inches=0.5, dpi=300)

    if (showFig is True):
        fig.show()
    else:
        del fig
        plt.close();

# --------------------------------------
# WETF
# --------------------------------------

'''
    Create values for Tday used in detailed figures
    Estimate TDay and Temperature response type II function (WETF) values from heading to maturity
'''
def createFigure_Type_II_WETF(sites, cmb=None, isVPDStress=False, roundVal=3, maxTDay=50, saveTable=False, 
                             path_to_save_results='./', fmt='parquet', verbose=False):
    #
    # Setup table with all values for displaying all locations in one figure 
    df = cmb.copy()
    nfig = 0
    df_WETF_TDays_vs_TempResponse = pd.DataFrame()
    
    def generating_TDayValues_WETF(site=None, RUE=3.0, TminFactor=0.25, Tmin=6, Topt=18, Tmax=24, 
                                   isVPDStress=False, roundVal=1, maxTDay=40, verbose=False ):
        if (TminFactor is None or TOpt is None):
            print("Parameters not found. Please check the inputs")
            return

        if (verbose):
            print("Getting TDay values for site")
        # Get TDays values
        tDay_site0 = tday.estimate_TDay(site.inputWPN['TMIN'].to_numpy(), 
                                site.inputWPN['TMAX'].to_numpy(), TminFactor )
        #print("tDay_site0", tDay_site0)
        tDay_site = [x/10 for x in range(0, maxTDay*10)] # Create curve values
        tDay_site = sorted(np.unique(np.concatenate([tDay_site, tDay_site0]))) #.round(roundVal)
        
        # Get WETF
        WETFTMAX = gpp.apply_WETF(np.asarray(tDay_site), Tmin, Topt, Tmax)
        
        #tn, tx, ndvi, solrad, VPDx, ipar, GPP = gpp.prepareVectors_dataset([site])
        #df_GYield, data_input, cols = model.setup_dataInput_forCombinations(sites) # Setup input data
        #array_params_to_run, array_results = model.getCombinations(functype='WETF', cols=cols, RUE=[RUE], 
        #                                                           Tmin=[Tmin], Topt=[Topt], Tmax=[Tmax],
        #                                                           TminFactor=[TminFactor], 
        #                                                           isVPDStress=isVPDStress)
        #
        #GPP = gpp.estimate( data_input, array_params_to_run, [isVPDStress], [2], array_results)
        #print(GPP)
        arrWETF = []
        for x in WETFTMAX: # Clean values lower than zero
            if x > 0.0:
                arrWETF.append(x)
            else:
                arrWETF.append(0.0)

        return np.array(tDay_site).round(roundVal), np.array(arrWETF).round(roundVal)
    
    #
    def get_mean_TDays_vs_TempResponse_WETF(df, isVPDStress=False, TminFactor=0.25, Tmin=9.0, Topt=15, Tmax=32):
        df_f = pd.DataFrame()
        for uid in df['UID'].unique():
            df_1 = df[( (df['UID']==uid)  & (df['TminFactor']==TminFactor) 
                       & (df['Tmin']==Tmin) & (df['Topt']==Topt) & (df['Tmax']==Tmax)
            )][['UID', 'TDay', 'TempResponse']].reset_index(drop=True)
            df_f = pd.concat([df_f, df_1], axis=1)

        TDay_avg = df_f['TDay'].mean(axis=1) #.values
        TempResponse_avg = df_f['TempResponse'].mean(axis=1) #.values
        df_avg = pd.concat([TDay_avg, TempResponse_avg], axis=1)
        df_avg.columns = ['TDay', 'TempResponse']
        df_avg['TminFactor'] = TminFactor
        df_avg['Tmin'] = Tmin
        df_avg['Topt'] = Topt
        df_avg['Tmax'] = Tmax
        df_f = None
        del df_f
        _ = gc.collect()
        return df_avg
    # ------
    # Process each site
    if (verbose):
        print("Processing all sites...")
    for uid in df['UID'].unique(): #[:1]:
        df2 = df[df['UID']==uid].reset_index(drop=True)
        country = df2['country'].unique()[0]
        loc = df2['location'].unique()[0]
        loc_code = df2['loc_code'].unique()[0]
        cycle = df2['cycle'].unique()[0]
        rue = df2['RUE'].unique()[0]
        obsYield = df2['ObsYield'].mean()
        simYield = df2['SimYield'].mean()
        #print("processing {} - {}".format(uid, loc))
        #
        sel_cmb_wetf = df2.groupby(['UID', 'TminFactor', 'Tmin', 'Topt', 'Tmax'], as_index=False).agg({'first'})
        for idx in sel_cmb_wetf.index:
            uid=sel_cmb_wetf['UID'][idx]
            TminFactor=sel_cmb_wetf['TminFactor'][idx]
            Tmin=sel_cmb_wetf['Tmin'][idx]
            TOpt=sel_cmb_wetf['Topt'][idx]
            Tmax=sel_cmb_wetf['Tmax'][idx]
            site=sites[uid-1]
            # Display charts for each combinations in one figure
            xVals, yVals = generating_TDayValues_WETF(site=site, TminFactor=TminFactor, Tmin=Tmin, Topt=TOpt, Tmax=Tmax, 
                                                      isVPDStress=isVPDStress, roundVal=roundVal, 
                                                      maxTDay=maxTDay, verbose=verbose )
            # Create new DF
            df3 = pd.DataFrame({'TDay': xVals, 'TempResponse': yVals})
            df3['UID'] = uid
            df3['country'] = country
            df3['location'] = loc
            df3['loc_code'] = loc_code
            df3['cycle'] = cycle
            df3['RUE'] = rue
            df3['TminFactor'] = TminFactor
            df3['Tmin'] = Tmin
            df3['Topt'] = TOpt
            df3['Tmax'] = Tmax
            df3['ObsYield'] = obsYield
            df3['SimYield'] = simYield
            df_WETF_TDays_vs_TempResponse = pd.concat([df_WETF_TDays_vs_TempResponse, df3])
    #
    del df3
    df_WETF_TDays_vs_TempResponse['TDay'] = df_WETF_TDays_vs_TempResponse['TDay'].astype(float).round(2)
    df_WETF_TDays_vs_TempResponse['TempResponse'] = df_WETF_TDays_vs_TempResponse['TempResponse'].astype(float).round(3)
    
    # ---------------------------
    #if (verbose):
    #print("Calculating average value for all simulations...")    
    df_WETF_TDays_vs_TempResponse_mean_allSites = pd.DataFrame()
    sel_cmb_wetf = df_WETF_TDays_vs_TempResponse.groupby(['TminFactor', 'Tmin', 'Topt', 'Tmax'], as_index=False).agg({'ObsYield':'mean'})
    for idx in sel_cmb_wetf.index:
        TminFactor=sel_cmb_wetf['TminFactor'][idx]
        Tmin=sel_cmb_wetf['Tmin'][idx]
        Topt=sel_cmb_wetf['Topt'][idx]
        Tmax=sel_cmb_wetf['Tmax'][idx]
        df_avg_tmfctor_topt = get_mean_TDays_vs_TempResponse_WETF(df=df_WETF_TDays_vs_TempResponse, 
                                                                  isVPDStress=isVPDStress, TminFactor=TminFactor, 
                                                                  Tmin=Tmin, Topt=Topt, Tmax=Tmax)
        df_WETF_TDays_vs_TempResponse_mean_allSites = pd.concat([df_WETF_TDays_vs_TempResponse_mean_allSites, 
                                                                 df_avg_tmfctor_topt], axis=0)

    if (saveTable):
        hoy = dt.datetime.now().strftime('%Y%m%d')
        #path_to_save_results = os.path.join(path_to_save_results, '{}_{}'.format('Data_for_Figures', hoy))
        if not os.path.isdir(path_to_save_results):
            os.makedirs(path_to_save_results)
        if (fmt=='csv'):
            df_WETF_TDays_vs_TempResponse.to_csv(os.path.join(path_to_save_results , f'WETF_dataForCharting_TDays_vs_TempResponse_{hoy}.csv'), index=False)
            df_WETF_TDays_vs_TempResponse_mean_allSites.to_csv(os.path.join(path_to_save_results , f'WETF_dataForCharting_TDays_vs_TempResponse_mean_allSites_{hoy}.csv'), index=False)
        elif (fmt=='parquet'):
            df_WETF_TDays_vs_TempResponse.to_parquet(os.path.join(path_to_save_results, f'WETF_dataForCharting_TDays_vs_TempResponse_{hoy}.parquet'), index=False)
            df_WETF_TDays_vs_TempResponse_mean_allSites.to_parquet(os.path.join(path_to_save_results , f'WETF_dataForCharting_TDays_vs_TempResponse_mean_allSites_{hoy}.parquet'), index=False)
    
    return df_WETF_TDays_vs_TempResponse, df_WETF_TDays_vs_TempResponse_mean_allSites


def diplay_Figure_Type_II_WETF_InOneFig(cmb=None, fnct='WETF', df_tdays=None, df_tdays_mean=None,
                                       saveFig=True, showFig=True, fmt='jpg', leg_ncol=3,
                                       path_to_save_results='./', fname=None):

    def createFigTDay_vs_TempResponseBySiteYrs_InOneFig_WETF(df_dailyVals=None, uid=1, fnct='WETF', ax=None):
        '''
            Temperature response type II function for heading to maturity
            Display TDay (°C) vs Relative temperature response for all sites
        '''
        # Get combinations parameters
        df = df_dailyVals[df_dailyVals['UID']==uid].reset_index(drop=True).sort_values(['TDay'])
        # Unique combinations by site
        # Unique combinations by site
        sel_cmb = df.groupby(['UID', 'TminFactor', 'Tmin', 'Topt', 'Tmax'], as_index=False).agg({'first'})
        # Create figure
        #print("Processing curves for each site...")
        for idx in sel_cmb.index:
            uid=sel_cmb['UID'][idx]
            TminFactor=sel_cmb['TminFactor'][idx]
            Tmin=sel_cmb['Tmin'][idx]
            TOpt=sel_cmb['Topt'][idx]
            Tmax=sel_cmb['Tmax'][idx]
            df_t = df[((df['TminFactor']==TminFactor) & (df['Tmin']==Tmin) & (df['Topt']==TOpt) & (df['Tmax']==Tmax) )] 
            x=df_t['TDay'].to_numpy()
            y=df_t['TempResponse'].to_numpy()
            lbl = None if (len(df['UID'].unique())==1) else 'Combinations for all sites'
            ax.plot(np.arange(len(y)), y, marker='o', ls='--', linewidth=0.15, markersize=0.01, 
                    color='gainsboro', zorder=0, label=lbl)
            if (idx == (len(sel_cmb)-1)):
                # Just to get the corrected x-axis labels 
                sns.pointplot(x=x, y=y, color=None, linestyles='', ax=ax, markers='')
        
    #
    def createFigAverageBySiteYrs_InOneFig_WETF(df_dailyVals=None, fnct='WETF', ax=None):
        '''
            Temperature response type II function for heading to maturity
            Display Average TDay (°C) vs Average Relative temperature response for all sites
        '''
        # Get combinations parameters
        df = df_dailyVals.copy()
        # Unique combinations by site
        sel_cmb = df.groupby(['TminFactor', 'Tmin', 'Topt', 'Tmax'], as_index=False).agg({'first'})
        # Create figure
        palette = sns.color_palette("Set1", len(sel_cmb))
        line = itertools.cycle(['-', '--', '-.', ':', 'solid', 'dashed', 'dashdot', 'dotted'])
        lines = [next(line) for i in range(len(sel_cmb))]
        marker = itertools.cycle(['o', '^', '*', '8', 's', 'p', 'd', 'v'])
        markers = [next(marker) for i in range(len(sel_cmb))]
        #print("Processing average curves for all site...")
        for idx in tqdm(sel_cmb.index):
            TminFactor=sel_cmb['TminFactor'][idx]
            Tmin=sel_cmb['Tmin'][idx]
            TOpt=sel_cmb['Topt'][idx]
            Tmax=sel_cmb['Tmax'][idx]
            df_t = df[((df['TminFactor']==TminFactor) & (df['Tmin']==Tmin) & (df['Topt']==TOpt) 
                       & (df['Tmax']==Tmax) )].reset_index(drop=True).sort_values(['TDay'])
            x=df_t['TDay'].to_numpy().round(3)
            y=df_t['TempResponse'].to_numpy()
            ax.plot(np.arange(len(y)), y, marker=markers[idx], ls=lines[idx], linewidth=0.25, markersize=0.01, 
                    color=palette[idx], zorder=0, 
                    label='TminFactor: {:.2f} - Tmin: {:.0f}°C - TOpt: {:.0f}°C - Tmax: {:.0f}°C'.format(TminFactor, Tmin, TOpt, Tmax))
            if (idx == (len(sel_cmb)-1)):
                # Just to get the corrected x-axis labels 
                sns.pointplot(x=x, y=y,  color=None, linestyles='', ax=ax, markers='')
        # 
    # -----------------------

    # Create Figure
    df = cmb.copy()
    
    # ----------
    fig, ax = plt.subplots(1, 1, figsize=(8,6), facecolor='white')
    handout=[]
    lablout=[]

    # All combinations
    # Average combinations
    if (len(df['UID'].unique())>1):
        for nfig in tqdm(range(len(df['UID'].unique()))):
            df2 = df[df['UID']==nfig+1].reset_index(drop=True)
            createFigTDay_vs_TempResponseBySiteYrs_InOneFig_WETF(df_tdays, uid=nfig+1, fnct=fnct, ax=ax)
    #
    # Average combinations
    createFigAverageBySiteYrs_InOneFig_WETF(df_tdays_mean, fnct=fnct, ax=ax)
    #
    plt.setp(ax.get_lines(),linewidth=0.75)  # set lw for all lines of g axes
    plt.setp(ax.collections, sizes=[0.01])
    
    # Agregar Optimun temperature guide lines 
    v_lines = []
    v_lines_values = ','.join([str(x) for x in df['Topt'].unique()]).split(',')
    for label in ax.get_xticklabels(): #minor=True
        for v in v_lines_values:
            if ( round(float(v),1) == round(float(label.get_text()),1) and '{:,.0f}'.format(float(v)) not in v_lines):
                if (len(v_lines) < len(v_lines_values)):
                    v_lines.append(label.get_position()[0])
                    
    ax.vlines(v_lines, 0, 1, ls='-.', color='lightgray', linewidth=0.75) #, label='Optimum Temperature')

    ax.set_xticks(ax.get_xticks()[::50]) # Number of TDays
    ax.tick_params(axis='x', labelsize=10, color='lightgray', rotation=90)
    ax.tick_params(axis='y', labelsize=10, color='lightgray')
    ax.grid(color='gainsboro', linestyle='-', linewidth=0.5, zorder=-1)
    ax.set_axisbelow(True)
    
    # Temperature response of type II (three cardinal temperature) function for heading to maturity across all sites.
    # Validate if only one site
    fname = "Figure_WETF_SiteYrs_Comparison" if fname is None else fname
    if (len(df['UID'].unique())==1):
        country = str(df.loc[0, 'country'])
        loc = str(df.loc[0, 'location'])
        cycle = str(df.loc[0, 'cycle'])
        title1='Temperature response function Type II ({})\n{} - {} - {}'.format(fnct, country, loc, cycle)
        fname = "Figure_WETF_Site_{}_{}_{}_Comparison".format(loc, country, cycle)
    else:
        title1='Average temperature response function Type II ({})'.format(fnct)
    # Title for the complete figure
    fig.suptitle('{}'.format(title1), fontsize='x-large', fontweight='bold' ) #, fontsize=20, y=1.005)
    ax.set_xlabel('TDay (°C)', fontweight='bold', fontsize=12)
    ax.set_ylabel('Relative response', fontweight='bold', fontsize=12)
    
    def getLegend_HandlesLabels(ax, handout, lablout):
        handles, labels = ax.get_legend_handles_labels()
        for h,l in zip(handles,labels):
            if l not in lablout:
                lablout.append(l)
                handout.append(h)
        return handout, lablout
    
    handout, lablout = getLegend_HandlesLabels(ax, handout, lablout)
    fig.legend(handout,lablout, bbox_to_anchor=(0.5, -0.25), loc="lower center", ncol=leg_ncol, 
               numpoints=1, markerscale=2,
               borderaxespad=0,bbox_transform=plt.gcf().transFigure, fontsize=10 ) #, fancybox=True, shadow=True)

    fig.tight_layout()
    hoy = dt.datetime.now().strftime('%Y%m%d')
    #figures_path = os.path.join(path_to_save_results, 'Figures') #fnct) #'{}_{}'.format(dirname, hoy))
    figures_path = os.path.join(path_to_save_results, '{}_{}'.format('Figures', hoy))
    if not os.path.isdir(figures_path):
        os.makedirs(figures_path)
    if (saveFig==True and (fmt=='jpg' or fmt=='png')):
        fig.savefig(os.path.join(figures_path, "{}_{}.{}".format(fname, hoy,fmt)), 
                    bbox_inches='tight', transparent=False, dpi=300)
    if (saveFig is True and fmt=='pdf'):
        fig.savefig(os.path.join(figures_path, "{}_{}.pdf".format(fname, hoy)), 
                    bbox_inches='tight', orientation='portrait', edgecolor='none', transparent=False, 
                    pad_inches=0.5, dpi=300)

    if (showFig is True):
        fig.show()
    else:
        del fig
        plt.close();
#

def display_FigTDay_vs_TempResponseBySiteYrs_WETF(cmb=None, df_tdays=None, roundVal=2, maxTDay=50, fnct='WETF', 
                                             saveFig=True, showFig=True, fmt='jpg', cols=5, leg_ncol=4,
                                             path_to_save_results='./'):
    df = cmb.copy()
                                       
    def createFigTDay_vs_TempResponseBySiteYrs_WETF(df_cmb=None, df_dailyVals=None, uid=1, fnct='WETF', ax=None):
        '''
            Temperature response type II function for heading to maturity
            Display TDay (°C) vs Relative temperature response by site
        '''
        
        # Get site data
        df0 = df_cmb[df_cmb['UID']==uid].reset_index(drop=True)
        country = df0['country'].unique()[0]
        loc = df0['location'].unique()[0]
        loc_code = df0['loc_code'].unique()[0]
        cycle = df0['cycle'].unique()[0]
        obsYield = df0['ObsYield'].mean()
        simYield = df0['SimYield'].mean()

        # Get combinations parameters
        df = df_dailyVals[df_dailyVals['UID']==uid].reset_index(drop=True).sort_values(['TDay'])
        #df = df.drop_duplicates(subset=['TDay', 'TempResponse']) # this produce strange plots
        # Unique combinations by site
        sel_cmb = df.groupby(['UID', 'TminFactor', 'Tmin', 'Topt', 'Tmax'], as_index=False).agg({'ObsYield':'mean'})
        # Create figure
        #fig, ax = plt.subplots(figsize=(8,6))
        palette = sns.color_palette("Set1", len(sel_cmb))
        line = itertools.cycle(['-', '--', '-.', ':', 'solid', 'dashed', 'dashdot', 'dotted'])
        lines = [next(line) for i in range(len(sel_cmb))]
        marker = itertools.cycle(['o', '^', '*', '8', 's', 'p', 'd', 'v'])
        markers = [next(marker) for i in range(len(sel_cmb))]
        for idx in sel_cmb.index:
            uid=sel_cmb['UID'][idx]
            TminFactor=sel_cmb['TminFactor'][idx]
            Tmin=sel_cmb['Tmin'][idx]
            TOpt=sel_cmb['Topt'][idx]
            Tmax=sel_cmb['Tmax'][idx]
            df_t = df[((df['TminFactor']==TminFactor) & (df['Tmin']==Tmin) & (df['Topt']==TOpt) & (df['Tmax']==Tmax) )].reset_index(drop=True).sort_values(['TDay'])
            x=df_t['TDay'].to_numpy().round(3)
            y=df_t['TempResponse'].to_numpy()
            ax.plot(np.arange(len(y)), y, marker=markers[idx], ls=lines[idx], linewidth=0.25, markersize=0.01, 
                    color=palette[idx], zorder=0, 
                    label='TminFactor: {:.2f} - Tmin: {:.0f}°C - TOpt: {:.0f}°C - Tmax: {:.0f}°C'.format(TminFactor, Tmin, TOpt, Tmax))
            if (idx == (len(sel_cmb)-1)):
                # Just to get the corrected x-axis labels 
                sns.pointplot(x=x, y=y,  color=None, linestyles='', ax=ax, markers='')
        # 
        plt.setp(ax.get_lines(),linewidth=0.75)  # set lw for all lines of g axes
        plt.setp(ax.collections, sizes=[0.01])
        v_lines = []
        v_lines_values = ','.join([str(x) for x in df0['Topt'].unique()]).split(',')
        for label in ax.get_xticklabels(): #minor=True
            for v in v_lines_values:
                if ( round(float(v),1) == round(float(label.get_text()),1) and '{:,.0f}'.format(float(v)) not in v_lines):
                    if (len(v_lines) < len(v_lines_values)):
                        v_lines.append(label.get_position()[0])

        ax.vlines(v_lines, 0, 1, ls='-.', color='lightgray', linewidth=0.75) #, label='Optimum Temperature (15, 16, 17°C)')
        # ------------
        ax.set_xticks(ax.get_xticks()[::50])
        ax.tick_params(axis='x', labelsize=6, color='lightgray', rotation=90)
        ax.tick_params(axis='y', labelsize=6, color='lightgray') #, rotation=90)

        ax.set_title('{} - {} - {} - {}'.format(uid, loc, loc_code, cycle)) #, fontweight='bold', fontsize=16)
        #ax.set_xlabel('TDay (°C)', fontweight='bold', fontsize=14)
        #ax.set_ylabel('Relative temperature response (°C)', fontweight='bold', fontsize=14)
        ax.set_xlabel('')  
        ax.set_ylabel('')
        ax.grid(color='gainsboro', linestyle='-', linewidth=0.5, zorder=-1)
        ax.set_axisbelow(True)


    def getLegend_HandlesLabels(ax, handout, lablout):
        handles, labels = ax.get_legend_handles_labels()
        for h,l in zip(handles,labels):
            if l not in lablout:
                lablout.append(l)
                handout.append(h)
        return handout, lablout
    
    # ------------------------------
    # Create Figure
    fig = plt.figure(figsize=(14, 20), facecolor='white', tight_layout=True) #, constrained_layout=True) 
    #cols = 5
    rows = int(np.ceil(len(df['UID'].unique()) / cols))
    gs = gridspec.GridSpec(nrows=rows, ncols=cols) #, wspace=0, hspace=0)

    handout=[]
    lablout=[]
    nfig = 0
    for nr in tqdm(range(rows)):
        for nc in range(cols):
            #if (nfig>0):
            #    ax = fig.add_subplot(gs[nr, nc], sharey=ax, sharex=ax)
            #else:
            ax = fig.add_subplot(gs[nr, nc])
            # Create figure
            #print("Processing site-year: ", nfig+1)
            createFigTDay_vs_TempResponseBySiteYrs_WETF(df, df_tdays, uid=nfig+1, fnct=fnct, ax=ax)
            # Legend
            handout, lablout = getLegend_HandlesLabels(ax, handout, lablout)
            nfig = nfig + 1

    # Temperature response of type I (one cardinal temperature) function for heading to maturity across all sites.
    fig.text(0.5, -0.01, 'TDay (°C)', ha='center', va='center', fontweight='bold', fontsize=14)
    fig.text(-0.02, 0.5, 'Relative temperature response', ha='center', va='center', 
                     rotation='vertical', fontweight='bold', fontsize=14)
    #title1='Temperature response WETF function for 50 site-years around the globe' #.format()
    # Title for the complete figure
    #fig.suptitle('{}'.format(title1), fontsize='x-large', fontweight='bold' ) #, fontsize=20, y=1.005)

    def updatescatter(handle, orig):
        handle.update_from(orig)
        handle.set_sizes([4])

    def updateline(handle, orig):
        handle.update_from(orig)
        handle.set_markersize(1)

    fig.legend(handout,lablout, bbox_to_anchor=(0.5, -0.08), loc="lower center", ncol=leg_ncol, 
               numpoints=1, markerscale=2,
               borderaxespad=0,bbox_transform=plt.gcf().transFigure, fontsize=12, 
               handler_map={PathCollection : HandlerPathCollection(update_func=updatescatter),
                            plt.Line2D : HandlerLine2D(update_func = updateline)}
              ) #, fancybox=True, shadow=True)

    fig.tight_layout()
    hoy = dt.datetime.now().strftime('%Y%m%d')
    #figures_path = os.path.join(path_to_save_results, fnct) #'{}_{}'.format(dirname, hoy))
    figures_path = os.path.join(path_to_save_results, '{}_{}'.format('Figures', hoy))
    if not os.path.isdir(figures_path):
        os.makedirs(figures_path)
    if (saveFig==True and (fmt=='jpg' or fmt=='png')):
        fig.savefig(os.path.join(figures_path, "Figure_WETF_BySiteYrs_{}.{}".format(hoy,fmt)), bbox_inches='tight',
                    transparent=False, dpi=300)
    if (saveFig is True and fmt=='pdf'):
        fig.savefig(os.path.join(figures_path, "Figure_WETF_BySiteYrs_{}.pdf".format(hoy)), 
                    bbox_inches='tight', orientation='portrait', edgecolor='none', transparent=False, 
                    pad_inches=0.5, dpi=300)

    if (showFig is True):
        fig.show()
    else:
        del fig
        plt.close();
        
#

# --------------------------------------
# TPF
# --------------------------------------
'''
    Create values for Tday used in detailed figures
    Estimate TDay and Temperature response type III function (TPF) values from heading to maturity
'''
def createFigure_Type_III_TPF(sites, cmb=None, isVPDStress=False, roundVal=3, maxTDay=50, saveTable=False, 
                             path_to_save_results='./', fmt='parquet', verbose=False):
    #
    # Setup table with all values for displaying all locations in one figure 
    df = cmb.copy()
    nfig = 0
    df_TPF_TDays_vs_TempResponse = pd.DataFrame()
    
    def generating_TDayValues_TPF(site=None, RUE=3.0, TminFactor=0.25, Tmin=6, Toptmin=10, Toptmax=20, Tmax=24, 
                                   isVPDStress=False, roundVal=1, maxTDay=40, verbose=False ):
        if (TminFactor is None or Toptmin is None or Toptmax is None):
            print("Parameters not found. Please check the inputs")
            return

        if (verbose):
            print("Getting TDay values for site")
        # Get TDays values
        tDay_site0 = tday.estimate_TDay(site.inputWPN['TMIN'].to_numpy(), 
                                site.inputWPN['TMAX'].to_numpy(), TminFactor )
        #print("tDay_site0", tDay_site0)
        tDay_site = [x/10 for x in range(0, maxTDay*10)] # Create curve values
        tDay_site = sorted(np.unique(np.concatenate([tDay_site, tDay_site0]))) #.round(roundVal)
        
        # Get TPF
        TPFTMAX = tpf.apply_TPF(tDay_site, Tmin, Toptmin, Toptmax, Tmax)
        arrTPF = []
        for x in TPFTMAX: # Clean values lower than zero
            if ((x > 0.0) and (x <= 1.0) ):
                arrTPF.append(x)
            elif x > 1.0:
                arrTPF.append(1.0)
            else:
                arrTPF.append(0.0)

        return np.array(tDay_site).round(roundVal), np.array(arrTPF).round(roundVal)
    
    #
    def get_mean_TDays_vs_TempResponse_TPF(df, isVPDStress=False, TminFactor=0.25, Tmin=9.0, Toptmin=10, 
                                           Toptmax=20, Tmax=32):
        df_f = pd.DataFrame()
        for uid in df['UID'].unique():
            df_1 = df[( (df['UID']==uid)  & (df['TminFactor']==TminFactor) 
                       & (df['Tmin']==Tmin) & (df['Toptmin']==Toptmin) & (df['Toptmax']==Toptmax) 
                       & (df['Tmax']==Tmax)
            )][['UID', 'TDay', 'TempResponse']].reset_index(drop=True)
            df_f = pd.concat([df_f, df_1], axis=1)

        TDay_avg = df_f['TDay'].mean(axis=1) #.values
        TempResponse_avg = df_f['TempResponse'].mean(axis=1) #.values
        df_avg = pd.concat([TDay_avg, TempResponse_avg], axis=1)
        df_avg.columns = ['TDay', 'TempResponse']
        df_avg['TminFactor'] = TminFactor
        df_avg['Tmin'] = Tmin
        df_avg['Toptmin'] = Toptmin
        df_avg['Toptmax'] = Toptmax
        df_avg['Tmax'] = Tmax
        df_f = None
        del df_f
        _ = gc.collect()
        return df_avg
    # ------
    # Process each site
    if (verbose):
        print("Processing all sites...")
    for uid in df['UID'].unique(): #[:1]:
        df2 = df[df['UID']==uid].reset_index(drop=True)
        country = df2['country'].unique()[0]
        loc = df2['location'].unique()[0]
        loc_code = df2['loc_code'].unique()[0]
        cycle = df2['cycle'].unique()[0]
        rue = df2['RUE'].unique()[0]
        obsYield = df2['ObsYield'].mean()
        simYield = df2['SimYield'].mean()
        #print("processing {} - {}".format(uid, loc))
        #
        sel_cmb_tpf = df2.groupby(['UID', 'TminFactor', 'Tmin', 'Toptmin', 'Toptmax', 'Tmax'], 
                                  as_index=False).agg({'first'})
        for idx in sel_cmb_tpf.index:
            uid=sel_cmb_tpf['UID'][idx]
            TminFactor=sel_cmb_tpf['TminFactor'][idx]
            Tmin=sel_cmb_tpf['Tmin'][idx]
            Toptmin=sel_cmb_tpf['Toptmin'][idx]
            Toptmax=sel_cmb_tpf['Toptmax'][idx]
            Tmax=sel_cmb_tpf['Tmax'][idx]
            site=sites[uid-1]
            # Display charts for each combinations in one figure
            xVals, yVals = generating_TDayValues_TPF(site=site, TminFactor=TminFactor, Tmin=Tmin, 
                                                     Toptmin=Toptmin, Toptmax=Toptmax, Tmax=Tmax, 
                                                     isVPDStress=isVPDStress, roundVal=roundVal, 
                                                     maxTDay=maxTDay, verbose=verbose )
            # Create new DF
            df3 = pd.DataFrame({'TDay': xVals, 'TempResponse': yVals})
            df3['UID'] = uid
            df3['country'] = country
            df3['location'] = loc
            df3['loc_code'] = loc_code
            df3['cycle'] = cycle
            df3['RUE'] = rue
            df3['TminFactor'] = TminFactor
            df3['Tmin'] = Tmin
            df3['Toptmin'] = Toptmin
            df3['Toptmax'] = Toptmax
            df3['Tmax'] = Tmax
            df3['ObsYield'] = obsYield
            df3['SimYield'] = simYield
            df_TPF_TDays_vs_TempResponse = pd.concat([df_TPF_TDays_vs_TempResponse, df3])
    #
    del df3
    df_TPF_TDays_vs_TempResponse['TDay'] = df_TPF_TDays_vs_TempResponse['TDay'].astype(float).round(2)
    df_TPF_TDays_vs_TempResponse['TempResponse'] = df_TPF_TDays_vs_TempResponse['TempResponse'].astype(float).round(3)
    
    # ---------------------------
    #if (verbose):
    #print("Calculating average value for all simulations...")    
    df_TPF_TDays_vs_TempResponse_mean_allSites = pd.DataFrame()
    sel_cmb_tpf = df_TPF_TDays_vs_TempResponse.groupby(['TminFactor', 'Tmin', 'Toptmin', 'Toptmax', 'Tmax'], 
                                                       as_index=False).agg({'ObsYield':'mean'})
    for idx in sel_cmb_tpf.index:
        TminFactor=sel_cmb_tpf['TminFactor'][idx]
        Tmin=sel_cmb_tpf['Tmin'][idx]
        Toptmin=sel_cmb_tpf['Toptmin'][idx]
        Toptmax=sel_cmb_tpf['Toptmax'][idx]
        Tmax=sel_cmb_tpf['Tmax'][idx]
        df_avg_tmfctor_topt = get_mean_TDays_vs_TempResponse_TPF(df=df_TPF_TDays_vs_TempResponse, 
                                                                  isVPDStress=isVPDStress, TminFactor=TminFactor, 
                                                                  Tmin=Tmin, Toptmin=Toptmin, Toptmax=Toptmax, 
                                                                 Tmax=Tmax)
        df_TPF_TDays_vs_TempResponse_mean_allSites = pd.concat([df_TPF_TDays_vs_TempResponse_mean_allSites, 
                                                                 df_avg_tmfctor_topt], axis=0)

    if (saveTable):
        hoy = dt.datetime.now().strftime('%Y%m%d')
        #path_to_save_results = os.path.join(path_to_save_results, '{}_{}'.format('Data_for_Figures', hoy))
        if not os.path.isdir(path_to_save_results):
            os.makedirs(path_to_save_results)
        if (fmt=='csv'):
            df_TPF_TDays_vs_TempResponse.to_csv(os.path.join(path_to_save_results , f'TPF_dataForCharting_TDays_vs_TempResponse_{hoy}.csv'), index=False)
            df_TPF_TDays_vs_TempResponse_mean_allSites.to_csv(os.path.join(path_to_save_results , f'TPF_dataForCharting_TDays_vs_TempResponse_mean_allSites_{hoy}.csv'), index=False)
        elif (fmt=='parquet'):
            df_TPF_TDays_vs_TempResponse.to_parquet(os.path.join(path_to_save_results, f'TPF_dataForCharting_TDays_vs_TempResponse_{hoy}.parquet'), index=False)
            df_TPF_TDays_vs_TempResponse_mean_allSites.to_parquet(os.path.join(path_to_save_results , f'TPF_dataForCharting_TDays_vs_TempResponse_mean_allSites_{hoy}.parquet'), index=False)
    
    return df_TPF_TDays_vs_TempResponse, df_TPF_TDays_vs_TempResponse_mean_allSites


# -----------
def diplay_Figure_Type_II_TPF_InOneFig(cmb=None, fnct='TPF', df_tdays=None, df_tdays_mean=None,
                                       saveFig=True, showFig=True, fmt='jpg', leg_ncol=3,
                                       path_to_save_results='./', fname=None):

    def createFigTDay_vs_TempResponseBySiteYrs_InOneFig_TPF(df_dailyVals=None, uid=1, fnct='TPF', ax=None):
        '''
            Temperature response type III function for heading to maturity
            Display TDay (°C) vs Relative temperature response for all sites
        '''
        # Get combinations parameters
        df = df_dailyVals[df_dailyVals['UID']==uid].reset_index(drop=True).sort_values(['TDay'])
        # Unique combinations by site
        sel_cmb = df.groupby(['UID', 'TminFactor', 'Tmin', 'Toptmin', 'Toptmax', 'Tmax'], as_index=False).agg({'first'})
        # Create figure
        #print("Processing curves for each site...")
        for idx in sel_cmb.index:
            uid=sel_cmb['UID'][idx]
            TminFactor=sel_cmb['TminFactor'][idx]
            Tmin=sel_cmb['Tmin'][idx]
            Toptmin=sel_cmb['Toptmin'][idx]
            Toptmax=sel_cmb['Toptmax'][idx]
            Tmax=sel_cmb['Tmax'][idx]
            df_t = df[((df['TminFactor']==TminFactor) & (df['Tmin']==Tmin) & (df['Toptmin']==Toptmin) 
                       & (df['Toptmax']==Toptmax) & (df['Tmax']==Tmax) )] 
            x=df_t['TDay'].to_numpy()
            y=df_t['TempResponse'].to_numpy()
            lbl = None if (len(df['UID'].unique())==1) else 'Combinations for all sites'
            ax.plot(np.arange(len(y)), y, marker='o', ls='--', linewidth=0.15, markersize=0.01, 
                    color='gainsboro', zorder=0, label=lbl)
            if (idx == (len(sel_cmb)-1)):
                # Just to get the corrected x-axis labels 
                sns.pointplot(x=x, y=y, color=None, linestyles='', ax=ax, markers='')
        
    #
    def createFigAverageBySiteYrs_InOneFig_TPF(df_dailyVals=None, fnct='TPF', ax=None):
        '''
            Temperature response type II function for heading to maturity
            Display Average TDay (°C) vs Average Relative temperature response for all sites
        '''
        # Get combinations parameters
        df = df_dailyVals.copy()
        # Unique combinations by site
        sel_cmb = df.groupby(['TminFactor', 'Tmin', 'Toptmin', 'Toptmax', 'Tmax'], as_index=False).agg({'first'})
        # Create figure
        palette = sns.color_palette("Set1", len(sel_cmb))
        line = itertools.cycle(['-', '--', '-.', ':', 'solid', 'dashed', 'dashdot', 'dotted'])
        lines = [next(line) for i in range(len(sel_cmb))]
        marker = itertools.cycle(['o', '^', '*', '8', 's', 'p', 'd', 'v'])
        markers = [next(marker) for i in range(len(sel_cmb))]
        #print("Processing average curves for all site...")
        for idx in tqdm(sel_cmb.index):
            TminFactor=sel_cmb['TminFactor'][idx]
            Tmin=sel_cmb['Tmin'][idx]
            Toptmin=sel_cmb['Toptmin'][idx]
            Toptmax=sel_cmb['Toptmax'][idx]
            Tmax=sel_cmb['Tmax'][idx]
            df_t = df[((df['TminFactor']==TminFactor) & (df['Tmin']==Tmin) & (df['Toptmin']==Toptmin) 
                       & (df['Toptmax']==Toptmax) & (df['Tmax']==Tmax) )] .reset_index(drop=True).sort_values(['TDay'])
            x=df_t['TDay'].to_numpy().round(3)
            y=df_t['TempResponse'].to_numpy()
            ax.plot(np.arange(len(y)), y, marker=markers[idx], ls=lines[idx], linewidth=0.25, markersize=0.01, 
                    color=palette[idx], zorder=0, 
                    label='TminFactor: {:.2f} - Tmin: {:.0f}°C - Toptmin: {:.0f}°C - Toptmax: {:.0f}°C - Tmax: {:.0f}°C'.format(TminFactor, Tmin, Toptmin, Toptmax, Tmax))
            if (idx == (len(sel_cmb)-1)):
                # Just to get the corrected x-axis labels 
                sns.pointplot(x=x, y=y,  color=None, linestyles='', ax=ax, markers='')
        # 
    # -----------------------

    # Create Figure
    df = cmb.copy()
    
    # ----------
    fig, ax = plt.subplots(1, 1, figsize=(8,6), facecolor='white')
    handout=[]
    lablout=[]

    # All combinations
    # Average combinations
    if (len(df['UID'].unique())>1):
        for nfig in tqdm(range(len(df['UID'].unique()))):
            df2 = df[df['UID']==nfig+1].reset_index(drop=True)
            createFigTDay_vs_TempResponseBySiteYrs_InOneFig_TPF(df_tdays, uid=nfig+1, fnct=fnct, ax=ax)
    #
    # Average combinations
    createFigAverageBySiteYrs_InOneFig_TPF(df_tdays_mean, fnct=fnct, ax=ax)
    #
    plt.setp(ax.get_lines(),linewidth=0.75)  # set lw for all lines of g axes
    plt.setp(ax.collections, sizes=[0.01])
    
    # Agregar Optimun temperature guide lines 
    v_lines = []
    v_lines_values = ','.join([str(x) for x in df['Toptmin'].unique()]).split(',')
    for label in ax.get_xticklabels(): #minor=True
        for v in v_lines_values:
            if ( round(float(v),1) == round(float(label.get_text()),1) and '{:,.0f}'.format(float(v)) not in v_lines):
                if (len(v_lines) < len(v_lines_values)):
                    v_lines.append(label.get_position()[0])
                    
    ax.vlines(v_lines, 0, 1, ls='-.', color='lightgray', linewidth=0.75) #, label='Optimum Temperature')

    ax.set_xticks(ax.get_xticks()[::50]) # Number of TDays
    ax.tick_params(axis='x', labelsize=10, color='lightgray', rotation=90)
    ax.tick_params(axis='y', labelsize=10, color='lightgray')
    ax.grid(color='gainsboro', linestyle='-', linewidth=0.5, zorder=-1)
    ax.set_axisbelow(True)
    
    # Temperature response of type II (three cardinal temperature) function for heading to maturity across all sites.
    # Validate if only one site
    fname = "Figure_TPF_SiteYrs_Comparison" if fname is None else fname
    if (len(df['UID'].unique())==1):
        country = str(df.loc[0, 'country'])
        loc = str(df.loc[0, 'location'])
        cycle = str(df.loc[0, 'cycle'])
        title1='Temperature response function Type III ({})\n{} - {} - {}'.format(fnct, country, loc, cycle)
        fname = "Figure_TPF_Site_{}_{}_{}_Comparison".format(loc, country, cycle)
    else:
        title1='Average temperature response function Type III ({})'.format(fnct)
    # Title for the complete figure
    fig.suptitle('{}'.format(title1), fontsize='x-large', fontweight='bold' ) #, fontsize=20, y=1.005)
    ax.set_xlabel('TDay (°C)', fontweight='bold', fontsize=12)
    ax.set_ylabel('Relative response', fontweight='bold', fontsize=12)
    
    def getLegend_HandlesLabels(ax, handout, lablout):
        handles, labels = ax.get_legend_handles_labels()
        for h,l in zip(handles,labels):
            if l not in lablout:
                lablout.append(l)
                handout.append(h)
        return handout, lablout
    
    handout, lablout = getLegend_HandlesLabels(ax, handout, lablout)
    fig.legend(handout,lablout, bbox_to_anchor=(0.5, -0.25), loc="lower center", ncol=leg_ncol, 
               numpoints=1, markerscale=2,
               borderaxespad=0,bbox_transform=plt.gcf().transFigure, fontsize=10 ) #, fancybox=True, shadow=True)

    fig.tight_layout()
    hoy = dt.datetime.now().strftime('%Y%m%d')
    #figures_path = os.path.join(path_to_save_results, 'Figures') #fnct) #'{}_{}'.format(dirname, hoy))
    figures_path = os.path.join(path_to_save_results, '{}_{}'.format('Figures', hoy))
    if not os.path.isdir(figures_path):
        os.makedirs(figures_path)
    if (saveFig==True and (fmt=='jpg' or fmt=='png')):
        fig.savefig(os.path.join(figures_path, "{}_{}.{}".format(fname, hoy,fmt)), 
                    bbox_inches='tight', transparent=False, dpi=300)
    if (saveFig is True and fmt=='pdf'):
        fig.savefig(os.path.join(figures_path, "{}_{}.pdf".format(fname, hoy)), 
                    bbox_inches='tight', orientation='portrait', edgecolor='none', transparent=False, 
                    pad_inches=0.5, dpi=300)

    if (showFig is True):
        fig.show()
    else:
        del fig
        plt.close();
# ---------
def display_FigTDay_vs_TempResponseBySiteYrs_TPF(cmb=None, df_tdays=None, roundVal=2, maxTDay=50, fnct='TPF', 
                                             saveFig=True, showFig=True, fmt='jpg', cols=5, leg_ncol=4,
                                             path_to_save_results='./'):
    df = cmb.copy()
                                       
    def createFigTDay_vs_TempResponseBySiteYrs_TPF(df_cmb=None, df_dailyVals=None, uid=1, fnct='TPF', ax=None):
        '''
            Temperature response type III function for heading to maturity
            Display TDay (°C) vs Relative temperature response by site
        '''
        
        # Get site data
        df0 = df_cmb[df_cmb['UID']==uid].reset_index(drop=True)
        country = df0['country'].unique()[0]
        loc = df0['location'].unique()[0]
        loc_code = df0['loc_code'].unique()[0]
        cycle = df0['cycle'].unique()[0]
        obsYield = df0['ObsYield'].mean()
        simYield = df0['SimYield'].mean()

        # Get combinations parameters
        df = df_dailyVals[df_dailyVals['UID']==uid].reset_index(drop=True).sort_values(['TDay'])
        # Unique combinations by site
        sel_cmb = df.groupby(['UID', 'TminFactor', 'Tmin', 'Toptmin', 'Toptmax', 'Tmax'], as_index=False).agg({'ObsYield':'mean'})
        # Create figure
        #fig, ax = plt.subplots(figsize=(8,6))
        palette = sns.color_palette("Set1", len(sel_cmb))
        line = itertools.cycle(['-', '--', '-.', ':', 'solid', 'dashed', 'dashdot', 'dotted'])
        lines = [next(line) for i in range(len(sel_cmb))]
        marker = itertools.cycle(['o', '^', '*', '8', 's', 'p', 'd', 'v'])
        markers = [next(marker) for i in range(len(sel_cmb))]
        for idx in sel_cmb.index:
            uid=sel_cmb['UID'][idx]
            TminFactor=sel_cmb['TminFactor'][idx]
            Tmin=sel_cmb['Tmin'][idx]
            Toptmin=sel_cmb['Toptmin'][idx]
            Toptmax=sel_cmb['Toptmax'][idx]
            Tmax=sel_cmb['Tmax'][idx]
            df_t = df[((df['TminFactor']==TminFactor) & (df['Tmin']==Tmin) & (df['Toptmin']==Toptmin) 
                       & (df['Toptmax']==Toptmax) & (df['Tmax']==Tmax) )].reset_index(drop=True).sort_values(['TDay'])
            x=df_t['TDay'].to_numpy().round(3)
            y=df_t['TempResponse'].to_numpy()
            ax.plot(np.arange(len(y)), y, marker=markers[idx], ls=lines[idx], linewidth=0.25, markersize=0.01, 
                    color=palette[idx], zorder=0, 
                    label='TminFactor: {:.2f} - Tmin: {:.0f}°C - Toptmin: {:.0f}°C - Toptmax: {:.0f}°C - Tmax: {:.0f}°C'.format(TminFactor, Tmin, Toptmin, Toptmax, Tmax))
            
            if (idx == (len(sel_cmb)-1)):
                # Just to get the corrected x-axis labels 
                sns.pointplot(x=x, y=y,  color=None, linestyles='', ax=ax, markers='')
        # 
        plt.setp(ax.get_lines(),linewidth=0.75)  # set lw for all lines of g axes
        plt.setp(ax.collections, sizes=[0.01])
        v_lines = []
        v_lines_values = ','.join([str(x) for x in df0['Toptmin'].unique()]).split(',')
        for label in ax.get_xticklabels(): #minor=True
            for v in v_lines_values:
                if ( round(float(v),1) == round(float(label.get_text()),1) and '{:,.0f}'.format(float(v)) not in v_lines):
                    if (len(v_lines) < len(v_lines_values)):
                        v_lines.append(label.get_position()[0])

        ax.vlines(v_lines, 0, 1, ls='-.', color='lightgray', linewidth=0.75) #, label='Optimum Temperature')
        # ------------
        ax.set_xticks(ax.get_xticks()[::50])
        ax.tick_params(axis='x', labelsize=6, color='lightgray', rotation=90)
        ax.tick_params(axis='y', labelsize=6, color='lightgray') #, rotation=90)

        ax.set_title('{} - {} - {} - {}'.format(uid, loc, loc_code, cycle)) #, fontweight='bold', fontsize=16)
        #ax.set_xlabel('TDay (°C)', fontweight='bold', fontsize=14)
        #ax.set_ylabel('Relative temperature response (°C)', fontweight='bold', fontsize=14)
        ax.set_xlabel('')  
        ax.set_ylabel('')
        ax.grid(color='gainsboro', linestyle='-', linewidth=0.5, zorder=-1)
        ax.set_axisbelow(True)


    def getLegend_HandlesLabels(ax, handout, lablout):
        handles, labels = ax.get_legend_handles_labels()
        for h,l in zip(handles,labels):
            if l not in lablout:
                lablout.append(l)
                handout.append(h)
        return handout, lablout
    
    # ------------------------------
    # Create Figure
    fig = plt.figure(figsize=(14, 20), facecolor='white', tight_layout=True) #, constrained_layout=True) 
    #cols = 5
    rows = int(np.ceil(len(df['UID'].unique()) / cols))
    gs = gridspec.GridSpec(nrows=rows, ncols=cols) #, wspace=0, hspace=0)

    handout=[]
    lablout=[]
    nfig = 0
    for nr in tqdm(range(rows)):
        for nc in range(cols):
            #if (nfig>0):
            #    ax = fig.add_subplot(gs[nr, nc], sharey=ax, sharex=ax)
            #else:
            ax = fig.add_subplot(gs[nr, nc])
            # Create figure
            #print("Processing site-year: ", nfig+1)
            createFigTDay_vs_TempResponseBySiteYrs_TPF(df, df_tdays, uid=nfig+1, fnct=fnct, ax=ax)
            # Legend
            handout, lablout = getLegend_HandlesLabels(ax, handout, lablout)
            nfig = nfig + 1

    # Temperature response of type I (one cardinal temperature) function for heading to maturity across all sites.
    fig.text(0.5, -0.01, 'TDay (°C)', ha='center', va='center', fontweight='bold', fontsize=14)
    fig.text(-0.02, 0.5, 'Relative temperature response', ha='center', va='center', 
                     rotation='vertical', fontweight='bold', fontsize=14)
    #title1='Temperature response TPF function for 50 site-years around the globe' #.format()
    # Title for the complete figure
    #fig.suptitle('{}'.format(title1), fontsize='x-large', fontweight='bold' ) #, fontsize=20, y=1.005)

    def updatescatter(handle, orig):
        handle.update_from(orig)
        handle.set_sizes([4])

    def updateline(handle, orig):
        handle.update_from(orig)
        handle.set_markersize(1)

    fig.legend(handout,lablout, bbox_to_anchor=(0.5, -0.15), loc="lower center", ncol=leg_ncol, 
               numpoints=1, markerscale=2,
               borderaxespad=0,bbox_transform=plt.gcf().transFigure, fontsize=12, 
               handler_map={PathCollection : HandlerPathCollection(update_func=updatescatter),
                            plt.Line2D : HandlerLine2D(update_func = updateline)}
              ) #, fancybox=True, shadow=True)

    fig.tight_layout()
    hoy = dt.datetime.now().strftime('%Y%m%d')
    #figures_path = os.path.join(path_to_save_results, fnct) #'{}_{}'.format(dirname, hoy))
    figures_path = os.path.join(path_to_save_results, '{}_{}'.format('Figures', hoy))
    if not os.path.isdir(figures_path):
        os.makedirs(figures_path)
    if (saveFig==True and (fmt=='jpg' or fmt=='png')):
        fig.savefig(os.path.join(figures_path, "Figure_TPF_BySiteYrs_{}.{}".format(hoy,fmt)), bbox_inches='tight',
                    transparent=False, dpi=300)
    if (saveFig is True and fmt=='pdf'):
        fig.savefig(os.path.join(figures_path, "Figure_TPF_BySiteYrs_{}.pdf".format(hoy)), 
                    bbox_inches='tight', orientation='portrait', edgecolor='none', transparent=False, 
                    pad_inches=0.5, dpi=300)

    if (showFig is True):
        fig.show()
    else:
        del fig
        plt.close();

# ---------
