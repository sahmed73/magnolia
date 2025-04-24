# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Nov 13 18:17:27 2023
"""
import matplotlib.pyplot as plt
import pandas as pd
from collections.abc import Iterable
import matplotlib as mpl
import random
import numpy as np

    
def mean_range_plot(df, column=None, ax=None, alpha=0.3,
                    label=None, meanplot=True, shade=True, **kwargs):
    
    ## if df is a dataframe:
    ##      x = df.index, y1=df[column1], y2=df[column2]...
    
    ## if df is a dictionary of dataframes with same indecies:
    ##      we need 'column' parameter
    ##      x = df[any_keys].index, y1=df[key1][column],...
    
    ## if df is a iteration of dataframes with same indices:
    ##      we need 'column' parameter
    ##      x = df[any_index].index, y1=df[index1][column],...
    
    ## plot: y_mean vs x AND range=(y_max-y_min) as shade
    
    ## The range shows the span from the lowest to the highest value,
    ## giving a clear picture of the total spread of the data. It is the
    ## simplest measure of variability. Range as shade. Mean as line.    
    
    ##
    if isinstance(df, pd.DataFrame):
        dff = df.copy()
        
    elif isinstance(df, dict):
        if column is None:
            raise ValueError("'column' is not defined!")
        
        dff = pd.DataFrame()
        for index, (key, value) in enumerate(df.items()):
            dff[f'column_{index}'] = value[column]
    
    elif isinstance(df, Iterable):
        if column is None:
            raise ValueError("'column' is not defined!")
        
        dff = pd.DataFrame()
        for index, value in enumerate(df):
            dff[f'column_{index}'] = value[column]
        
    else:
        raise TypeError("'df' type must be pandas OR dict of pandas")
    
    mean_ = dff.mean(axis=1)
    max_  = dff.max(axis=1)
    min_  = dff.min(axis=1)
    x_    = dff.index
    
    if ax is None: fig, ax = plt.subplots()
    
    ## ploting
    if meanplot:
        line,  = ax.plot(x_,mean_,label=label, **kwargs)
    else:
        line,  = ax.plot([],[], **kwargs)
        
    line_color = line.get_color()
    if shade:
        ax.fill_between(x_,min_,max_,alpha=alpha,color=line_color)
    ## fake plot for legend
    # ax.plot([],[],label=label,**kwargs)
    
def lineplot_mean_with_std_error(df, ax=None):
    """
    data is a dataframe like this:
             cutoff   Sim-1  Sim-2  Sim-3
              0.30     29     44     22
              0.35     27     43     21
              0.40     28     39     22
              0.45     26     38     22
              0.50     25     38     21
              0.55     24     37     25
              0.60     27     45     28
              0.65     41     55     42
              0.70     50     77     54
              0.75     78    120     85
    """
    if ax is None:
        fig, ax = plt.subplots()
    
    ax.errorbar(df.index, df.mean(axis=1), yerr=df.std(axis=1),
                fmt='-o', capsize=5)
    
    return ax

def categorized_bar_plot(categories, series_means, series_errors=None,
                         series_labels=None, bar_width=0.35,
                         series_gap=0.1, category_gap=0.5,
                         ax=None, **kwargs):
    """
    Plots a categorized bar plot with optional error bars for an indefinite number of series,
    with customizable gaps between series and categories.

    Parameters:
    - categories: List of category names on the x-axis.
    - series_means: List of lists, where each inner list contains the means for a series.
    - series_errors: List of lists, where each inner list contains the errors for a series (optional).
    - series_labels: List of series labels (optional). If None, labels will be 'Series 1', 'Series 2', etc.
    - bar_width: Width of the bars (optional, default is 0.35).
    - series_gap: Gap between each series of bars within a category (optional, default is 0.1).
    - category_gap: Gap between each category on the x-axis (optional, default is 0.5).
    - ax: Matplotlib axes object to plot on. If None, a new figure and axes are created.
    """
    n_categories = len(categories)
    n_series = len(series_means)
    
    # Set the positions of the bars on the x-axis, adjusting for the series and category gaps
    total_series_width = n_series * bar_width + (n_series - 1) * series_gap
    indices = np.arange(n_categories) * (total_series_width + category_gap)

    # Create the bar plots
    if ax is None:
        fig, ax = plt.subplots()
    
    for i in range(n_series):
        # Only add error bars if series_errors is provided
        if series_errors is not None:
            ax.bar(indices + i * (bar_width + series_gap), series_means[i], 
                   bar_width, yerr=series_errors[i],
                   label=series_labels[i] if series_labels else f'Series {i + 1}',
                   capsize=5, **kwargs)
        else:
            ax.bar(indices + i * (bar_width + series_gap), series_means[i], bar_width, 
                   label=series_labels[i] if series_labels else f'Series {i + 1}',
                   **kwargs)

    # Add labels and title
    ax.set_xlabel('Category')
    ax.set_ylabel('Scores')
    ax.set_xticks(indices + (total_series_width - bar_width) / 2)
    ax.set_xticklabels(categories)
    ax.legend()
    
    return ax

    
def R_square(y_exact,y_pred):
    y_mean    = y_exact.mean()
    SST       = ((y_exact-y_mean)**2).sum()
    SSR       =((y_exact-y_pred)**2).sum()
    r_squared = 1-SSR/SST
    return r_squared

def custom_plot_features(minorticks=False,font='arial',
                         fontsize_incr = 0):
    plt.style.use('classic')
    plt_width = 7
    aspect_ratio = 1.333333
    plt.figure(figsize=[plt_width,plt_width/aspect_ratio])
    plt.rcParams['font.family'] = font
    mpl.rcParams['axes.labelpad'] = 15
    plt.rc('axes', grid=True)
    plt.rc('font', size=20+fontsize_incr) # Default text sizes
    plt.rc('axes', titlesize=20+fontsize_incr) # Axes title size
    plt.rc('axes', labelsize=20+fontsize_incr) # Axes label size
    plt.rc('xtick', labelsize=20+fontsize_incr) # X-axis tick label size
    plt.rc('ytick', labelsize=20+fontsize_incr) # Y-axis tick label size
    plt.rc('legend', fontsize=20) # Legend fontsize
    if minorticks: plt.minorticks_on()
    else: plt.minorticks_off()
    

def saveplot(folder='figures', name='plot', filetype='png',
             dpi=350, trans=False):
    rand = random.randint(0, 99999999999999)
    dirr = r'C:\Users\arup2\OneDrive - University of California Merced\Desktop\LAMMPS\Post Processing\lammps_processing\python_outputs'
    plt.savefig(dirr+'\\'+folder+'\\'+name+f'-{rand}.{filetype}',
                dpi=dpi, bbox_inches='tight', transparent=trans)