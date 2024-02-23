# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Mon Nov 13 18:17:27 2023
"""
import matplotlib.pyplot as plt
import pandas as pd
from collections.abc import Iterable

def mean_range_plot(df,column=None, ax=None, alpha=0.3, label=None, **kwargs):
    
    ## if df is a dataframe:
    ##      x = df.index, y1=df[column1], y2=df[column2],...
    
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
    line,  = ax.plot(x_,mean_,label=label, **kwargs)
    line_color = line.get_color()
    ax.fill_between(x_,min_,max_,alpha=alpha,color=line_color)
    ## fake plot for legend
    # ax.plot([],[],label=label,**kwargs)

