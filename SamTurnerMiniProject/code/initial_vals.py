""" Generate initial values for Non Linear Least Squares model fitting """

__appname__ = 'initial_vals.py'
__author__ = 'Sam Turner (sat19@ic.ac.uk)'
__version__ = '0.0.1'
__license__ = 'GNU public' 


# IMPORTS

import pandas as pd
import os
import matplotlib.pyplot as plt
from collections import Counter
import scipy.stats as stats
import numpy as np
import math

# load data
data = pd.read_csv("../data/LogisticGrowthDataLogClean.csv",index_col=0, sep="\t")
# get unique IDs
IDs = set(data["ID"])

# get unique parameter combinations
combos = list(set(zip(data["ID"],data["Species"],data["Medium"],data["Temp"],data["Rep"])))

# ID --> parameter combination dictionary
combos_dic = {l[0]:l[1:] for l in combos}

# get r for gompertz, baranyi, and buchanan
# this is the maximum gradient of the log10 time series
def log_inits(id, plot_fig):
    """Calculates initial parameter estimates by rolling regression in logarithmic space for specified ID. Plots a figure displaying the rolling regression if specified
    
    PARAMETERS
    ----------
    id : int
        id number of data for which parameter estiamtes are needed.

    plot_fig : bool
        boolean indiacting whether to plot rolling regression figure

    RETURNS
    -------
    list
        list of intial parameter estimates [Nmin, Nmax, mu_max, t_lag]"""
    
    # get time series for this ID
    
    time_series = data[data["ID"] == id].sort_values("Time")
    time_series.index = range(0,time_series.shape[0]) 

    # get length of time series to calculate number of rolling regression lines
    n_rolls = time_series.shape[0]

    # initialise valeus
    r_max = 0
    tlag_est = 0
    N0 = min(time_series["logPopBio"])

    # plot sccatter of time series
    if plot_fig:
        fig, ax = plt.subplots()
        ax.scatter(time_series["Time"], time_series["logPopBio"])

    # calculate a smaller chunk size for rolling regression if there are few points
    chunk_size = min(int(np.ceil(time_series.shape[0]/4)),5)
    chunk_size = max(chunk_size,3)

    # specify number of regressions
    for i in range(0,n_rolls-chunk_size-1):
        # calculate appropriate subset
        subset = time_series.iloc[i:i+chunk_size,:]
        # fit linear model
        fit = stats.linregress(subset["Time"],subset["logPopBio"])
        # if we have found a new highest gradient, store it and the tlag estimate
        if fit[0] > r_max :
            r_max = fit[0]
            tlag_est = (N0-fit[1]) /  fit[0]

        if plot_fig:
            ax.plot(subset["Time"], subset["Time"]*fit[0]+fit[1], color = "red")


    if plot_fig:
        plt.show()

    return [min(time_series["logPopBio"]), max(time_series["logPopBio"]), r_max*math.log(10), max(0,tlag_est)]


def lin_inits(id, plot_fig):
    """
    Calculates initial parameter estimates by rolling regression in linear space for specified ID. Plots a figure displaying the rolling regression if specified
    
    PARAMETERS
    ----------
    id : int
        id number of data for which parameter estiamtes are needed.

    plot_fig : bool
        boolean indiacting whether to plot rolling regression figure

    RETURNS
    -------
    list
        list of intial parameter estimates [Nmin, Nmax, mu_max, t_lag]
    """

    # get time series for this ID
    time_series = data[data["ID"] == id].sort_values("Time")
    time_series.index = range(0,time_series.shape[0]) 

    # get length of time series to calculate number of rolling regression lines
    n_rolls = time_series.shape[0]

    # initialise values
    r_max = 0
    tlag_est = 0

    # plot time series
    if plot_fig:
        fig, ax = plt.subplots()
        ax.scatter(time_series["Time"], time_series["PopBio"])

    # specify number of regressions
    for i in range(0,n_rolls-4):
        # calculate appropriate subset
        subset = time_series.iloc[i:i+5,:]
        # fit linear regression
        fit = stats.linregress(subset["Time"],subset["PopBio"])
        # if new max gradient found
        if fit[0] > r_max :
            # get best estimate for tlag and r_max from this line
            P = np.mean(subset["PopBio"])
            K = np.max(time_series["PopBio"])
            r_max = 1.5*fit[0] / (P * (1 - P / K)) 
            tlag_est = -fit[1] /  fit[0]

        if plot_fig:
            ax.plot(subset["Time"], subset["Time"]*fit[0]+fit[1], color = "red")


    if plot_fig:
        plt.show()
    return [min(time_series["PopBio"]), max(time_series["PopBio"]), r_max, max(0,tlag_est)]



def main():

    # inital value data frame
    inits = pd.DataFrame(columns = ["ID", "loglin", "param", "value"], index = range(len(IDs)*2))

    # calculate log and linear inital values for each ID
    i = 0
    params = ["N0", "Nmax", "rmax", "tlag"]
    for id in IDs:
        init_val_id = {"log" : log_inits(id,0), "lin" : lin_inits(id,0)}
        for loglin in init_val_id.keys():
            for j in range(len(params)):
                inits.loc[i] = [id, loglin, params[j], init_val_id[loglin][j]]
                i+=1

    # save initial values
    inits.to_csv("../data/inits.csv")

    return 0

if __name__ == "__main__":
    main()

