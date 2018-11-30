"""
Development of FTIR processing algorithm

To Do:
  - Improve baseline and area normalization
  - Improve gaussian peak fitting (maybe restrict peak widths and movement)
  - Improve file open to enable opening in any directory, not just current
    directory


B.Kendrick modifications to file in April 2018 include:
  - updated to read file names with underscores (or any special character)
  - added dataframe sorting function to ensure FTIR data input gets arranged in
    descending wavenumber
  - included output of individual Gaussian fit curves to csv output file
  - removed redundant curve fitting function call for csv output
  - modified the matplotlib plots to stack the gaussian fit plot and the
    residual plot together
  - added extra gaussian peak at 1615 cm-1 for initial fit due to extra
    non-structural peak creating poor fit
  - added calculation of secondary structure elements (helix, turn, etc.) and
    output to csv

Notes:
 Program will throw the following error if any of the initial guess peaks in
 clist has a zero y-value in the FTIR dataset:
    TypeError: Improper input: N=30 must not exceed M=1

 Program will throw a key error if the initial guess peaks in clist lie
 outside the x-data range

 Sometimes you need to tweak the initial guess peaks (in clist) to get the fit
 to work

 Sometimes you need to tweak the height constant in guess_heights function as
 needed for difficult to fit spectra
"""

import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize


def underlying_gaussian(df, col):
    """ Finds the gaussian curves that make up the FTIR signals

    Returns a tuple of the xdata (freq), ydata (summed gaussian)
    and a list of ydata for each underlying gaussian"""
    # Creates an array of x, y data
    data = np.array(pd.concat([df['freq'], df[col]], axis=1))
    # print(data)

    def errfunc(p, x, y):
        """Does this to ensure that parameters are > 0"""
        if min(p) > 0:
            return gaussian_sum(x, *p) - y
        else:
            return 10000000000000000

    # Deconvoluted amide I band frequencies for proteins in water
    clist = [1694, 1691, 1687, 1676, 1672, 1660, 1656, 1650, 1642, 1638, 1634,
             1627, 1621]
    # Gets list of relevant heights corresponding to clist freqs
    hlist = guess_heights(df, col, clist)

    # creates a list of peak widths to use in fitting
    wlist = [5 for i in clist]

    # creates a tuple, e.g. [(0.012, 1700, 5), (0.032, 1688, 5), ...]
    tmp = list(zip(hlist, clist, wlist))

    # print('zipped list of hlist, clist, wlist')
    # print(tmp)

    # unpacks the tmp nested list/tuple into a 1-D array
    guess = np.array([item for sublist in tmp for item in sublist])
    # print('guess np.array from zipped list')
    # print(guess)
    optim, cov, infodict, mesg, ier = optimize.leastsq(
        errfunc, guess, args=(data[:, 0], data[:, 1]), full_output=True)
    xdata = data[:, 0]
    ydata = gaussian_sum(data[:, 0], *optim)
    gausslist = gaussian_list(data[:, 0], *optim)
    resid = infodict['fvec']
    ss_err = (resid**2).sum()
    ss_tot = ((data[:, 1] - data[:, 1].mean())**2).sum()
    rsquared = 1 - (ss_err/ss_tot)
    optim = list(optim)
    heights = optim[0::3]
    centers = optim[1::3]
    widths = optim[2::3]
    areas = []
    for a, b, c in zip(heights, centers, widths):
        area = gaussian_integral(a, c)
        areas.append(area)
    return xdata, ydata, gausslist, resid, rsquared, centers, areas


def guess_heights(df, col, center_list):
    """ Determines guesses for the heights
    """
    heights = []
    freq_map = {}
    for i in df.freq:
        # returns an integer of the frequency with the decimal values truncated
        j = math.floor(i)
        # creates library of floor freq vs intensity
        freq_map[j] = float(df[col].get(df.freq == i))
    for i in center_list:
        # gives intensity value for each frequency in underlying_gaussian clist
        height = freq_map[i]
        # Initial guess is 0.95*actual height at x=freq*
        # The constant 0.95 is arbitrary but seemed to give the best fit
        # Mess with the constant as needed for difficult to fit spectra
        # fills the empty heights list with 0.95*corresponding freqs.
        heights.append(0.95*height)
    return heights


def gaussian(x, height, center, width):
    """ Function defining a gaussian distribution
    """
    return height*np.exp(-(x - center)**2/(2*width**2))


def gaussian_sum(x, *args):
    """ Returns the sum of the gaussian function inputs
    """
    if len(args) % 3 != 0:
        raise ValueError('Args must divisible by 3')
    gausslist = []
    count = 0
    for i in range(int(len(args)/3)):
        gausstemp = gaussian(x, args[count], args[count+1], args[count+2])
        gausslist.append(gausstemp)
        count += 3
    return sum(gausslist)


def gaussian_list(x, *args):
    """ Returns the sum of the gaussian function inputs
    """
    if len(args) % 3 != 0:
        raise ValueError('Args must divisible by 3')
    gausslist = []
    count = 0
    for i in range(int(len(args)/3)):
        gausstemp = gaussian(x, args[count], args[count+1], args[count+2])
        gausslist.append(gausstemp)
        count += 3
    return gausslist


def gaussian_integral(height, width):
    """ Returns the integral of a gaussian curve with the given height, width
    and center
    """
    return height*width*math.sqrt(2*math.pi)


def gaussian_fit(ftir_df, d, folder_path):
    # Calls gaussian peak fitting function
    fit_data = underlying_gaussian(ftir_df, d)
    xdata, ydata, gausslist, resid, rsquared, centers, areas = fit_data

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(211)
    freqs = list(ftir_df['freq'])
    signal = list(ftir_df[d])
    ax.plot(freqs, signal, label='$2^{nd}$ derivative')
    ax.plot(xdata, ydata, label='Gaussian fit')
    for i in range(len(gausslist)):
        ax.plot(xdata, gausslist[i], ls='--')
        gauss_df = pd.DataFrame(gausslist[i], columns=['gauss_peak' + str(i)])
        ftir_df = pd.concat([ftir_df, gauss_df], axis=1)
        
    ax.set_xlim([1705, 1600])
    ax.legend(loc=1)

    ax = fig.add_subplot(212)

    ax.plot(xdata, resid, label='residuals')
    ax.set_xlim([1705, 1600])
    ax.set_xlabel('Wavenumber ($cm^{-1}$)', fontsize=11)
    ax.set_xlim([1705, 1600])
    ax.legend(loc=2)

    plt.savefig(folder_path + '/' + d+'_GaussianFit and Residuals.png')
    plt.close()
    
    ftir_df[''] = np.nan
    ftir_df['centers'] = np.nan
    ftir_df['areas'] = np.nan
    # adjust based on number of initial peak guesses in clist
    ftir_df['centers'][0:13] = centers
    ftir_df['areas'][0:13] = areas
    
    s = ftir_df[['centers', 'areas']]
    s = s[(s.centers <= 1700) & (s.centers >= 1620)]
    print(s)
    turn = s[(s.centers <= 1700) & (s.centers >= 1666)]
    helix = s[(s.centers < 1666) & (s.centers >= 1650)]
    unordered = s[(s.centers < 1650) & (s.centers >= 1644)]
    sheet = s[(s.centers < 1644) & (s.centers >= 1620)]
    print(turn)
    print(helix)
    print(unordered)
    print(sheet)
    
    struct_area = s['areas'].sum()

    pct_turn = round((((turn['areas'].sum())/struct_area)*100), 4)
    pct_helix = round((((helix['areas'].sum())/struct_area)*100), 4)
    pct_unordered = round((((unordered['areas'].sum())/struct_area)*100), 4)
    pct_sheet = round((((sheet['areas'].sum())/struct_area)*100), 4)
    
    ftir_df[''] = np.nan
    ftir_df['% Turn'] = np.nan
    ftir_df['% Turn'][:1] = pct_turn
    ftir_df['% Helix'] = np.nan
    ftir_df['% Helix'][:1] = pct_helix
    ftir_df['% Unordered'] = np.nan
    ftir_df['% Unordered'][:1] = pct_unordered
    ftir_df['% Sheet'] = np.nan
    ftir_df['% Sheet'][:1] = pct_sheet
    
    # ftir_df[['% Turn', '% Helix', '% Unordered', '% Sheet']].
    # plot(kind = 'bar', use_index=False)
    
    ftir_df.to_csv(folder_path + '/' + 'Gaussian fits for '+d+'.csv')
    
    return pct_turn, pct_helix, pct_unordered, pct_sheet


def create_plots(raw_df, folder_path):
    """Creates the following plots for each protein sample:"""
    
    d = raw_df.columns[1]
    print(d)
    structs = gaussian_fit(raw_df, d, folder_path)
    pct_turn, pct_helix, pct_unordered, pct_sheet = structs
    return d, pct_turn, pct_helix, pct_unordered, pct_sheet
