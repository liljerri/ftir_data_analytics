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
from ftir.modeling.peak_definitions import yang_h20_2015
from scipy import optimize
from scipy.spatial import ConvexHull


def _split_result_array(res):
    """ Test
    """
    # Pull out meaningful data triplets
    centers = list()
    width = list()
    height = list()
    for i in range(0, len(res.x), 3):
        height.append(res.x[i])
        centers.append(res.x[i+1])
        width.append(res.x[i+2])
    return centers, width, height,


def sd_baseline_correction(df, cols=None, freq=0, flip=False, 
                           method='min', bounds=[1550,1750], inplace=False):
    """ Performs a baseline subtraction on second derivative spectra

    Returns a dataframe of the baseline subtracted data. The dataframe can be
    modified in place, or unmodified with a new dataframe being returned. There
    are two methods for performing the buffer subtraction, a minimum
    value subtraction where zero is set to the smallest value in the defined 
    range, or a rubberband method which uses a convexhull approach to baseline
    correction. 

    Parameters
    ----------
    df : Dataframe
        pandas dataframe containing the FTIR data. The data must be contain a
        column of the wavenumber data, and a column of the spectral data.

    cols : list (default: None)
        List of column names which define the the absorbance data to be fit. 
        These values are column headers, not column indecies. Integer values 
        can be used as column names and are thus ambiguous and not allowed for 
        defining column indecies. 

    freq : Int or Str (default: 0)
        Column index or name for the wavelength data. Defaults to the first 
        column in the dataframe, but can be changed to a different index or 
        a different name if a different column is used for the wavelength 
        column. If an integer is passed, then it is first checked if the 
        integer is a column name. If the integer is not a column name, it is 
        assumed to be the index of the frequency range. 
    
    method :  Str (default: 'min')
        Method used for baseline subtraction. Can be `min` or `rubberband`. 
        `min` subtracts by the minimum value in the defined range. `rubberband`
        applies a convexhull fit of the baseline around the defined range. 
    
    flip : bool (default: False)
        A boolean to flip the data over the x-axis (i.e. muliply by -1)
        
    bounds : iterable of two numbers (default: [1600, 1700])
        Defines the range to use for baseline subtraction. The default values 
        are set around the Amide I band. These values can be expanded to
        include the Amide II or other FTIR features. The max and min value of 
        the interable are used

    Returns
    -------
    Dataframe
        Baseline corrected dataframe across the specified range. 
    """
  
    def minimum(spec):
        """ `minimum` value subtraction function applied to the dataframe """
        return spec - spec.min()

    def rubberband(x, y):
        """ `rubberband` subtraction function """
        v = ConvexHull(np.column_stack([x, y])).vertices
        
        ascending = True if x[0] < x[1] else False
        
        if ascending:
            # rotate vertices until they start from the lowest one
            v = np.roll(v, -v.argmin())
            v = v[:v.argmax() + 1]
        else: 
            # rotate vertices until they start from the highest one
            v = np.roll(v, -v.argmax())
            v = v[:v.argmin()+1]
    
        # Create baseline using linear interpolation between vertices
        return y - np.interp(x, x[v], y[v])

    # get the frequency column name
    if freq not in df.columns and isinstance(freq, int):
        # get column name if an integer and not a column header
        freq = df.columns[freq] 

    # filter around the defined bounds
    if not bounds:
        filtered_df = df
    else:       
        filtered_df = df[(df.iloc[:,0] > min(bounds)) & (df.iloc[:,0] < max(bounds))]
    if len(filtered_df) == 0:
        raise ValueError('Bounds or frequency column definition returned an '
                         'empty frequeny range')
    
    # determine which colums to apply corrections
    if not cols:
        cols = [x for x in df.columns if x != freq]

    # flip over the x-axis if needed   
    if flip:
        preprocessed_df = filtered_df[cols].apply(lambda x: x*-1)
    else:
        preprocessed_df = filtered_df[cols]

    # apply the baseline subtraction method
    if method == 'min':
        corrected_spectra = preprocessed_df.apply(minimum)
        
    elif method == 'rubberband':
        freqCol = filtered_df.iloc[:,0].values
        vals = dict()
        for colName, colData in preprocessed_df.iteritems():
            vals[colName] = rubberband(freqCol, colData.values)
        corrected_spectra = pd.DataFrame(data=vals)

    else:
        raise NameError('name {0} is not a supported baseline method'
                        ''.format(method))
        
    # create the final dataframe with a clean index and return
    filtered_df.reset_index(drop=True, inplace=True)
    corrected_spectra.reset_index(drop=True, inplace=True)
    return pd.concat([filtered_df.iloc[:,0], corrected_spectra], 
                     axis=1, sort=False)



def gaussian_least_squares(df, col, peaks=yang_h20_2015,
                           peak_width=5, params=dict()):

    def fun(p, x, y):
        """ Minimizing across parameter space p, for a given range, x"""
        return gaussian_sum(x, *p) - y

    data = np.array(pd.concat([df.iloc[:,0], df[col]], axis=1))
    heights = guess_heights(df, col, peaks['means'], gain=1.0)
    width = peak_width
    lb = list()
    ub = list()
    guess = list()

    # Make 1-D array for optimization func definition above
    for mean, bound, height in zip(peaks['means'], peaks['uncertainties'],
                                   heights):
        lb.extend([0, bound[0], 0])
        ubh = np.inf if height <= 0 else height
        ub.extend([ubh, bound[1], peak_width*1])
        guess.extend([height*0.95, mean, peak_width])

    args = [fun, np.array(guess)]
    params['args'] = (data[:, 0], data[:, 1])
    params['bounds'] = (np.array(lb), np.array(ub))
    res = optimize.least_squares(*args, **params)

    areas = list()
    for i in range(0, len(res.x), 3):
        height = res.x[i]
        width = res.x[i+2]
        area = gaussian_integral(height, width)
        areas.append(area)
    return areas, res


def gaussian_minimize(
        df, col, peaks=yang_h20_2015, peak_width=5,
        params={'method': 'L-BFGS-B'}):
    """
    Gradient based minimization implementation of the FTIR peak fitting

    Uses the Scipy `optimize.minimize` function for  minimization. Different
    solvers can be specified in the method parameters. This method is likely
    to converge upon a local minimum rather than the global minimum, but
    will converge MUCH faster than the different evolution solution.

    Parameters
    ----------
    df : DataFrame
        pandas dataframe containing the FTIR data. The data must be contain a
        column of the wavenumber data, and a column of the spectral data.

    col : Int or Str
        Column index for the absorbance data to be fit. This value is used to
        reference the `df` column using the standard pandas api, either integer
        or string values are permitted.

    freq : Int or Str (optional)
        Column index or name for the frequency data. Defaults to `freq`.
        Can be changed if a different name is used for the the wavenumber
        (or frequency) column.

    peak_width : Int (optional)
        Maximum peak width. Defaults to 5

    peaks : Peak Definitions (optional)
        A dictionary containing peak definitions to be used. Three kwargs are
        necessary:
            * `peaks` which should provide a list of the the peak means,
            * `uncertainties` which should provide a list of tuple bounds
               around the peak means. These must be ordered the same as the
               peak means list.
            * `assignments` which should provide a list of the peak secondary
               structure assignments. These must be ordered the same as the
               peak means list.
        Defaults to the Yang et. al, Nature Protocol 2015. peak definitions.

    params : Dict (optional)
        A dictionary of kwargs passed to the scipy differential evolution
        optimization algorithm. If `None`, the Default settings within scipy
        are used.


    Returns
    -------
    TBD
    """

    def func(p, x, y):
        """Function to find the local minimum in the boundary space

        Parameters
        ----------
        p : 1-D array
            Inputs are a 1-D array that follow the sequence:
            (peak_height, peak_mean, peak_width)_{n}

        x : 1-D array
            Frequency range for evaluation

        y : 1-D array
            Measured absorbance data corresponding to the frequency range
            provided
        """
        return np.sum((gaussian_sum(x, *p) - y)**2)

    data = np.array(pd.concat([df.iloc[:,0], df[col]], axis=1))
    heights = guess_heights(df, col, peaks['means'], gain=1.0)
    width = peak_width*1
    bounds = list()
    guess = list()

    # Make 1-D array for optimization func definition above
    for mean, bound, height in zip(peaks['means'], peaks['uncertainties'],
                                   heights):
        bounds.append((0, height))
        bounds.append(bound)
        bounds.append((0, width))
        guess.append(height*0.95)
        guess.append(mean)
        guess.append(peak_width)

    args = [func, np.array(guess)]
    params['args'] = (data[:, 0], data[:, 1])
    params['bounds'] = bounds
    res = optimize.minimize(*args, **params)

    areas = list()
    for i in range(0, len(res.x), 3):
        height = res.x[i]
        width = res.x[i+2]
        area = gaussian_integral(height, width)
        areas.append(area)
    return areas, res


def gaussian_differential_evolution(
        df, col, peaks=yang_h20_2015, peak_width=5,
        params=dict()):
    """
    Differential evolution minimization implementation of the FTIR peak fitting

    Uses the Scipy implementation of the Storn and Price differential evolution
    minimization technique. This optimization approach does not use gradient
    methods to find the global minimum across the defined space, and often
    requires a larger number of function evaluations to converge to the local
    minimum than other approaches, e.g. least squared optimization. The
    advantage of this approach is that we can define the bounds of the peak
    positions, and search of the global minima within this defined bounds
    without the worry of converging on a local minimum. The disadvantage is
    that is takes a really long time to run. Like hit go and walk away for a
    couple hours long.

    Parameters
    ----------
    df : DataFrame
        pandas dataframe containing the FTIR data. The data must be contain a
        column of the wavenumber data, and a column of the spectral data.

    col : Int or Str
        Column index for the absorbance data to be fit. This value is used to
        reference the `df` column using the standard pandas api, either integer
        or string values are permitted.

    freq : Int or Str (optional)
        Column index or name for the frequency data. Defaults to `freq`.
        Can be changed if a different name is used for the the wavenumber
        (or frequency) column.

    peak_width : Int (optional)
        Maximum peak width. Defaults to 5

    peaks : Peak Definitions (optional)
        A dictionary containing peak definitions to be used. Three kwargs are
        necessary:
            * `peaks` which should provide a list of the the peak means,
            * `uncertainties` which should provide a list of tuple bounds
               around the peak means. These must be ordered the same as the
               peak means list.
            * `assignments` which should provide a list of the peak secondary
               structure assignments. These must be ordered the same as the
               peak means list.
        Defaults to the Yang et. al, Nature Protocol 2015. peak definitions.

    params : Dict (optional)
        A dictionary of kwargs passed to the scipy differential evolution
        optimization algorithm. If `None`, the Default settings within scipy
        are used.


    Returns
    -------
    TBD
    """

    def func(p, x, y):
        """Function to find the local minimum in the boundary space

        Parameters
        ----------
        p : 1-D array
            Inputs are a 1-D array that follow the sequence:
            (peak_mean, peak_height, peak_width)_{n}

        x : 1-D array
            Frequency range for evaluation

        y : 1-D array
            Measured absorbance data corresponding to the frequency range
            provided
        """
        return np.sum((gaussian_sum(x, *p) - y)**2)

    data = np.array(pd.concat([df.iloc[:,0], df[col]], axis=1))
    heights = guess_heights(df, col, peaks['means'], gain=1.0)
    width = peak_width
    bounds = list()
    # Make 1-D array for optimization to match the func definition above
    for bound, height in zip(peaks['uncertainties'], heights):
        bounds.append((0, height))
        bounds.append(bound)
        bounds.append((0, width))

    args = [func, np.array(bounds)]
    params['args'] = (data[:, 0], data[:, 1])
    res = optimize.differential_evolution(*args, **params)

    areas = []
    centers, width, height = _split_result_array(res)
    for a, b in zip(heights, width):
        area = gaussian_integral(a, b)
        areas.append(area)
    return areas, res


def guess_heights(df, col, center_list, gain=0.95):
    """ Determines guesses for the heights based on measured data.

    Function creates an integer mapping to the measured frequencies, and then
    creates an initial peak height guess of gain*actual height at x=freq*. A
    Default of 0.95 seems to work best for most spectra, but can be change to
    improve convergence.

    Parameters
    ----------
    df : Dataframe
        Dataframe containing the measured absorbance data

    col : string or integer
        Column index for the absorbance data being fit. Accepts either index
        or string convention.

    center_list : iterable of integers
        An iterable of integer peak positions used to find the experiment
        absorbance at a given wavenumber. I.e, the heights are returned at the
        center values in this iterable

    gain : number (optional)
        Fraction of the measured absorbance value to use determine the initial
        guess for the peak height. The value Default value is 0.95, and thus
        by default, all initial peak guesses are 95% of the peak max.

    """
    heights = []
    freq_map = {}
    for i in df.iloc[:,0]:
        j = math.floor(i)
        freq_map[j] = float(df[col].get(df.iloc[:,0] == i))
    for i in center_list:
        height = freq_map[i]
        heights.append(gain*height)
    return heights


def gaussian(x, height, center, width):
    """ Function defining a gaussian distribution
    """
    return height*np.exp(-(x - center)**2/(2*width**2))


def gaussian_sum(x, *args):
    """ Returns the sum of the gaussian function inputs
    """
    return sum(gaussian_list(x, *args))


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


def secondary_structure(areas, peaks):
    """ Returns secondary structure content

    Takes ordered area definitions and peak definitions and returns secondary
    structure quantities.

    Parameters
    ----------

    areas : list
        An ordered list of peak areas. Areas should be ordered from lowest
        wavenumber to the highest wavenumber

    peaks : Dict
        Peak definition dictionary containing a `means` and `assignments`
        items.

    Returns
    -------
    Dict
        Dictionary of secondary structure content

    """

    # check area length
    if not len(areas) == len(peaks['means']):
        raise ValueError('Area definitions do not match the number of peak'
                         'definitions')
    structures = {i: 0 for i in set(peaks['assignments'])}
    for i, assignment in enumerate(peaks['assignments']):
        structures[assignment] += areas[i]/sum(areas)

    return structures


def create_fit_plots(raw_df, col, peak_list):
    """ Creates the following plots for each protein sample:

    Parameters
    ----------
    raw_df : Dataframe
        Pandas dataframe containing the raw FTIR second derivative data.
        Frequency and absorbance values will be used for plotting and residual
        calculations.

    col : String
        Dataframe column name to be used for plotting.

    peak_list : Iterable
        Iterable containing the peak fit data. Each value in the iterable will
        be plotted to show the quality of the model fit. The peaks are also sum
        to show the composite fit and calculate the residuals.

    Returns
    -------
    Plot
        `matplotlib.pyplot` handle is returned. This plot can be displayed via
        `plot.show()` or saved as a file image via `plot.save()`
    """

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(211)

    xdata = raw_df.iloc[:,0]
    y_fit = sum(peak_list)
    ax.plot(xdata, raw_df[col], label='$2^{nd}$ derivative')
    ax.plot(xdata, y_fit, label='Model fit')
    for i in range(len(peak_list)):
        ax.plot(xdata, peak_list[i], ls='--', label='')

    ax.set_xlim([1705, 1600])
    ax.legend(loc=2)

    resid = raw_df[col] - y_fit
    ax = fig.add_subplot(212)
    ax.plot(xdata, resid, label='residuals')
    ax.set_xlim([1705, 1600])
    ax.set_xlabel('Wavenumber ($cm^{-1}$)', fontsize=11)
    ax.set_xlim([1705, 1600])
    ax.legend(loc=2)
    return plt
