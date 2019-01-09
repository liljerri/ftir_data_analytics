"""
coding: utf-8

Algorithm to Buffer subtract, truncate, and area normalize FTIR spectra

5/25/18 Updated directory structure handling

5/29/18 Revised buffer subtraction to use least squares regression in place of
manual input of acceptable minima for iteration.  - replaced the buffer
subtraction-constant iteration from a static list to a modifyable numpy arange
iterator
 - modified the buffer subtraction mean value to obtain better buffer
 subtraction output

To Do:
  - Investigate ways to better subtract the buffer, especially if air bubbles
  are present/slight changes in cell pathlength (peak area based normalization,
  with/without band narrowing, baseline correction, etc.)
"""
import pandas as pd
import numpy as np
from scipy import optimize


WINDOW_SIZE = 12
WINDOW_MIN = 1710
WINDOW_MAX = 2200
MIN_PERIODS = 1


def find_buffer_subtraction_constant(df, sample, buffer, **params):
    """ Returns the constant to use for the buffer signal subtraction

    Parameters
    ----------
    df : Dataframe
        Dataframe containing both the sample and buffer spectra

    sample : String
        Column name in the dataframe for the sample spectra

    buffer : String
        Column name in the dataframe for the buffer spectra

    params : Optional
        Parameters for the buffer subtraction optimization.
            * `window_min` specifies the minimum wavenumber to be used for
            optimization
            * `window_max` specifies the maximum wavenumber to be used for
            optimization
            * `rolling_window` specifies the size of the de-noising window
            * `min_period` specifies the minimum number of observations in the
            de-noising window to have a value
    """

    # Get parameters for optimization function
    window = params.pop('rolling_window', WINDOW_SIZE)
    window_min = params.pop('window_min', WINDOW_MIN)
    window_max = params.pop('window_max', WINDOW_MAX)
    min_periods = params.pop('min_periods', MIN_PERIODS)

    def func(c, dataframe):
        """ Optimization function for buffer subtraction

        Takes a constant `c` and dataframe, and minimizes the maximum range
        across the specified window range. A rolling average is applied to
        reduce the impact of noise.
        """
        sub_series = dataframe[sample] - c * dataframe[buffer]
        trunc = sub_series[(dataframe.index < window_max) &
                           (dataframe.index > window_min)]
        # minimize impact of noise
        smoothed = trunc.rolling(min_periods=min_periods, center=True,
                                 window=window).mean()
        return abs(smoothed.max() - smoothed.min())

    res = optimize.minimize(func, 0.99, args=df)
    return res.x


def buffer_subtract(df, buffer=1, baseline_min=1729, baseline_max=1731,
                    constant=find_buffer_subtraction_constant,
                    constant_params=dict()):
    """ Returns a DataFrame of the subtracted signal data

    Parameters
    ----------
    df : Dataframe
        Dataframe of the absorbance data for the buffer and sample data.

    buffer : Column position (optional)
        Buffer column position. Defaults to zero, but can be specified if
        required.

    baseline_min : Float (optional)
        Wavenumber minimum value used for eliminated any offset. The spectra
        is zeroed to the average of all spectral values between the
        `baseline_min` and `baseline_max` values

    baseline_max : Float (optional)
        Wavenumber maximum value used for eliminated any offset. The spectra
        is zeroed to the average of all spectral values between the
        `baseline_min` and `baseline_max` values

    constant : Callable or Int (optional)
        Can be a callable that takes the dataframe, sample location and buffer
        location and returns an integer. Can also specify an integer constant

    constant_params : Dict (optional)
        Parameters passed to the buffer subtraction `constant` function

    Returns
    -------
    Dataframe
        Returns a dataframe with a the frequency data and all subtracted
        spectra.
    """
    buffer_col = df.columns[buffer]
    samples = df.columns.drop(buffer_col)

    subtracted_spectra = {'freq': df.index}
    for sample in samples:
        if callable(constant):
            scaling_constant = constant(df, sample, buffer_col,
                                        **constant_params)
        else:
            scaling_constant = constant
        offset_result = df[sample] - scaling_constant * df[buffer_col]

        baseline = offset_result[(df.index > baseline_min) &
                                 (df.index < baseline_max)].mean()
        if np.isnan(baseline):
            raise ValueError(
                'Could not determine a baseline value for the buffer '
                'subtraction. Attempted to average the baseline values from '
                '{0} to {1}.'.format(baseline_min, baseline_max))

        result = offset_result - baseline
        subtracted_spectra[sample] = result.values
    final_subtracted = pd.DataFrame(subtracted_spectra)
    return final_subtracted
