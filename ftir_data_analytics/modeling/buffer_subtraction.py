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

from scipy import optimize


def fun(c, df):
    
    df1 = df.copy()
    df1['subtr'] = df1[df1.columns[2]] - c*df1[df1.columns[1]]

    # orig 2500 - 1720,new 2200, 1710
    df2 = df1[(df1[df1.columns[0]] < 2200) & (df1[df1.columns[0]] > 1710)]

    # minimize impact of noise
    df3 = df2['subtr'].rolling(min_periods=1, center=True, window=12).mean()
    
    return abs(df3.max()-df3.min())


def get_constant(df):
    """ Returns the constant to use for the buffer signal subtraction
    """

    min_params = optimize.minimize(fun, 0.99, args=df)
    d = min_params.x
    print('buffer subtraction factor = ', d)
    result = df[df.columns[2]] - d*df[df.columns[1]]
    return result


def buffer_subtract(df):
    """ Updates the DataFrame to have subtracted signal data
    """
    result = get_constant(df)
    df[df.columns[2] + '_subtracted'] = result
    df1 = df[(df[df.columns[0]] > 1729) & (df[df.columns[0]] < 1731)]
    baseline_value = df1[df1.columns[3]].values[0]
    df[df.columns[3]] = result - baseline_value
    return df


def crop(df):
    """ Returns a DataFrame with only the data in the frequency range
    [1705-1600]cm-1
    """
    df = df[(df[df.columns[0]] < 1706) & (df[df.columns[0]] > 1599)]
    df.reset_index(drop=True, inplace=True)
    return df
