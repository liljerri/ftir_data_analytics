from scipy.spatial import ConvexHull
from scipy import sparse
from scipy.sparse.linalg import spsolve
import numpy as np
import pandas as pd

def find_deriv(df, flip, window_length=5):
    """Adds the 2nd derivative of the chosen signal to the DataFrame
        Window_length=5 is recommended from our studies"""
    from scipy.signal import savgol_filter

    proteins = list(df.columns)[1:]
    drub_df = df[[df.columns[0]]].copy()
    for i in proteins:
        colname = str(i) + '_deriv'

        if flip:
            # negative if want flipped
            dd = -1 * savgol_filter(df[i], deriv=2,
                                    window_length=window_length, polyorder=3
                                    )

            x = df[df.columns[0]].values

            drub = asym_baseline(dd) #uses asymmetric baseline algorithm

            #drub = rubberband(x, dd)
            drub_df[i] = pd.Series(drub, name=colname)
            dd = dd - drub

        else:
            dd = savgol_filter(
                df[i], deriv=2, window_length=window_length, polyorder=3)

        dd = pd.Series(dd, name=colname)

        df = pd.concat([df, dd], axis=1)
        df = df[df.columns.drop(str(i))]
    df.columns = df.columns.str.replace(r"_deriv", "")
    return df

def find_deriv2(deriv_df, flip, window_length=5):
    """Adds the 2nd derivative of the chosen signal to the DataFrame
        Window_length=5 is recommended from our studies"""
    from scipy.signal import savgol_filter

    proteins = list(deriv_df.columns)[1:]
    drub_df = deriv_df[[deriv_df.columns[0]]].copy()
    for i in proteins:
#         colname = str(i) + '_deriv'

        if flip:
            # negative if want flipped
            dd = -1 * savgol_filter(deriv_df[i], deriv=2,
                                    window_length=window_length, polyorder=3
                                    )

            x = deriv_df[deriv_df.columns[0]].values

#             drub = asym_baseline(dd) #uses asymmetric baseline algorithm

            #drub = rubberband(x, dd)
#             drub_df[i] = pd.Series(drub, name=str(i))
#             dd = dd - drub

        else:
            dd = savgol_filter(
                deriv_df[i], deriv=2, window_length=window_length, polyorder=3)

        dd = pd.Series(dd, name=str(i))
        deriv_df = pd.concat([deriv_df.iloc[:,0], dd], axis=1)

#         deriv_df[i] = dd
#         df = df[df.columns.drop(str(i))]

    return deriv_df

def asym_baseline(y):
    lam = 10**2.5
    p = 0.007
    niter = 10
    L = len(y)
    D = sparse.diags([1, -2, 1], [0, -1, -2], shape=(L, L - 2))
    w = np.ones(L)
    for i in range(niter):
        W = sparse.spdiags(w, 0, L, L)
        Z = W + lam * D.dot(D.transpose())
        z = spsolve(Z, w * y)
        w = p * (y > z) + (1 - p) * (y < z)
    return z    

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
    
    def asym_baseline(y):
        lam = 10**2.5
        p = 0.007
        niter = 10
        L = len(y)
        D = sparse.diags([1, -2, 1], [0, -1, -2], shape=(L, L - 2))
        w = np.ones(L)
        for i in range(niter):
            W = sparse.spdiags(w, 0, L, L)
            Z = W + lam * D.dot(D.transpose())
            z = spsolve(Z, w * y)
            w = p * (y > z) + (1 - p) * (y < z)
        return z    

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
        freq = df.columns[0] 

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
        
    elif method == 'asym':
        corrected_spectra = preprocessed_df.apply(asym_baseline)
        
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
