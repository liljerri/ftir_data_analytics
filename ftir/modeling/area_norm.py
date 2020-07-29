def area_norm(df):
    """Normalize to area of 1"""
    from scipy import integrate
    import pandas as pd
    pd.options.mode.chained_assignment = None  # default='warn'
    proteins = list(df.columns)[1:]

    for i in proteins:
        dy = integrate.simps(df.loc[:,i], df.iloc[:,0]) 
        df.loc[:,i] = df.loc[:,i]/abs(dy) #absolute value because decreasing x-values give negative area
    return df.copy()
