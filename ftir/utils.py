def crop_amide_one(df):
    """ Returns a DataFrame with only the data in the frequency range
    [1705-1600]cm-1
    """
    df = df[(df[df.columns[0]] < 1706) & (df[df.columns[0]] > 1599)]
    df.reset_index(drop=True, inplace=True)
    return df
