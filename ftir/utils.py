
def crop_spectra(df, min_freq, max_freq):
    """ Crop a FTIR spectrum dataframe based upon max and min frequency values.

    Params
    ------
    df : Dataframe
        Pandas dataframe including spectral data. Frequency data must be the
        dataframe index.

    min_freq : Int
        Integer value for the minimum value to use for truncation. Inclusive,
        meaning that the `min_freq` value will be including in the returned
        dataframe.

    max_freq : Int
        Integer value for the maximum value to use for truncation. Inclusive,
        meaning that the `max_freq` value will be including in the returned
        dataframe.

    Returns
    -------

    Dataframe
        Truncated dataframe. This is a new dataframe, now a slice of the
        original dataframe.

    """
    df = df.truncate(before=min_freq, after=max_freq)
    df.reset_index(drop=True, inplace=True)
    return df


def crop_amide_one(df, min_freq=1599, max_freq=1706):
    """ Simple wrapper function with default values for amide one region.
    """
    return crop_spectra(df, min_freq=min_freq, max_freq=max_freq)
