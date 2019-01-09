import os
import pandas as pd


def create_df_from_single_file(data_filename,  max_freq=3999,
                               min_freq=1000):
    """ Creates a DataFrame with protein formulations for the given input data
    files
    """

    df = pd.read_csv(data_filename, header=None)

    col_names = ['freq']
    title = os.path.basename(data_filename).split('.')[0]
    title = title.split('.')[0]
    col_names.append(title)
    df.columns = col_names

    # Ensures the dataframe is sorted descending wavenumber
    df.set_index(df.columns[0], inplace=True)
    df.sort_index(ascending=True, inplace=True)

    # Ensures the dataframe is truncated to maximum 1000 - 3999 wavenumber
    # range to start
    df = df.truncate(before=min_freq, after=max_freq)

    return df


def create_df_from_multiple_files(
        data_filenames, max_freq=3999, min_freq=1000):
    """Creates a DataFrame with protein formulations for the given input data
    files

    data_file_names : Iterable
        Iterable of file path strings. Must be ordered such that the buffer
        file path is first. All frequency data must be the same to be analyzed
        together. The frequency data of the buffer file is taken and, and all
        other frequency values are checked against the buffer frequency.
    """
    initial_file = data_filenames[0]
    file_names = []
    df = pd.read_csv(initial_file, header=None)
    title = os.path.basename(initial_file).split('.')[0]
    df.columns = ['freq', title]

    for f in data_filenames[1:]:
        d = pd.read_csv(f, header=None)
        title = os.path.basename(f).split('.')[0]
        d.columns = ['freq', title]
        if not df['freq'].equals(d['freq']):
            raise ValueError(
                'Frequency ranges do not match. Sample data file: {0}\n'
                'Buffer data file: {1}'.format(f, data_filenames[0]))
        df = pd.concat([df, d[title]], axis=1)
        file_names.append(title)

    # Ensures the dataframe is sorted descending wavenumber
    df.set_index('freq', inplace=True)
    df.sort_index(ascending=True, inplace=True)
    df = df.truncate(before=min_freq, after=max_freq)
    return df
