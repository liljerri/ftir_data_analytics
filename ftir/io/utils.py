import os
import pandas as pd


def create_df_from_single_file(data_filename, folder_path):
    """ Creates a DataFrame with protein formulations for the given input data
    files
    """

    df = pd.read_csv(folder_path + '/' + data_filename)

    df2 = df.columns.get_values()
    filename_list = df2.tolist()
    file_names = ['freq']
    for f in filename_list[1:]:
        title = f.split('.')[0]
        file_names.append(title)

    df.columns = file_names

    # Ensures the dataframe is sorted descending wavenumber
    df = df.sort_values(by=df.columns[0], ascending=False)
    df = df.reset_index(drop=True)

    # Ensures the dataframe is truncated to maximum 1000 - 3999 wavenumber
    # range to start
    df = df[(df[df.columns[0]] < 3999) & (df[df.columns[0]] > 999)]
    df.reset_index(drop=True, inplace=True)

    return df, file_names


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
    df = pd.read_csv(data_filenames[0], header=None)
    df.columns = ['freq', 'buffer']
    file_names = ['freq', 'buffer']
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
    df = df.sort_values(by=['freq'], ascending=False)
    df = df[(df['freq'] <= max_freq) & (df['freq'] >= min_freq)]
    df.reset_index(drop=True, inplace=True)

    return df, file_names
