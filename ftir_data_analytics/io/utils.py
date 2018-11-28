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


def create_df_from_multiple_files(data_filenames, folder_path):
    """Creates a DataFrame with protein formulations for the given input data
    files
    """

    df = pd.read_csv(folder_path + '/' + data_filenames[0], header=None)
    df.columns = ['freq', 'buffer']
    file_names = ['freq', 'buffer']
    for f in data_filenames[1:]:
        d = pd.read_csv(folder_path + '/' + f, usecols=[1], header=None)
        title = f.split('.')[0]
        d.columns = [title]
        df = pd.concat([df, d], axis=1)
        file_names.append(title)

    # Ensures the dataframe is sorted descending wavenumber
    df = df.sort_values(by=['freq'], ascending=False)
    df = df.reset_index(drop=True)

    # Ensures the dataframe is truncated to maximum 1000 - 3999 wavenumber
    # range to start
    df = df[(df[df.columns[0]] < 3999) & (df[df.columns[0]] > 999)]
    df.reset_index(drop=True, inplace=True)

    return df, file_names
