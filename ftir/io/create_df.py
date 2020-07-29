import pandas as pd
import sys, os, math, warnings, glob, shutil
from pathlib import Path

def excel_df(filetype, data_path, file_name, max_pk_normalize=True, file_delete=False):
    datafile = data_path+'/'+file_name
    empty_meta = {1: {"Message": "No metadata for xlsx files!", 'sequence_line_or_injection':"none", 'sample':'none', 'operator':'none', 'date':'none', 'method':'none', 'detector':'none'}}
    if filetype == 'excel':
        df = pd.read_excel(datafile, header=0)
    else: df = pd.read_csv(datafile, header=0)

    if file_delete: os.remove(datafile)
    filenames = df.columns.to_list()[1:]
    for name in filenames:
        if max_pk_normalize: 
            df[name] = df[name]/df[name].max()
        else:
            df[name] = df[name]
    return df, filenames, empty_meta
