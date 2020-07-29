import tempfile
import zipfile as zf
import ipywidgets as widgets
import os
from ftir.io.create_df import excel_df

# make class for file_upload to allow separate instances of file handling (may not be necessary)
class file_select(): 
    def __init__(self):
        self.files = widgets.FileUpload(
            accept= '.csv, .xlsx, xls', 
            multiple=False 
        )
        display(self.files)

def zip_importer(file_binary):
    with tempfile.TemporaryFile() as tmp:
        tmp.write(file_binary.files.data[0])
        unzip_file = zf.ZipFile(tmp, 'r')
        unzip_file.extractall(path='/home/jovyan/ftir_data_analytics/data/')
        unzip_file.close()

def xls_importer(file_binary, file_name):    
    with open('/home/jovyan/ftir_data_analytics/data/'+file_name, 'wb') as f: f.write(file_binary.files.data[0])
        
def xls_importer2(file_binary, file_name):    
    with open('/home/jovyan/ftir_data_analytics/data/'+file_name, 'wb') as f: f.write(file_binary.data[0])

def file_upload(file_binary):
    if file_binary.files.metadata[0]['name'].split('.')[1] == 'xlsx':
        print("Excel file imported")  
        file_name = file_binary.files.metadata[0]['name']
        xls_importer(file_binary, file_name)
    else:
        print("csv file imported")
        file_name = file_binary.files.metadata[0]['name']
        xls_importer(file_binary, file_name)

class graphical_VM_file_importer():
    
    def __init__(self):
        
        from IPython.display import display
        import ipywidgets as widgets

        import matplotlib.pyplot as plt
        from pathlib import Path

        import pandas as pd
        import numpy as np
            
        sel_file = widgets.FileUpload(
            accept= '.csv, .xlsx, xls', 
            multiple=False 
        )
        sel_file_output = widgets.Output()
        self.plot_output = widgets.Output()
        
        display(sel_file, sel_file_output)

        def select_file_eventhandler(sel_file):  
            sel_file_output.clear_output()
            self.plot_output.clear_output()
            
            with sel_file_output:       
                print(sel_file.metadata)
                if sel_file.metadata[0]['name'].split('.')[1] == 'xlsx':
                    print("Excel file imported")  
                    file_name = sel_file.metadata[0]['name']
                    xls_importer2(sel_file, file_name)
                else:
                    print("csv file imported")
                    file_name = sel_file.metadata[0]['name']
                    
                    f = sel_file.data[0]
                    import tempfile
                    TEMPDIR=tempfile.TemporaryFile()
                    TEMPDIR.write(f)
                    TEMPDIR.seek(0)
                    ef = pd.read_csv(TEMPDIR)
                    self.df = pd.DataFrame(ef) 
                    self.filenames = self.df.columns.to_list()[1:]

            with self.plot_output:
                self.df.plot(kind='line',x=self.df.columns[0],y=self.df.columns[1:], figsize=(15,10))
                plt.show()
            
        button = widgets.Button(description="Show Me The Data!")
        output = widgets.Output()

        display(button, output)

        def on_button_clicked(b):
            with output:
                sel_file.observe(select_file_eventhandler(sel_file), names = 'ALL')

        button.on_click(on_button_clicked)
        display(self.plot_output)
        