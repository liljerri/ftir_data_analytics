import re
import os
from tkinter import Tk
from tkinter.filedialog import askopenfilename, askopenfilenames

numbers = re.compile(r'(\d+)')


def select_file():
    Tk().withdraw()

    file_types = [
        ('csv or text files', '*.csv *.txt')
    ]

    filename = str(askopenfilename(
        filetypes=file_types,
        title='Select the datafile containing wavenumber in column 1, buffer '
              'spectra in column 2, and sample spectra in remaining columns.'))
    folder_path = os.path.split(filename)[0]
    filename = filename.split('/')[-1]
    return filename, folder_path


def select_files():
    Tk().withdraw()

    file_types = [
        ('csv or text files', '*.csv *.txt')
    ]

    buffer_file = str(askopenfilename(
        filetypes=file_types, title='Choose the buffer file.'))
    reference_file = str(askopenfilename(
        filetypes=file_types, title='Choose the reference file.'))
    data_files = [str(i) for i in askopenfilenames(
        filetypes=file_types, title='Choose the sample files.')]

    folder_path = os.path.split(buffer_file)[0]

    def numerical_sort(value):
        """ Sort function for sorting file names numerically by the first
        number appearing in the file name
        """
        parts = numbers.split(value)
        parts[1::2] = map(int, parts[1::2])
        return parts

    data_files = sorted(data_files, key=numerical_sort)
    file_names = [buffer_file] + [reference_file] + data_files
    file_names = [i.split('/')[-1] for i in file_names]

    return file_names, folder_path
