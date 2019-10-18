import pandas as pd
from ftir.modeling.buffer_subtraction import find_buffer_subtraction_constant, buffer_subtract
from ftir.modeling.peak_fitting import gaussian_minimize, gaussian_least_squares, secondary_structure, create_fit_plots, gaussian_list
from ftir.modeling.peak_definitions import yang_h20_2015
from ftir.io.utils import create_df_from_single_file

# get some data
raw_data_filename = "ExampleBSA_IgG1_2ndDer_AmideI.csv"
directory = "/home/jyoung/Devel/ftir_data_analytics/tests/data/"
rawData_df = pd.read_csv(directory + raw_data_filename)
proteins = list(rawData_df.columns)[1:]

structures = {}
structures['File'] = []
structures['Turn'] = []
structures['α-Helix'] = []
structures['Unordered'] = []
structures['β-Sheet'] = []

num_files = len(raw_data_filename)

for i in proteins:
    current_df = rawData_df[['freq', i]].copy()
    area, res = gaussian_least_squares(current_df, current_df.columns[1])

    structs = secondary_structure(area, yang_h20_2015)
    print(i)
    print(structs)
    gaussian_list_data = gaussian_list(rawData_df['freq'], *res.x)
    #plt = create_fit_plots(current_df, i, gaussian_list_data)
    #plt.show()
