from unittest import TestCase
from ftir.io.utils import create_df_from_multiple_files
from ftir.modeling.buffer_subtraction import *
from ftir.modeling.peak_fitting import *
from scipy.signal import savgol_filter


class TestSecondaryStructure(TestCase):

    def setUp(self):
        # Load Data
        self.files = ['./data/H2O Buffer_0.txt',
                      './data/Protein in H2O_0.txt']

        data = create_df_from_multiple_files(self.files)

        file_path = './data/OVA in H2O_0.txt'
        self.ovalbumin = create_df_from_multiple_files([file_path])
        self.ovalbumin_der = savgol_filter(data['OVA in H2O_0'], 9, 2, deriv=2)

        # Buffer Subtract Protein
        # subtracted = buffer_subtract(data, buffer=0)
        # derivative = savgol_filter(
        #     subtracted['Protein in H2O_0'], 9, 2, deriv=2)

    def test_baseline_correction(self):
        sd_baseline_correction(self.ovalbumin_der)

    def test_albumin_control(self):
        file_path = './data/OVA in H2O_0.txt'
        data = create_df_from_multiple_files([file_path])
        deriv = savgol_filter(data['OVA in H2O_0'], 9, 2, deriv=2)
        pass
