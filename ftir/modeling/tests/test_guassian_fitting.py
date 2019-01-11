from unittest import TestCase
from ftir.io.utils import create_df_from_multiple_files
from ftir.modeling.buffer_subtraction import buffer_subtract

class TestSecondaryStructure(TestCase):

    def setUp(self):
        # Load Data
        self.files = ['./data/H2O Buffer_0.txt',
                      './data/OVA in H2O_0.txt']
        data = create_df_from_multiple_files(self.files)

        # Buffer Subtract
        subtracted = buffer_subtract(data, buffer=0)
        deriviative =

    def test_albumin_control(self):
        pass
