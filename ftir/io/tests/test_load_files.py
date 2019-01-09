from unittest import TestCase

from ftir.io.utils import create_df_from_multiple_files, \
    create_df_from_single_file


class BaseDataChecks(TestCase):
    def check_protein_full(self, protein_df, protein_col):
        self.assertEqual(protein_df.columns[protein_col], 'Protein in H2O_0')
        self.assertEqual(999.9028550560, protein_df.index[0])
        self.assertEqual(2137.6900961034, protein_df.index[1180])
        self.assertEqual(3998.6471937487, protein_df.index[3110])
        self.assertEqual(0.2717286944,
                         protein_df.at[999.9028550560, 'Protein in H2O_0'])
        self.assertEqual(0.1141722724,
                         protein_df.at[2137.6900961034, 'Protein in H2O_0'])
        self.assertEqual(0.0036336286,
                         protein_df.at[3998.6471937487, 'Protein in H2O_0'])

    def check_protein_trunc(self, protein_df, protein_col):
        self.assertEqual(protein_df.columns[protein_col], 'Protein in H2O_0')
        self.assertEqual(2000.7699365875, protein_df.index[0])
        self.assertEqual(2499.2750244024, protein_df.index[-1])
        self.assertNotIn(3998.6471937487, protein_df.index)
        self.assertEqual(0.0887036771,
                         protein_df.at[2000.7699365875, 'Protein in H2O_0'])
        self.assertEqual(0.0517123379,
                         protein_df.at[2366.2117707883, 'Protein in H2O_0'])
        self.assertEqual(0.0304686241,
                         protein_df.at[2499.2750244024, 'Protein in H2O_0'])

    def check_buffer_full(self, buffer_df, buffer_col):
        self.assertEqual(buffer_df.columns[buffer_col], 'H2O Buffer_0')
        self.assertEqual(999.9028550560, buffer_df.index[0])
        self.assertEqual(2137.6900961034, buffer_df.index[1180])
        self.assertEqual(3998.6471937487, buffer_df.index[3110])
        self.assertEqual(0.2718908489,
                         buffer_df.at[999.9028550560, 'H2O Buffer_0'])
        self.assertEqual(0.1189842746,
                         buffer_df.at[2137.6900961034, 'H2O Buffer_0'])
        self.assertEqual(0.0086131096,
                         buffer_df.at[3998.6471937487, 'H2O Buffer_0'])


class TestSingleFile(BaseDataChecks):
    def setUp(self):
        """ Run on each test """
        self.buffer_file = './data/H2O Buffer_0.txt'
        self.protein_file = './data/Protein in H2O_0.txt'

    def test_load_single_file(self):
        """ Spot check values upon loading files """
        protein = create_df_from_single_file(self.protein_file, min_freq=0,
                                             max_freq=5000)
        self.check_protein_full(protein, 0)

        buffer = create_df_from_single_file(self.buffer_file, min_freq=0,
                                            max_freq=5000)
        self.check_buffer_full(buffer, 0)

    def test_load_multiple_files(self):
        with self.assertRaises(ValueError):
            protein = create_df_from_single_file(
                [self.protein_file, self.buffer_file])

    def test_frequency_range(self):
        protein = create_df_from_single_file(self.protein_file, min_freq=2000,
                                             max_freq=2500)
        self.check_protein_trunc(protein, 0)

        buffer = create_df_from_single_file(self.buffer_file, min_freq=2000,
                                            max_freq=2000)
        self.assertEqual(buffer.columns[0], 'H2O Buffer_0')
        self.assertTrue(buffer.index.empty)


class TestMultipleFiles(BaseDataChecks):
    def setUp(self):
        self.file = ['./data/H2O Buffer_0.txt']
        self.files = ['./data/H2O Buffer_0.txt',
                      './data/Protein in H2O_0.txt']
        self.different_range = './data/OVA in H2O_0.txt'

    def test_load_single_file(self):
        buffer = create_df_from_multiple_files(
            self.file, min_freq=0, max_freq=5000)
        self.check_buffer_full(buffer, 0)

    def test_load_many_files(self):
        df = create_df_from_multiple_files(
            self.files, min_freq=0, max_freq=5000)
        self.assertEqual(['H2O Buffer_0', 'Protein in H2O_0'],
                         list(df.columns))
        self.check_buffer_full(df, 0)
        self.check_protein_full(df, 1)

    def test_frequency_range(self):
        df = create_df_from_multiple_files(
            self.files, min_freq=2000, max_freq=2000)
        self.assertEqual(['H2O Buffer_0', 'Protein in H2O_0'],
                         list(df.columns))
        self.assertTrue(df.index.empty)

        df = create_df_from_multiple_files(
            self.files, min_freq=2000, max_freq=2500)
        self.assertEqual(['H2O Buffer_0', 'Protein in H2O_0'],
                         list(df.columns))
        self.check_protein_trunc(df, 1)

    def test_different_frequency(self):
        bad_file_list = self.files + [self.different_range]
        with self.assertRaises(ValueError):
            create_df_from_multiple_files(bad_file_list)

