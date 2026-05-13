import unittest
import logging as log
import pandas as pd
import os
import numpy as np
import db

from bin.pseudobulkPlot import PseudobulkPlot
from bin.pseudobulkSample import PseudobulkSample
from bin.pseudobulkExpt import PseudobulkExpt
from bin.pseudobulkConfig import PseudobulkConfig

log.basicConfig(level=log.INFO, format="%(asctime)s [%(levelname)s] %(message)s")


class PseudobulkTest(unittest.TestCase):

    def setUp(self):
        log.info('PseudobulkTest')

    def tearDown(self):
        # take down anything we've specifically created for each test method
        self.app = None

    def test_create_experiment(self):
        log.info('test_create_experiment')

        name = "Tabula Muris: Transcriptomic characterization of 20 organs and tissues from Mus musculus at single cell resolution"
        description = "Single cell RNA sequencing of single cells across 20 tissues of 3 month aged mice"
        pseudobulkSample = PseudobulkSample("E-ENAD-15")
        
        experimentKey = pseudobulkSample.createExperiment(name, description)
        log.info(f'experimentKey: {experimentKey}')

        pseudobulkSample.createExperimentHTSamples(experimentKey)

        pseudobulkSample.insertRNASeqSetTables(experimentKey)

        pseudobulkSample.showHTSamples(experimentKey)

        log.info(f'experimentKey: {experimentKey}')

        db.commit()

    def test_create_pseduobulked_file(self):
        log.info('test_create_pseduobulked_file')
        db.useOneConnection(1)

        pseudobulkExpt = PseudobulkExpt('')
        pseudobulkExpt.config = PseudobulkConfig(PseudobulkConfig.findBulk("Bladder"), PseudobulkConfig.findOption("A"))

        pseudobulkDataFile = pseudobulkExpt.createPseudobulkFile(db)

        log.info(f'createPseudobulkFile: {pseudobulkDataFile }')

    def test_check_qn(self):
        # check if all columns have same distribution after QN
        log.info('test_check_qn')

        dir = "/data/loads/mgi/rnaseqload/input/"
        false_file_names = ["sex_2_group_QN_input_female_3.csv", "sex_2_group_QN_input_male_3.csv", 
                            "sex_4_group_QN_input_female_1_2.csv", "sex_4_group_QN_input_female_2_1.csv",
                            "sex_4_group_QN_input_male_1_2.csv", "sex_4_group_QN_input_male_2_1.csv"]
        true_file_names = ["sex_2_group_QN_output_female_3.csv", "sex_2_group_QN_output_male_3.csv",
                           "sex_4_group_QN_output_female_1_2.csv", "sex_4_group_QN_output_female_2_1.csv",
                           "sex_4_group_QN_output_male_1_2.csv", "sex_4_group_QN_output_male_2_1.csv"]

        for file_name in true_file_names:
            log.info(f"Check: {file_name}")
            file = os.path.join(dir, file_name)
            is_qn = self.is_quantile_normalized(file)
            if not is_qn:
                self.check_qn_detailed(file)

            log.info(f"  QN check {file_name}: {is_qn}")
            #self.assertEqual(result, True)  
 
    def is_quantile_normalized(self, file_path, tol=1e-6):
        # Load data
        df = pd.read_csv(file_path)

        # If first column is gene IDs, set it as index
        if df.columns[0].lower() in ["gene", "id"]:
            df = df.set_index(df.columns[0])

        # Sort each column
        sorted_cols = np.sort(df.values, axis=0)

        # Compare all columns to the first column
        reference = sorted_cols[:, 0]

        for i in range(1, sorted_cols.shape[1]):
            if not np.allclose(reference, sorted_cols[:, i], atol=tol):
                log.info(f"  Column {df.columns[i]} does NOT match distribution")
                return False

        return True
    
    def check_qn_detailed(self, file_path, tol=1e-6):
        df = pd.read_csv(file_path)

        if df.columns[0].lower() in ["gene", "id"]:
            df = df.set_index(df.columns[0])

        sorted_df = pd.DataFrame(
            np.sort(df.values, axis=0),
            columns=df.columns
        )

        ref = sorted_df.iloc[:, 0]

        for col in sorted_df.columns[1:]:
            diff = np.abs(ref - sorted_df[col])
            max_diff = diff.max()
            log.info(f"  {col}: max difference = {max_diff}")

