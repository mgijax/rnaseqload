import unittest
import logging as log
import pandas as pd
import os
import pandas as pd
import numpy as np
import db

from bin.scPlot import SCPlot
from bin.scRNASeq import SCRNASeq

log.basicConfig(level=log.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

# Test Strains APIs

class SCRNASeqTest(unittest.TestCase):

    def setUp(self):
        log.info('SCRNASeqTest')

    def tearDown(self):
        # take down anything we've specifically created for each test method
        self.app = None

    def test_create_pseduobulked_file(self):
        log.info('test_create_pseduobulked_file')
        db.useOneConnection(1)

        output_dir = '/data/loads/mgi/rnaseqload/input/'
        tissue = "Bladder"
        organism_part = "urinary bladder"
        cell_type = "bladder cell"

        scRNASeq = SCRNASeq('')
        scRNASeq.create_pseudobulk_file(db, output_dir, 'A', tissue, organism_part, cell_type)
        scRNASeq.create_pseudobulk_file(db, output_dir, 'B', tissue, organism_part, cell_type)

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

