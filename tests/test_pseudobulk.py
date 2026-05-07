import unittest
import logging as log
import pandas as pd
import os
import numpy as np
import db

from bin.pseudobulkPlot import PseudobulkPlot
from bin.pseudobulkExpt import PseudobulkExpt
from bin.pseudobulkConfig import PseudobulkConfig

log.basicConfig(level=log.INFO, format="%(asctime)s [%(levelname)s] %(message)s")


class PseudobulkTest(unittest.TestCase):

    def setUp(self):
        log.info('PseudobulkTest')

    def tearDown(self):
        # take down anything we've specifically created for each test method
        self.app = None

    def test_create_ht_sample(self):
        log.info('test_create_ht_sample')
        db.useOneConnection(1)

        for bulkData in PseudobulkConfig.BULK_DATA_LIST:
            bulkShortName = bulkData["shortName"]
            tissue = bulkData["tissue"]
            organismPart = bulkData["organismPart"]
            self.createHTSamples(db, bulkShortName, "A", tissue, organismPart, bulkData["pivotAFields"])
            self.createHTSamples(db, bulkShortName, "B", tissue, organismPart, bulkData["pivotBFields"])
        db.commit()

    def findEmapaKey(self, db, tissue, organismPart):
        query = '''
            SELECT t._term_key
            FROM voc_term t, voc_term_emapa e
            WHERE t.term = '{organismPart}'
                AND t._term_key = e._term_key
        '''.format(organismPart = organismPart)
        log.info(query)
        results = db.sql(query, 'auto')
        for r in results:
            return r["_term_key"]
        
    def findHTSampleKey(self, db, sampleName):
        query = '''
            SELECT _sample_key FROM gxd_htsample WHERE _experiment_key = 99680 AND name = '{sampleName}'
        '''.format(sampleName = sampleName)
        # log.info(query)
        results = db.sql(query, 'auto')
        for r in results:
            return r["_sample_key"]        
    
    def findCreateHTSamples(self, db, tissue, organismPart, sex, sampleName):
        sampleKey = self.findHTSampleKey(db, sampleName)
        if sampleKey:
            log.info(f'Skip create: sampleKey: {sampleKey}')
            return

        emapaKey = self.findEmapaKey(db, tissue, organismPart)
        if not emapaKey: 
            log.info(f'No emapaKey: tissue={tissue}, organismPart={organismPart}')
            if tissue == 'Bladder':
                emapaKey = 18236413
            elif tissue == "Fat":
                emapaKey = 18236394
            else:
                emapaKey = 18236394
        if sex == 'Male':
            sexKey = 315165 
        else:
            sexKey = 315164
        cellTypeTermKey = 99536731

        log.info(f'sampleName: {sampleName}  emapaKey: {emapaKey}')
        query = '''
            INSERT INTO mgd.gxd_htsample(
                _sample_key, _experiment_key, _relevance_key, name, age, agemin, agemax, _organism_key, _sex_key,
                _emapa_key, _stage_key, _genotype_key, _celltype_term_key, _rnaseqtype_key)
            VALUES (
                nextval('gxd_htsample_seq'), 99680, 20475450, '{sampleName}',
                'postnatal week 10', -1, -1, 1, {sexKey},
                {emapaKey}, 27, 1, {cellTypeTermKey}, 114866225
            )
        '''.format(sampleName = sampleName,
                   sexKey = sexKey,
                   emapaKey = emapaKey,
                   cellTypeTermKey = cellTypeTermKey)
        log.info(query)
        results = db.sql(query, 'auto') 

    def createHTSamples(self, db, bulkShortName, option, tissue, organismPart, fields):
        for i, field in enumerate(fields):
            if isinstance(field, list):
                sampleName = ""
                for i, f in enumerate(field):
                    if i > 0:
                        sampleName += PseudobulkConfig.AND
                    sampleName += f
            else:
                sampleName = field

            if "M" in sampleName:
                sex = "Male"
            else:
                sex = "Female"

            log.info(f"{tissue} {organismPart}: {sampleName}")
            self.findCreateHTSamples(db, tissue, organismPart, sex, f'{bulkShortName}_{option}_{sampleName}')


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

