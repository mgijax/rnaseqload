import os
import gzip
import ftplib
import shutil
import pandas as pd
import random
from scipy.io import mmread
from scipy import sparse
import numpy as np
import logging as log
import psutil
import time
import csv
from decimal import Decimal

log.basicConfig(level=log.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

class SCRNASeq:
    def __init__(self, exp_id):
        self.exp_id = exp_id
        self.ftp_host = "ftp.ebi.ac.uk"
        self.ftp_path = f"/pub/databases/microarray/data/atlas/sc_experiments_202509/{exp_id}/"

        self.files = [
            f"{exp_id}.aggregated_counts.mtx.gz",
            f"{exp_id}.aggregated_counts.mtx_rows",
            f"{exp_id}.aggregated_counts.mtx_cols",
            f"{exp_id}.sdrf.txt",
        ]

        self.base_dir = os.getenv("RAW_INPUTDIR")

    @staticmethod
    def is_single_cell(exp_id):
        return exp_id.upper().startswith("E-ENAD")        

    @staticmethod
    def print_memory():
        process = psutil.Process(os.getpid())
        mem = process.memory_info().rss / 1024**2  # MB
        log.info(f"Memory usage: {mem:.2f} MB")    

    def download(self, force=False):
        """
        Download experiment files from EBI FTP.

        Behavior:
        - If file does not exist -> download
        - If file exists but size differs -> re-download
        - If file exists and size matches -> skip
        - force=True -> always re-download
        """

        ftp = ftplib.FTP(self.ftp_host)
        ftp.login("anonymous")
        ftp.cwd(self.ftp_path)

        
        os.makedirs(self.base_dir, exist_ok=True)

        for filename in self.files:

            local_path = os.path.join(self.base_dir, filename)

            # Get remote file size
            try:
                remote_size = ftp.size(filename)
            except:
                log.info(f"[WARNING] Could not determine remote size for {filename}")
                remote_size = None

            # Check local file
            if os.path.exists(local_path) and not force:
                local_size = os.path.getsize(local_path)

                if remote_size and local_size == remote_size:
                    log.info(f"[SKIP] {filename} already downloaded ({local_size} bytes)")
                    continue
                else:
                    log.info(
                        f"[REDOWNLOAD] {filename} size mismatch "
                        f"(local={local_size}, remote={remote_size})"
                    )

            else:
                log.info(f"[DOWNLOAD] {filename}")

            with open(local_path, "wb") as f:
                ftp.retrbinary(f"RETR {filename}", f.write)

        ftp.quit()
        log.info("Download step complete.")

    def isValidGene(self, geneID, multiMarkerEnsemblDict, multiEnsemblMarkerDict,
                         ensemblMarkerSet, ensemblSequenceSet):
        
        # multi marker per ensembl
        if geneID in multiMarkerEnsemblDict:
            return False
        # multi ensembl per marker
        if geneID in multiEnsemblMarkerDict:
            return False
        # not in MGI
        if geneID not in ensemblMarkerSet and geneID not in ensemblSequenceSet:        
            return False
        # assoc only with a sequence
        if geneID not in ensemblMarkerSet and geneID in ensemblSequenceSet:
            return False
        return True

    def createJoinedFile(self, joinedFile, multiMarkerEnsemblDict, multiEnsemblMarkerDict,
                         ensemblMarkerSet, ensemblSequenceSet):
        log.info("create joined file")

        joinedFile = os.path.join(
            os.getenv("INPUTDIR"), f"{self.exp_id}.joined.txt"
        )
        # Check if file already exists
        if os.path.exists(joinedFile):
            log.info(f"{joinedFile} already exists. Skipping.")
            return

        # ---- Load genes ----
        genes_path = os.path.join(
            self.base_dir, f"{self.exp_id}.aggregated_counts.mtx_rows"
        )
        genes = pd.read_csv(genes_path, header=None)[0].str.split().str[0].tolist()

        # ---- Load samples ----
        cols_path = os.path.join(
            self.base_dir, f"{self.exp_id}.aggregated_counts.mtx_cols"
        )
        samples = pd.read_csv(cols_path, header=None)[0].tolist()

        # ---- Load SDRF ----
        sdrf_path = os.path.join(self.base_dir, f"{self.exp_id}.sdrf.txt")
        sdrf = pd.read_csv(sdrf_path, sep="\t")

        sample_col = "Comment[technical replicate group]"
        run_col = "Comment[ENA_RUN]"
        sample_to_run = dict(zip(sdrf[sample_col], sdrf[run_col]))

        # ---- Matrix file ----
        mtx_path = os.path.join(
            self.base_dir, f"{self.exp_id}.aggregated_counts.mtx.gz"
        )

        with gzip.open(mtx_path, "rt") as f, open(joinedFile, "w") as out:

            out.write("geneID\trunID\ttpm\tsampleID\n")

            header_parsed = False
            line_count = 0

            for line in f:

                # Skip comments
                if line.startswith("%"):
                    continue

                # First non-comment line is matrix dimensions
                if not header_parsed:
                    header_parsed = True
                    continue

                # sample 1% data 
                if random.random() > 0.001:
                        continue

                gene_i, sample_j, val = line.strip().split()

                gene_id = genes[int(gene_i) - 1]
                if not self.isValidGene(gene_id, multiMarkerEnsemblDict, multiEnsemblMarkerDict,
                         ensemblMarkerSet, ensemblSequenceSet):
                    log.info(f"Skip {gene_id}")
                    continue

                sample_id = samples[int(sample_j) - 1]
                run_id = sample_to_run.get(sample_id, "NA")

                out.write(f"{gene_id}\t{run_id}\t{val}\t{sample_id}\n")

                line_count += 1
                if line_count % 1_000_000 == 0:
                    log.info(f"processed {line_count:,} rows")

        log.info(f"Processed file written to {joinedFile}")

    def getPseudobulkInfo(self, db, query):
        pseudobulkInfo = {}

        log.info(query)
        results = db.sql(query, 'auto')
        for r in results:
            replicateKey = r['replicateKey']
            sample_keys = r['_sample_keys'] 
            sample_key_names = r['sample_key_names'] 
            num_of_samples = r['num_of_samples'] 
            sample_key_set = {
                int(x) for x in sample_keys.split(',') if x.strip()
            } if sample_keys else set() 
            sample_key_names_set = {
                x for x in sample_key_names.split(',') if x.strip()
            } if sample_key_names else set() 

            if not replicateKey in pseudobulkInfo:
                pseudobulkInfo[replicateKey] = {}
            pseudobulkInfo[replicateKey]['sample_keys'] = sample_key_set
            pseudobulkInfo[replicateKey]['sample_key_names'] = sample_key_names_set
            pseudobulkInfo[replicateKey]['num_of_samples'] = num_of_samples
            log.info(f"Replicate: {replicateKey}: num_of_samples = {num_of_samples} ({list(sample_key_set)[:3]})")
        return pseudobulkInfo

    def toReplicates(self, pseudobulkInfo):
        replisetDict = {}
        for key in pseudobulkInfo:
            row = pseudobulkInfo[key]
            replisetDict[key] = row['sample_keys']
        return replisetDict

    def getBioReplicates(self, db, expID, tissue, organism_part, cell_type):
        query ='''
            SELECT distinct concat_ws('|', tissue, organism_part, cell_type, sex, 
                _organism_key, _sex_key, _celltype_term_key, _emapa_key, _genotype_key) AS replicateKey,
                string_agg(DISTINCT _sample_key::text, ',') _sample_keys, count(*) num_of_samples
            FROM gxd_htsample
            WHERE _experiment_key = {expID}
                AND tissue = {tissue}
                AND organism_part = {organism_part}
                AND cell_type = {cell_type}
            GROUP BY tissue, organism_part, cell_type, sex, _organism_key, _sex_key, _celltype_term_key, _emapa_key, _genotype_key
            '''.format(
                expID=99680,
                tissue=self.sql_escape(tissue),
                organism_part=self.sql_escape(organism_part),
                cell_type=self.sql_escape(cell_type)
            )
        
        log.info(query)
        results = db.sql(query, 'auto')
        replisetDict = {}

        for r in results:
            replicateKey = r['replicateKey']
            sample_keys = r['_sample_keys'] 
            sample_key_set = {
                int(x) for x in sample_keys.split(',') if x.strip()
            } if sample_keys else set()
            num_of_samples = r['num_of_samples'] 
            replisetDict[replicateKey] = sample_key_set
            log.info(f"Replicate: {replicateKey}: num_of_samples = {num_of_samples} ({list(sample_key_set)[:3]})")
        return replisetDict
    
    def sql_escape(self, val):
        if val is None:
            return 'NULL'
        return "'" + str(val).replace("'", "''") + "'"
    
    def find_sample_key(self, pseudobulkInfo, colName):
        for key in pseudobulkInfo:
            row = pseudobulkInfo[key]
            sample_key_names = row['sample_key_names']
            for key_name in sample_key_names:
                parts = key_name.split("|")
                if parts[1].lower() == colName.lower():
                    return int(parts[0])
                
    def find_sample_name_by_keys(self, pseudobulkInfo, sampleKeys):
        sampleNames = []
        for sampleKey in sampleKeys:
            sampleName = None
            for key in pseudobulkInfo:
                row = pseudobulkInfo[key]
                sample_key_names = row['sample_key_names']
                for key_name in sample_key_names:
                    parts = key_name.split("|")
                    if int(parts[0]) == int(sampleKey):
                        sampleName = parts[1]
            sampleNames.append(sampleName)
        return sampleNames          
                
    def find_num_of_samples(self, pseudobulkInfo, key):
        row = pseudobulkInfo[key]
        if row:
            return int(row['num_of_samples'])             


    def calcTPMAveSD(self, expID, scRNAFile, multiMarkerEnsemblDict, multiEnsemblMarkerDict,
                         ensemblMarkerSet, ensemblSequenceSet, ensemblMarkerDict, pseudobulkInfo, replisetDict):
        aveTPMDict = {}
        log.info('calcTPMAveSD: %s' % expID)

        with open(scRNAFile, newline='', encoding='utf-8') as f:
            reader = csv.reader(f)
            header = next(reader)

            indexSampleKeyDict = {}
            indexSampleNameDict = {}
            for idx, colName in enumerate(header):
                if idx == 0:
                    continue
                sampleKey = self.find_sample_key(pseudobulkInfo, colName)
                if sampleKey:
                    indexSampleKeyDict[idx] = sampleKey
                    indexSampleNameDict[idx] = colName

            for i, row in enumerate(reader):
                # if i >= 100:  
                #     break
                gene_id = row[0]
                if not self.isValidGene(gene_id, multiMarkerEnsemblDict, multiEnsemblMarkerDict,
                            ensemblMarkerSet, ensemblSequenceSet):
                        log.info(f"Skip {gene_id}")
                        continue
                markerKey = ensemblMarkerDict[gene_id]
                log.info(f"aveTPMDict {gene_id}")
                if gene_id == 'ENSMUSG00000000028':
                    t = 'dd'

                for idx, value in enumerate(row):
                    if idx == 0:
                        continue

                    if not idx in indexSampleKeyDict:
                        continue

                    sampleKey = indexSampleKeyDict[idx]
                    if sampleKey not in aveTPMDict:
                        aveTPMDict[sampleKey] = {}
                    aveTPMDict[sampleKey][markerKey] = Decimal(value)
                    #log.info(f"{gene_id} {sampleKey} {markerKey}: {value }") 
                
        return aveTPMDict        
    
    def create_pseudobulk_file(self, db, output_dir, pseduobulk_option, tissue, organism_part, cell_type):

        tissue_clean = tissue.lower().replace(" ", "_")
        organism_part_clean = organism_part.lower().replace(" ", "_")

        pseduobulk_table_name = f'tm_pseduobulk_{pseduobulk_option.lower()}__{tissue_clean}__{organism_part_clean}'
        rpk_sum_table = f'{pseduobulk_table_name}_rpk_sum'
        output_file = f'{output_dir}Pseudobulk_Option_{pseduobulk_option}__{tissue_clean}__{organism_part_clean}.csv'

        bioreplicate_name_field = ''
        select_fields = ''
        if pseduobulk_option == 'A':
            bioreplicate_name_field = 'individual'
            select_fields = '''
                COALESCE(MAX(tpm_round) FILTER (WHERE bioreplicate_name = '3_38_F'), 0)  AS "3_38_F",
                COALESCE(MAX(tpm_round) FILTER (WHERE bioreplicate_name = '3_56_F'), 0)  AS "3_56_F",
                COALESCE(MAX(tpm_round) FILTER (WHERE bioreplicate_name = '3_39_F'), 0)  AS "3_39_F",
                COALESCE(MAX(tpm_round) FILTER (WHERE bioreplicate_name = '3_8_M'), 0)   AS "3_8_M",
                COALESCE(MAX(tpm_round) FILTER (WHERE bioreplicate_name = '3_9_M'), 0)   AS "3_9_M",
                COALESCE(MAX(tpm_round) FILTER (WHERE bioreplicate_name = '3_10_M'), 0)  AS "3_10_M",
                COALESCE(MAX(tpm_round) FILTER (WHERE bioreplicate_name = '3_11_M'), 0)  AS "3_11_M"
            '''
        else:
            bioreplicate_name_field = '''
                CASE
                    WHEN individual IN ('3_38_F', '3_56_F') THEN '3_38_F_and_3_56_F'
                    WHEN individual IN ('3_8_M', '3_9_M') THEN '3_8_M_and_3_9_M'
                    ELSE individual
                END
            '''
            select_fields = '''
                COALESCE(MAX(tpm_round) FILTER (WHERE bioreplicate_name = '3_38_F_and_3_56_F'), 0)  AS "3_38_F_and_3_56_F",
                COALESCE(MAX(tpm_round) FILTER (WHERE bioreplicate_name = '3_39_F'), 0)             AS "3_39_F",
                COALESCE(MAX(tpm_round) FILTER (WHERE bioreplicate_name = '3_8_M_and_3_9_M'), 0)    AS "3_8_M_and_3_9_M",
                COALESCE(MAX(tpm_round) FILTER (WHERE bioreplicate_name = '3_10_M_and_3_11_M'), 0)  AS "3_10_M_and_3_11_M"
            '''



        query = '''
            DROP TABLE IF EXISTS {rpk_sum_table};
            DROP TABLE IF EXISTS {pseduobulk_table_name};

            WITH data AS (
                SELECT g.gene, grp.individual, SUM(g.value) AS count_sum,
                    {bioreplicate_name_field} AS bioreplicate_name
                FROM tm_gene_expression g
                JOIN tm_group grp 
                    ON g.group_id = grp.group_id
                    AND grp.cell_type is not null
                    AND LOWER(grp.tissue) = {tissue}
                    AND LOWER(grp.organism_part) = {organism_part}
                    -- AND grp.cell_type = 'bladder cell'         
                GROUP BY g.gene, grp.individual
            )
            SELECT d.gene, d.bioreplicate_name, SUM(d.count_sum)::bigint AS count_sum, g.gene_length, 
            0.0 rpk, 0.0 tpm, 0.0  tpm_round
            INTO {pseduobulk_table_name}
            FROM data d JOIN tm_gene_info g on g.gene = d.gene  
            GROUP BY d.gene, bioreplicate_name, g.gene_length
            ORDER BY d.gene;

            UPDATE {pseduobulk_table_name} SET rpk = count_sum / (gene_length / 1000.0);

            SELECT distinct bioreplicate_name, sum(rpk) rpk_sum
            into {rpk_sum_table}
            FROM {pseduobulk_table_name}
            GROUP by bioreplicate_name;

            UPDATE {pseduobulk_table_name} t
            SET tpm = (t.rpk / s.rpk_sum) * 1000000.0
            FROM {rpk_sum_table} s
            WHERE t.bioreplicate_name = s.bioreplicate_name;

            UPDATE {pseduobulk_table_name} t
            SET tpm_round = CASE
                WHEN tpm >= 1 THEN ROUND(tpm)
                WHEN tpm > 0 THEN ROUND(
                    tpm::numeric,
                    -CAST(FLOOR(LOG(NULLIF(tpm, 0))) AS int)
                )
                ELSE 0
            END;           

        '''.format(
                pseduobulk_table_name=pseduobulk_table_name,
                rpk_sum_table=rpk_sum_table,
                output_file=output_file,
                bioreplicate_name_field=bioreplicate_name_field,
                tissue=self.sql_escape(tissue.lower()),
                organism_part=self.sql_escape(organism_part.lower()),
                cell_type=self.sql_escape(cell_type) 
        )
        log.info(query)
        db.sql(query, 'auto')
        db.commit() 

        query = '''
            SELECT 
                gene,
                {select_fields}
            FROM {pseduobulk_table_name}
            GROUP BY gene
            ORDER BY gene
        '''.format(
                pseduobulk_table_name=pseduobulk_table_name,
                select_fields=select_fields
        )
        log.info(query)
        results = db.sql(query, 'auto')        
        df = pd.DataFrame([
            r.myDict if hasattr(r, "myDict") else dict(r)
            for r in results
        ])
        df.to_csv(output_file, index=False)
        log.info(f"Export: {output_file}")


#
# Main
#

# scrna = SCRNASeq("E-ENAD-15")
# scrna.download()

# scrna.createJoinedFile()