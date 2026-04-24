import os
import gzip
import ftplib
import pandas as pd
import random
import logging as log
import psutil
import csv
from decimal import Decimal
from .pseudobulkConfig import PseudobulkConfig

log.basicConfig(level=log.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

class PseudobulkExpt:

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
        self.config = None

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

    def calcTPMAveSD(self, expID, scRNAFile, multiMarkerEnsemblDict, multiEnsemblMarkerDict,
                         ensemblMarkerSet, ensemblSequenceSet, ensemblMarkerDict):
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
                sampleKey = self.config.findSampleKey(colName)
                if sampleKey:
                    indexSampleKeyDict[idx] = sampleKey
                    indexSampleNameDict[idx] = colName

            for i, row in enumerate(reader):
                # if i >= 100:  
                #     break
                gene_id = row[0]
                if not self.isValidGene(gene_id, multiMarkerEnsemblDict, multiEnsemblMarkerDict,
                            ensemblMarkerSet, ensemblSequenceSet):
                        #log.info(f"Skip {gene_id}")
                        continue
                markerKey = ensemblMarkerDict[gene_id]
                #log.info(f"aveTPMDict {gene_id}")
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

    def toPivotSql(self, pivotFields):
        pivotSql = ''
        for i, item in enumerate(pivotFields):
            if i > 0:
                pivotSql += ",\n"
            if isinstance(item, list):
                nameAs = ''
                for i, name in enumerate(item):
                    if i > 0:
                        nameAs += PseudobulkConfig.AND
                    nameAs += name
            else:
                nameAs = item
            pivotSql += f'COALESCE(MAX(tpm_round) FILTER (WHERE bioreplicate_name = \'{nameAs}\'), 0)  AS "{nameAs}"'
        return pivotSql

    def toCaseSql(self, pivotFields):
        caseSql = "CASE "
        for item in pivotFields:
            if isinstance(item, list):
                nameAs = ''
                caseSql += f"WHEN individual IN ("
                for i, name in enumerate(item):
                    if i > 0:
                        nameAs += PseudobulkConfig.AND
                        caseSql += ","
                    nameAs += name
                    caseSql += f"'{name}'"
                caseSql += f") THEN '{nameAs}' "

            else:
                caseSql += f"WHEN individual IN ('{item}') THEN '{item}' "
        caseSql += " END"
        return caseSql   
    
    def createPseudobulkFile(self, db):

        tissue = self.config.bulkData["tissue"]
        organismPart = self.config.bulkData["organismPart"]

        pseudobulkTableName = self.config.getPseudobulkTableName()
        rpkSumTable = f'{pseudobulkTableName}_rpk_sum'
        outputFile = self.config.getPseudobulkDataFile()

        if os.path.exists(outputFile):
            log.info(f"Existing, Skipping creating file: {outputFile}")
            # return outputFile

        pivotSql = ''
        caseSql = ''
        if self.config.runOption["optionName"] == 'A':
            caseSql = 'individual'
            pivotFields = self.config.bulkData["pivotAFields"]
            pivotSql = self.toPivotSql(pivotFields)          
        else:
            pivotFields = self.config.bulkData["pivotBFields"]
            pivotSql = self.toPivotSql(pivotFields)
            caseSql = self.toCaseSql(pivotFields)

        query = '''
            DROP TABLE IF EXISTS {rpkSumTable};
            DROP TABLE IF EXISTS {pseudobulkTableName};

            WITH data AS (
                SELECT g.gene, grp.individual, SUM(g.value) AS count_sum,
                    {caseSql} AS bioreplicate_name
                FROM tm_gene_expression g
                JOIN tm_group grp 
                    ON g.group_id = grp.group_id
                    AND grp.cell_type is not null
                    AND LOWER(grp.tissue) = {tissue}
                    AND LOWER(grp.organism_part) = {organismPart}    
                GROUP BY g.gene, grp.individual
            )
            SELECT d.gene, d.bioreplicate_name, SUM(d.count_sum)::bigint AS count_sum, g.gene_length, 
            0.0 rpk, 0.0 tpm, 0.0  tpm_round
            INTO {pseudobulkTableName}
            FROM data d JOIN tm_gene_info g on g.gene = d.gene  
            GROUP BY d.gene, bioreplicate_name, g.gene_length
            ORDER BY d.gene;

            UPDATE {pseudobulkTableName} SET rpk = count_sum / (gene_length / 1000.0);

            SELECT distinct bioreplicate_name, sum(rpk) rpk_sum
            into {rpkSumTable}
            FROM {pseudobulkTableName}
            GROUP by bioreplicate_name;

            UPDATE {pseudobulkTableName} t
            SET tpm = (t.rpk / s.rpk_sum) * 1000000.0
            FROM {rpkSumTable} s
            WHERE t.bioreplicate_name = s.bioreplicate_name;

            UPDATE {pseudobulkTableName} t
            SET tpm_round = CASE
                WHEN tpm >= 1 THEN ROUND(tpm)
                WHEN tpm > 0 THEN ROUND(
                    tpm::numeric,
                    -CAST(FLOOR(LOG(NULLIF(tpm, 0))) AS int)
                )
                ELSE 0
            END;           

        '''.format(
                pseudobulkTableName=pseudobulkTableName,
                rpkSumTable=rpkSumTable,
                caseSql=caseSql,
                tissue=self.sql_escape(tissue.lower()),
                organismPart=self.sql_escape(organismPart.lower())
        )
        log.info(query)
        db.sql(query, 'auto')
        db.commit() 

        query = '''
            SELECT 
                gene,
                {pivotSql}
            FROM {pseduobulk_table_name}
            GROUP BY gene
            ORDER BY gene
        '''.format(
                pseduobulk_table_name=pseudobulkTableName,
                pivotSql=pivotSql
        )
        log.info(query)
        results = db.sql(query, 'auto')        
        df = pd.DataFrame([
            r.myDict if hasattr(r, "myDict") else dict(r)
            for r in results
        ])
        df.to_csv(outputFile, index=False)
        log.info(f"Export: {outputFile}")

        return outputFile
    
    def writeQNInputOutputFile(self, key, qnInput, qnOutput, ensemblGeneDict):
        scRNAQNInputFile = self.config.getQNInputFile(key)
        scRNAQNOutputFile = self.config.getQNOutputFile(key)
        customHeader = self.config.findSampleNames(qnInput.columns)

        rowLabel = [f"{ensemblGeneDict[idx]}" for idx in qnInput.index]
        qnInput.to_csv(scRNAQNInputFile, index=rowLabel, header=customHeader, index_label='gene')
        qnOutput.to_csv(scRNAQNOutputFile, index=rowLabel, header=customHeader, index_label='gene')

    def writeMatrixFile(self, key, sampleSet, qnInput, ensemblGeneDict, matrix, pseudobulkMatrix):
        matrixLabel = f'{self.config.getShortLabel()}_{key}'
        columns = ["marker_key"]
        pseudobulkColumns = ["gene"]
        sampleNames = self.config.findSampleNames(sampleSet)
        for sampleName in sampleNames:
            columns.append(sampleName)
            pseudobulkColumns.append(f'{matrixLabel}_{sampleName}_bulk_avg')
            pseudobulkColumns.append(f'{matrixLabel}_{sampleName}_qn_avg')
        columns.append("replicate_group_avg")
        pseudobulkColumns.append(f'{matrixLabel}_replicate_group_avg')
        columns.append("num_of_bioreplicates")

        scRNAMatrixOutputFile = self.config.getReplicateGroupAvgFile(key)
        rowLabel = [f"{ensemblGeneDict[idx]}" for idx in qnInput.index]
        df = pd.DataFrame(matrix, columns=columns, index=rowLabel)
        df.to_csv(scRNAMatrixOutputFile, index_label='gene')

        return pd.DataFrame(pseudobulkMatrix, columns=pseudobulkColumns)
    
    def writeMergedMatrixFile(self, pseudobulkDataframeList):
        merged_df = None

        for df in pseudobulkDataframeList:
            if merged_df is None:
                merged_df = df
            else:
                merged_df = pd.merge(merged_df, df, on="gene", how="outer")

        merged_df = merged_df.sort_values("gene").reset_index(drop=True)
        print(merged_df)

        # write detail file
        output_file = os.path.join(self.config.dataDir, self.config.getAllMergedFile())
        merged_df.to_csv(output_file, index=False)

        # write brief file
        # keep gene + every column containing replicate_group_avg
        cols_to_keep = ["gene"] + [col for col in merged_df.columns if "replicate_group_avg" in col]
        df_brief = merged_df[cols_to_keep]
        output_file = os.path.join(self.config.dataDir, self.config.getAllMergedFileTpmOnly())
        df_brief.to_csv(output_file, index=False) 