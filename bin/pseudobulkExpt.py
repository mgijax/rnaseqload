import os
import gzip
import ftplib
import pandas as pd
import random
import logging as log
import psutil
import csv
import db
from decimal import Decimal
from bin.pseudobulkConfig import PseudobulkConfig
from bin.pseudobulkSample import PseudobulkSample

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
        self.isWriteDetailFiles = False
        self.experimentKey = None

    @staticmethod
    def isPseudobulkData(exp_id):
        return exp_id.upper().startswith("E-ENAD")        

    @staticmethod
    def print_memory():
        process = psutil.Process(os.getpid())
        mem = process.memory_info().rss / 1024**2  # MB
        log.info(f"Memory usage: {mem:.2f} MB") 

    def initExperimentKey(self):
        pseudobulkSample = PseudobulkSample(self.exp_id)
        self.experimentKey = pseudobulkSample.findExperimentKey()
        if self.experimentKey:
            log.info(f'Existing {self.exp_id}, experimentKey: {self.experimentKey}')
            return
        
        name = "Tabula Muris: Transcriptomic characterization of 20 organs and tissues from Mus musculus at single cell resolution"
        description = "Single cell RNA sequencing of single cells across 20 tissues of 3 month aged mice"
        self.experimentKey = pseudobulkSample.createExperiment(name, description)
        pseudobulkSample.createExperimentHTSamples(self.experimentKey)
        pseudobulkSample.insertRNASeqSetTables(self.experimentKey)
        pseudobulkSample.showHTSamples(self.experimentKey)
        log.info(f'{self.exp_id}, experimentKey: {self.experimentKey}')
        db.commit()

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
            return (False, "multi marker per ensembl")
        # multi ensembl per marker
        if geneID in multiEnsemblMarkerDict:
            return (False, f"multi ensembl per marker: {geneID}={multiEnsemblMarkerDict[geneID]}")
        # not in MGI
        if geneID not in ensemblMarkerSet and geneID not in ensemblSequenceSet:        
            return (False, "not in MGI")
        # assoc only with a sequence
        if geneID not in ensemblMarkerSet and geneID in ensemblSequenceSet:
            return (False, "assoc only with a sequence")
        return (True, "")

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
                isValid, errorMessage = self.isValidGene(gene_id, multiMarkerEnsemblDict, multiEnsemblMarkerDict,
                         ensemblMarkerSet, ensemblSequenceSet)
                if not isValid:
                    log.info(f"Skip {gene_id} {errorMessage}")
                    #continue

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
                expID=self.experimentKey,
                tissue=self.sqlEscape(tissue),
                organism_part=self.sqlEscape(organism_part),
                cell_type=self.sqlEscape(cell_type)
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
    
    def sqlEscape(self, val):
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

            isvalidGeneCount = 0
            for i, row in enumerate(reader):
                # if i >= 100:  
                #     break
                gene_id = row[0]
                isValid, errorMessage = self.isValidGene(gene_id, multiMarkerEnsemblDict, multiEnsemblMarkerDict,
                                        ensemblMarkerSet, ensemblSequenceSet)
                if not isValid:
                        isvalidGeneCount += 1
                        log.info(f"Skip {isvalidGeneCount}. {gene_id} {errorMessage}")
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

    def toPivotSql(self, columnNames, isUseCoalesce=True):
        pivotSql = ''
        for i, name in enumerate(columnNames):
            if i > 0:
                pivotSql += ",\n"

            if isUseCoalesce:
                pivotSql += f'COALESCE(MAX(tpm_round) FILTER (WHERE bioreplicate_name = \'{name}\'), 0)  AS "{name}"'
            else:
                pivotSql += f'MAX(tpm_round) FILTER (WHERE bioreplicate_name = \'{name}\')  AS "{name}"'
        return pivotSql

    def toCaseSql(self, pivotFields):
        caseSql = "CASE "
        for item in pivotFields:
            if isinstance(item, list):
                name = ''
                caseSql += f"WHEN individual IN ("
                for i, field in enumerate(item):
                    if i > 0:
                        name += PseudobulkConfig.AND
                        caseSql += ","
                    name += field
                    caseSql += f"'{field}'"
                caseSql += f") THEN '{name}' "

            else:
                caseSql += f"WHEN individual IN ('{item}') THEN '{item}' "
        caseSql += " END"
        return caseSql 

    def toCauseSql(self, columnNames):
        causeSql = ""
        for i, name in enumerate(columnNames):
            if i > 0:
                causeSql += " AND "
            causeSql += f'"{name}">0'
        return causeSql

    def toColumnNames(self, pivotFields):
        columnNames = []
        for field in pivotFields:
            if isinstance(field, list):
                name = ''
                for i, n in enumerate(field):
                    if i > 0:
                        name += PseudobulkConfig.AND
                    name += n
            else:
                name = field
            columnNames.append(name)
        return columnNames        
    
    def createPseudobulkFile(self, db):

        tissue = self.config.bulkData["tissue"]
        organismPart = self.config.bulkData["organismPart"]

        pseudobulkTableName = self.config.getPseudobulkTableName()
        rpkSumTable = f'{pseudobulkTableName}_rpk_sum'
        outputFile = self.config.getPseudobulkDataFile()

        if os.path.exists(outputFile):
            log.info(f"Existing, Skipping creating file: {outputFile}")
            # return outputFile

        caseSql = ''
        causeSql = ''
        if self.config.runOption["optionName"] == 'A':
            pivotFields = self.config.bulkData["pivotAFields"]
            caseSql = 'individual'
            columnNames = pivotFields
        else:
            pivotFields = self.config.bulkData["pivotBFields"]
            caseSql = self.toCaseSql(pivotFields)
            columnNames = self.toColumnNames(pivotFields)

        pivotSql = self.toPivotSql(columnNames)       
        causeSql = self.toCauseSql(columnNames)       

        organismPartCause = ''
        if organismPart and len(organismPart) > 0:
            organismPartCause = f"AND LOWER(grp.organism_part) = {self.sqlEscape(organismPart.lower())}"

        imput = self.config.imput
        dataQuery = '''
                SELECT g.gene, grp.individual, SUM(g.value) AS count_sum,
                    {caseSql} AS bioreplicate_name
                FROM tm_gene_expression g
                JOIN tm_group grp 
                    ON g.group_id = grp.group_id
                    AND grp.cell_type is not null
                    AND LOWER(grp.tissue) = {tissue}
                    {organismPartCause}    
                GROUP BY g.gene, grp.individual
        '''
        if imput:
            dataQuery = f"SELECT * FROM tm_im_lung_{imput}_{self.config.runOption["optionName"].lower()}"

        query = '''
            DROP TABLE IF EXISTS {rpkSumTable};
            DROP TABLE IF EXISTS {pseudobulkTableName};

            WITH data AS (
                {dataQuery}
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
                dataQuery=dataQuery,
                pseudobulkTableName=pseudobulkTableName,
                rpkSumTable=rpkSumTable,
                caseSql=caseSql,
                tissue=self.sqlEscape(tissue.lower()),
                organismPartCause=organismPartCause
        )
        log.info(query)
        db.sql(query, 'auto')
       
        query = '''
            WITH data as (
                SELECT 
                    gene,
                    {pivotSql}
                FROM {pseduobulk_table_name}
                GROUP BY gene
            )
            SELECT * FROM data
            --WHERE {causeSql}
            ORDER BY gene
        '''.format(
                pseduobulk_table_name=pseudobulkTableName,
                pivotSql=pivotSql,
                causeSql=causeSql
        )
        log.info(query)
        results = db.sql(query, 'auto')        
        df = pd.DataFrame([
            r.myDict if hasattr(r, "myDict") else dict(r)
            for r in results
        ])
        df.to_csv(outputFile, index=False)
        log.info(f"Export: {outputFile}")
        # only if want to preserve them for debug
        # db.commit()

        #self.runGeneSummary(db, self.config.runOption["optionName"], tissue, organismPart, pseudobulkTableName, columnNames)
        return outputFile
    
    def toGeneCountSelectField(self, option, columnNames, countColumnName, sampleCountName):
        if option == 'A':
            op = 'OR'
        else:
            op = 'AND'
        caseSql = 'CASE WHEN '
        for i, name in enumerate(columnNames):
            if i > 0:
                caseSql += f' {op} '
            caseSql += f'"{name}" IS NULL'
        caseSql += f'  THEN NULL ELSE '

        for i, name in enumerate(columnNames):
            if i > 0:
                caseSql += ' + '
            if option == 'A':
                caseSql += f'"{name}"'
            else:
                caseSql += f'COALESCE("{name}",0)'
        caseSql += f' END AS {countColumnName},\n'

        for i, name in enumerate(columnNames):
            if i > 0:
                caseSql += ' + '
            caseSql += f'("{name}" IS NOT NULL)::int'
        caseSql += f' AS {sampleCountName}'

        return caseSql    
    
    def runGeneSummary(self, db, option, tissue, organismPart, pseudobulkTableName, columnNames):
        pivotSql = self.toPivotSql(columnNames, False)
        femaleColumnName = f'"{option}_Female"'
        maleColumnName = f'"{option}_Male"'
        femaleSampleCountName = f'"{option}_Female_Sample_Count"'
        maleSampleCountName = f'"{option}_Male_Sample_Count"'

        femaleColumns = []
        maleColumns = []
        for name in columnNames:
            if name.endswith("F"):
                femaleColumns.append(name)
            elif name.endswith("M"):
                maleColumns.append(name)

        query = '''
            WITH summaryData AS (
                WITH geneData AS (
                    WITH data AS (
                        SELECT
                            gene,
                            {pivotSql}
                        FROM {pseduobulk_table_name}
                        GROUP BY gene
                    )
                    SELECT *,
                        {femaleSelectField},
                        {maleSelectField}
                    FROM data
                    ORDER BY gene
                )
                SELECT 
                    COUNT({femaleColumnName}) Female_Count, 
                    COUNT({maleColumnName}) Male_Count,
                    COUNT(*) FILTER (WHERE {femaleColumnName} IS NOT NULL AND {maleColumnName} IS NOT NULL) AS Female_and_Male_Count, 
                    COUNT(*) Total_Genes,
                    COUNT(*) FILTER (WHERE {femaleSampleCountName} < (SELECT MAX({femaleSampleCountName}) FROM geneData)) AS Female_Null_Count,
                    COUNT(*) FILTER (WHERE {femaleSampleCountName} = 0) AS Female_All_Null,
                    COUNT(*) FILTER (WHERE {maleSampleCountName} < (SELECT MAX({maleSampleCountName}) FROM geneData)) AS Male_Null_Count,
                    COUNT(*) FILTER (WHERE {maleSampleCountName} = 0) AS Male_All_Null
                    
                FROM geneData
            )
            SELECT Total_Genes, Female_Null_Count,
                Female_All_Null, (Total_Genes-Female_All_Null) AS Female_Adjusted_Total, 
                (Female_Null_Count-Female_All_Null) AS Female_Recoverable_Genes,
                ROUND(100.0 *  (Female_Null_Count-Female_All_Null)/(Total_Genes-Female_All_Null), 2) AS Female_Recoverable_Percentage,

                Male_Null_Count,
                Male_All_Null, (Total_Genes-Male_All_Null) AS Male_Adjusted_Total, 
                (Male_Null_Count-Male_All_Null) AS Male_Recoverable_Genes,
                ROUND(100.0 *  (Male_Null_Count-Male_All_Null)/(Total_Genes-Male_All_Null), 2) AS Male_Recoverable_Percentage
            FROM summaryData
        '''.format(pseduobulk_table_name=pseudobulkTableName,
                   pivotSql=pivotSql,
                   femaleSelectField=self.toGeneCountSelectField(option, femaleColumns, femaleColumnName, femaleSampleCountName),
                   maleSelectField=self.toGeneCountSelectField(option, maleColumns, maleColumnName, maleSampleCountName),
                   femaleColumnName=femaleColumnName,
                   maleColumnName=maleColumnName,
                   femaleSampleCountName=femaleSampleCountName,
                   maleSampleCountName=maleSampleCountName
                  )
        log.info(query)

        results = db.sql(query, 'auto')
        for r in results:
            self.insertGeneSummary(db, option, tissue, organismPart, r)

    def insertGeneSummary(self, db, option, tissue, organismPart, row):

        if len(organismPart) > 0:
            orgPart = organismPart
        else:
            orgPart = 'All'
        query = '''
            CREATE TABLE IF NOT EXISTS tm_gene_count_summary
            (
                tissue text COLLATE pg_catalog."default",
                organism_part text COLLATE pg_catalog."default",
                bulk_option text COLLATE pg_catalog."default",
                total_genes bigint,
                female_null_count bigint,
                female_all_null bigint,
                female_adjusted_total bigint,
                female_recoverable_genes bigint,
                female_recoverable_percentage numeric,
                male_null_count bigint,
                male_all_null bigint,
                male_adjusted_total bigint,
                male_recoverable_genes bigint,
                male_recoverable_percentage numeric
            );
            INSERT INTO tm_gene_count_summary(
	            tissue, organism_part, bulk_option, total_genes, 
                female_null_count, female_all_null, female_adjusted_total, female_recoverable_genes, female_recoverable_percentage, 
	            male_null_count, male_all_null, male_adjusted_total, male_recoverable_genes, male_recoverable_percentage)
	        VALUES ({tissue}, {organismPart}, {option}, {total_genes}, 
                {female_null_count}, {female_all_null}, {female_adjusted_total}, {female_recoverable_genes}, {female_recoverable_percentage}, 
                {male_null_count}, {male_all_null}, {male_adjusted_total}, {male_recoverable_genes}, {male_recoverable_percentage})
            '''.format(tissue=self.sqlEscape(tissue),
                       organismPart=self.sqlEscape(orgPart),
                       option=self.sqlEscape(option),
                       total_genes=row["total_genes"],
                       female_null_count=row["female_null_count"],
                       female_all_null=row["female_all_null"],
                       female_adjusted_total=row["female_adjusted_total"],
                       female_recoverable_genes=row["female_recoverable_genes"],
                       female_recoverable_percentage=row["female_recoverable_percentage"],
                       male_null_count=row["male_null_count"],
                       male_all_null=row["male_all_null"],
                       male_adjusted_total=row["male_adjusted_total"],
                       male_recoverable_genes=row["male_recoverable_genes"],
                       male_recoverable_percentage=row["male_recoverable_percentage"])
        log.info(query)
        results = db.sql(query, 'auto')
        db.commit()
    
    def writeQNInputOutputFile(self, key, qnInput, qnOutput, ensemblGeneDict):
        if not self.isWriteDetailFiles:
            return
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
        if self.isWriteDetailFiles:  
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
        #if self.isWriteDetailFiles:
        output_file = os.path.join(self.config.dataDir, self.config.getAllMergedFile())
        merged_df.to_csv(output_file, index=False)

        # write brief file
        # keep gene + every column containing replicate_group_avg
        cols_to_keep = ["gene"] + [col for col in merged_df.columns if "replicate_group_avg" in col]
        df_brief = merged_df[cols_to_keep]
        output_file = os.path.join(self.config.dataDir, self.config.getAllMergedFileTpmOnly())
        df_brief.to_csv(output_file, index=False) 