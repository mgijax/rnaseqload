import os
import logging as log
import db
from bin.pseudobulkConfig import PseudobulkConfig

log.basicConfig(level=log.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

class PseudobulkSample:

    def __init__(self, exp_id):
        self.exp_id = exp_id
        db.useOneConnection(1)

    def sqlEscape(self, val):
        if val is None:
            return 'NULL'
        return "'" + str(val).replace("'", "''") + "'" 

    def findExperimentKey(self):
        query = '''
            SELECT _object_key FROM acc_accession
            WHERE LOWER(accid) = LOWER({expId})
        '''.format(expId = self.sqlEscape(self.exp_id))
        results = db.sql(query, 'auto')
        for r in results:
            return int(r["_object_key"])             

    def createExperiment(self, name, description):
        accid = self.exp_id
        prefixpart, numericpart = accid.rsplit('-', 1)
        numericpart = int(numericpart)

        # insert experiment
        query = '''
        INSERT INTO gxd_htexperiment (
            _source_key,
            name,
            description,
            _studytype_key,
            _experimenttype_key,
            _evaluationstate_key,
            _curationstate_key
        )
        VALUES (
            20475431,
            {name},
            {description},
            20475462,
            20475437,
            20225942,
            20475421
        )
        RETURNING _experiment_key
        '''.format(name=self.sqlEscape(name),
                   description=self.sqlEscape(description)
                )
        experimentKey = db.sql(query, "auto")[0]['_experiment_key']

        # insert set
        query = '''
        INSERT INTO mgi_set (
            _set_key,
            _mgitype_key,
            name,
            sequencenum
        )
        SELECT
            COALESCE(MAX(_set_key), 0) + 1,
            42,
            'RNASeq Load Experiments',
            1
        FROM mgi_set
        RETURNING _set_key
        '''
        setKey = db.sql(query, "auto")[0]['_set_key']

        # insert set member
        query = '''
        INSERT INTO mgi_setmember (
            _setmember_key,
            _set_key,
            _object_key,
            label,
            sequencenum
        )
        SELECT
            COALESCE(MAX(_setmember_key), 0) + 1,
            {setKey},
            {experimentKey},
            NULL,
            1
        FROM mgi_setmember
        RETURNING _setmember_key
        '''.format(setKey=setKey,
                   experimentKey=experimentKey
                )
        db.sql(query, "auto")

        # insert accession
        query = '''
        INSERT INTO acc_accession (
            _accession_key,
            accid,
            prefixpart,
            numericpart,
            _logicaldb_key,
            _object_key,
            _mgitype_key,
            private,
            preferred
        )
        SELECT
            COALESCE(MAX(_accession_key), 0) + 1,
            {accid},
            {prefixpart},
            {numericpart},
            189,
            {experimentKey},
            42,
            0,
            1
        FROM acc_accession
        '''.format(
            accid=self.sqlEscape(accid),
            prefixpart=self.sqlEscape(prefixpart),
            numericpart=numericpart,
            experimentKey=experimentKey,
        )
        db.sql(query, "auto")
        
        # insert experiment variable
        query = f'''
        INSERT INTO gxd_htexperimentvariable (
            _experiment_key,
            _term_key
        )
        VALUES (
            {experimentKey},
            20475441
        )
        '''
        db.sql(query, "auto")

        return experimentKey
        
    def insertRNASeqSetTables(self, experimentKey):
        query = '''     
            DELETE FROM gxd_htsample_rnaseqsetmember
            WHERE _rnaseqset_key in (SELECT _rnaseqset_key FROM gxd_htsample_rnaseqset WHERE _experiment_key = {experimentKey})
        '''.format(experimentKey = experimentKey)
        db.sql(query, 'auto')

        query = '''     
            DELETE FROM gxd_htsample_rnaseqset
            WHERE _experiment_key = {experimentKey}
        '''.format(experimentKey = experimentKey)
        db.sql(query, 'auto')

        query = '''     
            SELECT DISTINCT _organism_key, _sex_key, _emapa_key, _stage_key, _genotype_key, age, ARRAY_AGG(_sample_key) AS sampleKeys
            FROM GXD_HTSample WHERE _experiment_key = {experimentKey}
            GROUP BY _organism_key, _sex_key, _emapa_key, _stage_key, _genotype_key, age
        '''.format(experimentKey = experimentKey)
        results = db.sql(query, 'auto')
        for row in results:
            self.insertRNASeqSetTable(experimentKey, row)

    def insertRNASeqSetTable(self, experimentKey, row):
        query = '''     
            INSERT INTO gxd_htsample_rnaseqset(
                _experiment_key, age, _organism_key, _sex_key, _emapa_key, _stage_key, _genotype_key, note)
            VALUES ({experimentKey}, {age}, {organismKey}, {sexKey}, {emapaKey}, {stageKey}, {genotypeKey}, '')
            RETURNING _rnaseqset_key
        '''.format(experimentKey = experimentKey,
                   age = self.sqlEscape(row["age"]),
                   organismKey = row["_organism_key"],
                   sexKey = row["_sex_key"],
                   emapaKey = row["_emapa_key"],
                   stageKey = row["_stage_key"],
                   genotypeKey = row["_genotype_key"])
        rnaseqsetKey = db.sql(query, "auto")[0]['_rnaseqset_key']
        for sampleKey in row["sampleKeys"]:
            self.insertRNASeqSetTableMember(rnaseqsetKey, sampleKey)

    def insertRNASeqSetTableMember(self, rnaseqsetKey, sampleKey):
        query = '''     
            INSERT INTO gxd_htsample_rnaseqsetmember(
	        _rnaseqset_key, _sample_key)
	        VALUES ({rnaseqsetKey}, {sampleKey})
        '''.format(rnaseqsetKey = rnaseqsetKey,
                   sampleKey = sampleKey
                   )
        db.sql(query, 'auto')

    def createExperimentHTSamples(self, experimentKey):
        log.info('test_create_ht_sample')
      
        for bulkData in PseudobulkConfig.BULK_DATA_LIST:
            bulkShortName = bulkData["shortName"]
            tissue = bulkData["tissue"]
            organismPart = bulkData["organismPart"]
            self.createPseudobulkSamples(experimentKey, bulkShortName, "A", tissue, organismPart, bulkData["pivotAFields"])
            self.createPseudobulkSamples(experimentKey, bulkShortName, "B", tissue, organismPart, bulkData["pivotBFields"])
            #break        

    def showHTSamples(self, experimentKey):
        query = '''
            SELECT _sample_key, name FROM gxd_htsample
            WHERE _experiment_key = {experimentKey}
        '''.format(experimentKey = experimentKey)
        results = db.sql(query, 'auto')
        for r in results:
            log.info(f'_sample_key={r["_sample_key"]}, name={r["name"]}')        

    def findEmapaKey(self, tissue, organismPart):
        query = '''
            SELECT t._term_key FROM voc_term t, voc_term_emapa e
            WHERE LOWER(t.term) = LOWER({organismPart}) AND t._term_key = e._term_key
        '''.format(organismPart = self.sqlEscape(organismPart))
        #log.info(query)
        results = db.sql(query, 'auto')
        for r in results:
            return r["_term_key"]
        
    def findCellTypeTermKey(self, tissue, organismPart):
        query = '''
            SELECT t._term_key FROM tm_group g, voc_term t 
            WHERE LOWER(g.tissue) = LOWER({tissue}) AND LOWER(g.organism_part) = LOWER({organismPart}) AND LOWER(t.term) = LOWER(g.cell_type)
        '''.format(tissue = self.sqlEscape(tissue),
                   organismPart = self.sqlEscape(organismPart)
                   )
        #log.info(query)
        results = db.sql(query, 'auto')
        for r in results:
            return r["_term_key"]        
        
    def findHTSampleKey(self, experimentKey, sampleName):
        query = '''
            SELECT _sample_key FROM gxd_htsample WHERE _experiment_key = {experimentKey} AND name = {sampleName}
        '''.format(experimentKey = experimentKey,
                   sampleName = self.sqlEscape(sampleName))
        #log.info(query)
        results = db.sql(query, 'auto')
        for r in results:
            return r["_sample_key"]        
    
    def findCreateHTSamples(self, experimentKey, tissue, organismPart, sex, sampleName):
        sampleKey = self.findHTSampleKey(experimentKey, sampleName)
        if sampleKey:
            log.info(f'Skip create: sampleKey: {sampleKey}')
            return

        emapaKey = self.findEmapaKey(tissue, organismPart)
        if not emapaKey: 
            log.info(f'No emapaKey: tissue={tissue}, organismPart={organismPart}')
            if tissue == 'Bladder':
                emapaKey = 36021744  # use urinary bladder muscle
            elif tissue == "Fat" and organismPart == 'interscapular brown adipose tissue':
                emapaKey = 18239620  # use 'brown fat'
            elif tissue == "Heart" and organismPart == '':
                emapaKey = 36418446  # use 'primary heart field'
            elif tissue == "Large_Intestine" and organismPart == '':
                emapaKey = 18238046  # use 'large intestine'             
            else:
                emapaKey = 18236394
        if sex == 'Male':
            sexKey = 315165 
        else:
            sexKey = 315164
        cellTypeTermKey = self.findCellTypeTermKey(tissue, organismPart)
        if not cellTypeTermKey:
            cellTypeTermKey = 99536731

        #log.info(f'sampleName: {sampleName}  emapaKey: {emapaKey}')
        query = '''
            INSERT INTO gxd_htsample(
                _sample_key, _experiment_key, _relevance_key, name, age, agemin, agemax, _organism_key, _sex_key,
                _emapa_key, _stage_key, _genotype_key, _celltype_term_key, _rnaseqtype_key)
            VALUES (
                nextval('gxd_htsample_seq'), {experimentKey}, 20475450, {sampleName},
                'postnatal week 10', -1, -1, 1, {sexKey},
                {emapaKey}, 27, 1, {cellTypeTermKey}, 114866225
            )
        '''.format(experimentKey=experimentKey,
                   sampleName = self.sqlEscape(sampleName),
                   sexKey = sexKey,
                   emapaKey = emapaKey,
                   cellTypeTermKey = cellTypeTermKey)
        #log.info(query)
        results = db.sql(query, 'auto') 

    def createPseudobulkSamples(self, experimentKey, bulkShortName, option, tissue, organismPart, fields):
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

            #log.info(f"{tissue} {organismPart}: {sampleName}")
            self.findCreateHTSamples(experimentKey, tissue, organismPart, sex, f'{bulkShortName}_{option}_{sampleName}')    


 