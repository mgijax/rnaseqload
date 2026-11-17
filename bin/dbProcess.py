import db
import os
import csv
import logging as log
from bin.pseudobulkConfig import PseudobulkConfig


log.basicConfig(level=log.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

class DBProcess:
    def __init__(self, dataDir):
        self.dataDir = dataDir

    def getTissueInfo(self, tissueKey):
        return PseudobulkConfig.getTissueInfo(tissueKey)

    def getOriginalTableName(self, tissueKey):
        return f'tm_temp_{tissueKey}'

    def getOptionAIndividuals(self, tissueKey):
        tissueInfo = PseudobulkConfig.findBulk(tissueKey)
        return tissueInfo["pivotAFields"]
        # if sex == 'male':
        #     return ["3_8_M", "3_9_M", "3_10_M", "3_11_M"]
        # else:
        #     return ["3_38_F", "3_39_F", "3_56_F"]

    def getOptionBIndividuals(self, tissueKey):
        tissueInfo = PseudobulkConfig.findBulk(tissueKey)
        return tissueInfo["pivotBFields"]        

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

    def runOriginal(self, tissueKey):
        log.info(f"Processing original: {tissueKey}")
        originalTableName = self.getOriginalTableName(tissueKey)

        tissues, organismPart, sex = self.getTissueInfo(tissueKey)
        tissueCause = ", ".join(f"'{v.lower()}'" for v in tissues)
        organismPartCause = ""
        if organismPart:
            organismPartCause =  f"AND LOWER(grp.organism_part) = '{organismPart.lower()}'"

        individuals = self.getOptionAIndividuals(tissueKey)
        groupFields = ", ".join(
            f"MAX(count_sum) FILTER (WHERE bioreplicate_name = '{v}') AS \"{v}\""
            for v in individuals
        )

        query = '''
            DROP TABLE IF EXISTS {originalTableName};
            with data as (
                SELECT g.gene, grp.cell_type, grp.individual, SUM(g.value) AS count_sum,
                            individual AS bioreplicate_name
                FROM tm_gene_expression g
                JOIN tm_group grp 
                    ON g.group_id = grp.group_id
                    AND grp.cell_type is not null and length(grp.cell_type) > 0
                    AND LOWER(grp.tissue) in ({tissueCause})
                    {organismPartCause}
                    AND grp.sex = '{sex}'
                Group by g.gene, grp.cell_type, grp.individual, grp.cell_type
                ORDER By grp.cell_type, length(grp.individual), grp.individual, grp.cell_type
            )
            SELECT gene, cell_type, {groupFields}
            into {originalTableName}
            FROM data 
            GROUP BY gene, cell_type
            ORDER BY gene, cell_type    
            '''.format(tissueCause=tissueCause, organismPartCause=organismPartCause,
                       sex=sex, originalTableName=originalTableName, groupFields=groupFields
                       )
        log.info(query)
        db.sql(query, 'auto')

    def exportOrignalCSVFile(self, tissueKey):
        log.info(f"Export csv file: {tissueKey}")
        originalTableName = self.getOriginalTableName(tissueKey)        
        query = f'SELECT * FROM {originalTableName} ORDER BY gene, cell_type'
        log.info(query)
        results = db.sql(query, 'auto')

        file = os.path.join(self.dataDir, tissueKey, f'{tissueKey}_original.csv')
        os.makedirs(os.path.dirname(file), exist_ok=True)

        log.info(f'write: {file}')
        with open(file, "w", newline="") as f:
            writer = csv.writer(f)
            if results:
                    writer.writerow(results[0].keys())
                    for row in results:
                        writer.writerow(row.myDict.values())

    def toSelectFields(self, individuals):
        result = []

        for field in individuals:
            if isinstance(field, list):
                expression = " + ".join(f'"{x}"' for x in field)
                alias = PseudobulkConfig.AND.join(field)
                result.append(f'{expression} AS "{alias}"')
            else:
                result.append(f'"{field}" AS "{field}"')
        return ", ".join(result)

    def toJoinFields(self, individuals):
        result = []

        for field in individuals:
            if isinstance(field, list):
                alias = PseudobulkConfig.AND.join(field)
            else:
                alias = field

            result.append(
                f"""('{alias}', "{alias}"::numeric)"""
            )

        return ",".join(result)                         

    def updateMethodTable(self, method, sex, indiviudalTableA, individuals, isDropIndividualTable=False):
        values = ", ".join(
            f"('{v}', \"{v}\"::numeric)"
            for v in individuals
        )
        methodTableA = f'tm_temp_{tissueKey}_{method}_a'
        query = f'''
            DROP TABLE IF EXISTS {methodTableA};
            SELECT
                gene,
                individual,
                SUM(count_sum) AS count_sum,
                individual AS bioreplicate_name
            INTO {methodTableA}
            FROM {indiviudalTableA}
            CROSS JOIN LATERAL (
                VALUES
                    {values}
            ) AS v(individual, count_sum)
            GROUP BY gene, individual
            ORDER BY gene, individual;
        '''
        log.info(query)
        db.sql(query, 'auto')

        optionBIndividuals = self.getOptionBIndividuals(tissueKey)
        indiviudalTableB = f'tm_temp_{tissueKey}_{method}_ind_b'

        selectFields = self.toSelectFields(optionBIndividuals)
        joinFields = self.toJoinFields(optionBIndividuals)
       
        query = f'''
            DROP TABLE IF EXISTS {indiviudalTableB};

            SELECT gene, cell_type, 
                {selectFields}
            INTO {indiviudalTableB}
            FROM {indiviudalTableA} 
            ORDER by gene, cell_type;
        '''
        log.info(query)
        db.sql(query, 'auto')

        methodTableB = f'tm_temp_{tissueKey}_{method}_b'
        query = f'''
            DROP TABLE IF EXISTS {methodTableB};
            SELECT
                gene,
                individual,
                SUM(count_sum) AS count_sum,
                individual AS bioreplicate_name
            INTO {methodTableB}
            FROM {indiviudalTableB}
            CROSS JOIN LATERAL (
                VALUES
                    {joinFields}
            ) AS v(individual, count_sum)
            GROUP BY gene, individual
            ORDER BY gene, individual
        '''
        log.info(query)
        db.sql(query, 'auto')

        dropSql = ""
        if isDropIndividualTable:
            dropSql = f'''
                DROP TABLE IF EXISTS {indiviudalTableA};
                DROP TABLE IF EXISTS {indiviudalTableB};
            '''
        query = f'''
            {dropSql}
            CREATE INDEX idx_{methodTableA}_gene_bioreplicate_name ON {methodTableA} (gene, bioreplicate_name);
            CREATE INDEX idx_{methodTableB}_gene_bioreplicate_name ON {methodTableB} (gene, bioreplicate_name);
        '''
        log.info(query)
        db.sql(query, 'auto')

    def runNullToZero(self, tissueKey, method):
        log.info(f"Processing {method}: {tissueKey}")
        originalTableName = self.getOriginalTableName(tissueKey)
        tissues, organismPart, sex = self.getTissueInfo(tissueKey)

        indiviudalTableA = f'tm_temp_{tissueKey}_{method}_ind_a'
        query = f'''
            DROP TABLE IF EXISTS {indiviudalTableA};
            SELECT * INTO {indiviudalTableA} FROM {originalTableName} order by gene;
        '''
        log.info(query)
        db.sql(query, 'auto')

        individuals = self.getOptionAIndividuals(tissueKey)
        for ind in individuals:
            query = f'UPDATE {indiviudalTableA} SET "{ind}" = 0 WHERE "{ind}" is null;'
            log.info(query)
            db.sql(query, 'auto')

        self.updateMethodTable(method, sex, indiviudalTableA, individuals)

    def runAverage(self, tissueKey, method):
        log.info(f"Processing {method}: {tissueKey}")
        originalTableName = self.getOriginalTableName(tissueKey)
        tissues, organismPart, sex = self.getTissueInfo(tissueKey)
        individuals =self.getOptionAIndividuals(tissueKey)
        values = ", ".join(
            f"(\"{v}\"::numeric)"
            for v in individuals
        )  

        averageTable = f'tm_temp_{tissueKey}_avg'
        query = f'''
            DROP TABLE IF EXISTS {averageTable};
            SELECT gene, cell_type, AVG(value) AS cell_type_avg
            INTO {averageTable}
            FROM {originalTableName}
            CROSS JOIN LATERAL (
                VALUES
                    {values}
            ) AS v(value)
            GROUP BY gene, cell_type
            ORDER BY gene, cell_type;

            CREATE INDEX idx_{averageTable}_gene_cell_type ON {averageTable} (gene, cell_type);        
        '''
        log.info(query)
        db.sql(query, 'auto')        

        indiviudalTableA = f'tm_temp_{tissueKey}_{method}_Ind_a'
        query = f'''
            DROP TABLE IF EXISTS {indiviudalTableA};
            SELECT * INTO {indiviudalTableA} FROM {originalTableName} order by gene;
        '''
        log.info(query)
        db.sql(query, 'auto')
      
        for ind in individuals:
            query = f'''
                UPDATE {indiviudalTableA} i
                SET "{ind}" = (SELECT cell_type_avg FROM {averageTable} a WHERE i.gene = a.gene AND i.cell_type = a.cell_type )
                WHERE "{ind}" is null;
            '''
            log.info(query)
            db.sql(query, 'auto')

        self.updateMethodTable(method, sex, indiviudalTableA, individuals)

        query = f'''
            DROP TABLE IF EXISTS {averageTable};
        '''
        log.info(query)
        db.sql(query, 'auto')

    def removeNull(self, tissueKey, method):
        log.info(f"Processing {method}: {tissueKey}")
        originalTableName = self.getOriginalTableName(tissueKey)
        tissues, organismPart, sex = self.getTissueInfo(tissueKey)

        indiviudalTableA = f'tm_temp_{tissueKey}_{method}_ind_a'
        individuals = self.getOptionAIndividuals(tissueKey)

        selectFields = " + ".join( f'("{v}" IS NOT NULL)::int' for v in individuals)
        selectIndividuals = ",".join( f'"{v}"' for v in individuals )

        query = f'''
            DROP TABLE IF EXISTS {indiviudalTableA};
            with data as (
                SELECT *,
                ({selectFields}) AS non_null_count
                FROM {originalTableName}
            )
            SELECT gene, cell_type, {selectIndividuals}
            INTO {indiviudalTableA}
            FROM data
            WHERE non_null_count = {len(individuals)}
            order by gene, cell_type;
        '''
        log.info(query)
        db.sql(query, 'auto')

        self.updateMethodTable(method, sex, indiviudalTableA, individuals)

    def metabRFCrateTable(self, tissueKey, method):
        #log.info(f"Processing {method}: {tissueKey}")
        indiviudalTableA = f'tm_temp_{tissueKey}_{method}_ind_a'
        tissues, organismPart, sex = self.getTissueInfo(tissueKey)
        individuals = self.getOptionAIndividuals(tissueKey)

        fields = ",".join( f'"{v}" real' for v in individuals )

        query = f'''
            CREATE TABLE IF NOT EXISTS {indiviudalTableA}
            (
                gene text COLLATE pg_catalog."default",
                cell_type character varying COLLATE pg_catalog."default",
                {fields}
            )
        '''
        #log.info(query)
        db.sql(query, 'auto')
        print(f"\copy {indiviudalTableA} FROM 'C:/Users/liangh/Desktop/Imput/{tissueKey}/results.csv' WITH (FORMAT csv,HEADER true);")

    def metabimpute_rf(self, tissueKey, method):
        log.info(f"Processing {method}: {tissueKey}")
        indiviudalTableA = f'tm_temp_{tissueKey}_{method}_ind_a'
        tissues, organismPart, sex = self.getTissueInfo(tissueKey)
        individuals = self.getOptionAIndividuals(tissueKey)
        self.updateMethodTable(method, sex, indiviudalTableA, individuals, isDropIndividualTable=False)      

if __name__ == "__main__":

    process = DBProcess("/data/loads/liangh/rnaseqload/original")
    
    tissueKeys = ["b_cerebellum", "b_cerebral", "b_hippocampus", "b_striatum", 
               "fat", "heart", "intestine", "lung", "spleen", "thymus", "trachea"]
    tissueKeys = ["lung"]
    for tissueKey in tissueKeys:
        # process.runOriginal(tissueKey)
        # process.exportOrignalCSVFile(tissueKey)
        # process.runNullToZero(tissueKey, "null_to_zero")
        # process.runAverage(tissueKey, "average")
        # process.removeNull(tissueKey, "remove_null")
        process.metabimpute_rf(tissueKey, "metabimpute_rf")

        # process.metabRFCrateTable(tissueKey, "metabimpute_rf")
        
        db.commit()
        #break