import db
import os
import csv
import logging as log

log.basicConfig(level=log.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

class DBProcess:
    def __init__(self, dataDir):
        self.dataDir = dataDir

    def getTissueInfo(self, tissueKey):
        if tissueKey == "b_cerebellum":
            return(["Brain_Myeloid", "Brain_Non-Myeloid"], "cerebellum", "female")
        elif tissueKey == "b_cerebral":
           return(["Brain_Myeloid", "Brain_Non-Myeloid"], "cerebral cortex", "female")
        elif tissueKey == "b_hippocampus":
            return(["Brain_Myeloid", "Brain_Non-Myeloid"], "hippocampus", "female")
        elif tissueKey == "b_striatum":
            return(["Brain_Myeloid", "Brain_Non-Myeloid"], "striatum", "female")
        elif tissueKey == "fat":
            return(["Fat"], "gonadal fat pad", "male")
        elif tissueKey == "heart":
            return(["Heart"], "", "male")
        elif tissueKey == "intestine":
            return(["Large_Intestine"], "ascending colon", "male")
        elif tissueKey == "lung":
            return(["Lung"], "", "male")
        elif tissueKey == "spleen":
            return(["Spleen"], "", "male")
        elif tissueKey == "thymus":
            return(["Thymus"], "thymus", "male")
        elif tissueKey == "trachea":
            return("Trachea", "trachea", "male")
        else:
            raise ValueError(f"Unknown tissueKey: {tissueKey}")

    def getOriginalTableName(self, tissueKey):
        return f'tm_temp_{tissueKey}'

    def getIndividuals(self, sex):
        if sex == 'male':
            return ["3_8_M", "3_9_M", "3_10_M", "3_11_M"]
        else:
            return ["3_38_F", "3_39_F", "3_56_F"]
    
    def runOriginal(self, tissueKey):
        log.info(f"Processing original: {tissueKey}")
        originalTableName = self.getOriginalTableName(tissueKey)

        tissues, organismPart, sex = self.getTissueInfo(tissueKey)
        tissueCause = ", ".join(f"'{v.lower()}'" for v in tissues)
        organismPartCause = ""
        if organismPart:
            organismPartCause =  f"AND LOWER(grp.organism_part) = '{organismPart.lower()}'"

        individuals = self.getIndividuals(sex);
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
                    -- AND g.gene in ('ENSMUSG00000000049', 'ENSMUSG00000001027')
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
        results = db.sql(query, 'auto')
        # db.commit()

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

        individuals = self.getIndividuals(sex);
        for ind in individuals:
            query = f'UPDATE {indiviudalTableA} SET "{ind}" = 0 WHERE "{ind}" is null;'
            log.info(query)
            results = db.sql(query, 'auto')

        values = ", ".join(
            f"('{v}', \"{v}\"::numeric)"
            for v in individuals
        )
        nullToZeroTableA = f'tm_temp_{tissueKey}_{method}_a'
        query = f'''
            DROP TABLE IF EXISTS {nullToZeroTableA};
            SELECT
                gene,
                individual,
                SUM(count_sum) AS count_sum,
                individual AS bioreplicate_name
            INTO {nullToZeroTableA}
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

        indiviudalTableB = f'tm_temp_{tissueKey}_{method}_ind_b'
        if sex == 'male':
            selectFields = '"3_8_M" + "3_9_M" AS "3_8_M_n_3_9_M", "3_10_M" + "3_11_M" AS "3_10_M_n_3_11_M"'
            joinFields = '(\'3_8_M_n_3_9_M\', "3_8_M_n_3_9_M"::numeric),(\'3_10_M_n_3_11_M\', "3_10_M_n_3_11_M"::numeric)'
        else:
            selectFields = '"3_38_F" + "3_56_F" AS "3_38_F_n_3_56_F", "3_39_F" AS "3_39_F"'
            joinFields = '(\'3_38_F_n_3_56_F\', "3_38_F_n_3_56_F"::numeric),(\'3_39_F\', "3_39_F"::numeric)'
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

        nullToZeroTableB = f'tm_temp_{tissueKey}_{method}_b'
        query = f'''
            DROP TABLE IF EXISTS {nullToZeroTableB};
            SELECT
                gene,
                individual,
                SUM(count_sum) AS count_sum,
                individual AS bioreplicate_name
            INTO {nullToZeroTableB}
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

        query = f'''
            DROP TABLE IF EXISTS {indiviudalTableA};
            DROP TABLE IF EXISTS {indiviudalTableB};
        '''
        log.info(query)
        db.sql(query, 'auto')
        db.commit()     

if __name__ == "__main__":

    process = DBProcess("/data/loads/liangh/rnaseqload/original")
    
    tissueKeys = ["b_cerebellum", "b_cerebral", "b_hippocampus", "b_striatum", 
               "fat", "heart", "intestine", "lung", "spleen", "thymus", "trachea"]
    for tissueKey in tissueKeys:
        process.runOriginal(tissueKey)
        process.runNullToZero(tissueKey, "null_to_zero")
        break