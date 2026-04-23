import logging as log
import os
from itertools import count

log.basicConfig(level=log.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

class PseudobulkConfig:
    
    AND = '_n_'   # used to join name for multiple fields

    def __init__(self, bulkData, runOption):
        self.bulkData = bulkData
        self.runOption = runOption

        self.dataDir = os.getenv('INPUTDIR')
        self.initReplicateGroup()

    def initReplicateGroup(self):
        counter = count(1)
        female = {"groupName": "Female",
                  "replicates": []
                  }
        male = {"groupName": "Male",
                  "replicates": []
                }
        replicateGroups = [female, male]

        if self.runOption["optionName"] == "A":
            fields = self.bulkData["pivotAFields"]
        else:
            fields = self.bulkData["pivotBFields"]
        for i, field in enumerate(fields):
            if isinstance(field, list):
                name = ""
                for i, f in enumerate(field):
                    if i > 0:
                        name += PseudobulkConfig.AND
                    name += f
            else:
                name = field
            
            if "F" in name:
                group = female
            else:
                group = male
            group["replicates"].append(
                {
                    "replicateName": name,
                    "sampleKey": next(counter)
                }
            )

        self.runOption["replicateGroups"] = replicateGroups

    def findReplicateGroup(self, groupName):
        for replicateGroup in self.runOption["replicateGroups"]:
            if groupName == replicateGroup["groupName"]:
                return replicateGroup
    
    def findNumberOfReplicates(self, groupName):
        replicateGroup = self.findReplicateGroup(groupName)
        if replicateGroup:
            return len( replicateGroup["replicates"])  

    def getReplicates(self):
        replisetDict = {}
        for group in self.runOption["replicateGroups"]:
            sample_keys = set() 
            for replicate in group["replicates"]:
                sample_keys.add(replicate["sampleKey"])
            replisetDict[group["groupName"]] = sample_keys
        return replisetDict   

    def findSampleKey(self, replicateName):
        for replicateGroup in self.runOption["replicateGroups"]:
            for replicate in replicateGroup["replicates"]:
                if replicateName.lower() == replicate["replicateName"].lower():
                    return replicate["sampleKey"] 

    def findSampleNames(self, sampleKeys):
        sampleKeyNameDict = {}
        for replicateGroup in self.runOption["replicateGroups"]:
            for replicate in replicateGroup["replicates"]:
                sampleKeyNameDict[int(replicate["sampleKey"])] = replicate["replicateName"]
        sampleNames = []
        for sampleKey in sampleKeys:
            sampleName = sampleKeyNameDict.get(int(sampleKey), None)
            sampleNames.append(sampleName)
        return sampleNames                         

    def getLabel(self):
        return f'{self.bulkData["name"]}_{self.runOption["optionName"]}'
    
    def getShortLabel(self):
        return f'{self.bulkData["shortName"]}_{self.runOption["optionName"]}'    
    
    def getPseudobulkTableName(self):
        return f'tm_temp_pseduobulk_{self.getLabel().lower()}'
    
    def toFile(self, fileName):
        return os.path.join(self.dataDir, fileName)
    
    def getFileNamePrefix(self):
        return self.getLabel()   

    def getPseudobulkDataFile(self):
        return self.toFile(f'{self.getFileNamePrefix()}_Pseudobulk_Data.csv')

    def getQNInputFile(self, groupName):
        numOfReplicates = self.findNumberOfReplicates(groupName)
        return self.toFile(f'{self.getFileNamePrefix()}_{groupName}_QN_Input_{numOfReplicates}_samples.csv')

    def getQNOutputFile(self, groupName):
        numOfReplicates = self.findNumberOfReplicates(groupName)
        return self.toFile(f'{self.getFileNamePrefix()}_{groupName}_QN_Output_{numOfReplicates}_samples.csv')

    def getReplicateGroupAvgFile(self, groupName):
        numOfReplicates = self.findNumberOfReplicates(groupName)
        return self.toFile(f'{self.getFileNamePrefix()}_{groupName}_Replicate_Group_Average_{numOfReplicates}_samples.csv')

    def getAllMergedFile(self):
        return self.toFile(f'All_Structure_Detail.csv')

    def getAllMergedFileTpmOnly(self):
        return self.toFile(f'All_Structure_Avg_Only.csv')

    @staticmethod
    def findBulk(name):
        for bData in PseudobulkConfig.BULK_DATA_LIST:
            if name == bData["name"] or name == bData["shortName"]:
                return bData
            
    @staticmethod
    def findOption(name):
        for option in PseudobulkConfig.RUN_OPTIONS:
            if name == option["optionName"]:
                return option

    BULK_DATA_LIST = []
    bData = {
        "name": "Bladder__urinary_bladder",
        "shortName": "Bladder",
        "tissue": "Bladder",
        "organismPart": "urinary bladder",
        "pivotAFields": ["3_38_F", "3_39_F", "3_56_F", "3_8_M", "3_9_M", "3_10_M"],
        "pivotBFields": [["3_38_F", "3_56_F"], "3_39_F", ["3_8_M", "3_9_M"], "3_10_M"]
    }
    BULK_DATA_LIST.append(bData)

    bData = {
        "name": "Pancreas__endocrine_pancreas",
        "shortName": "Pancreas",
        "tissue": "Pancreas",
        "organismPart": "endocrine pancreas",
        "pivotAFields": ["3_38_F", "3_39_F", "3_8_M", "3_10_M"],
        "pivotBFields": ["3_38_F", "3_39_F", "3_8_M", "3_10_M"]
    }
    BULK_DATA_LIST.append(bData)    

    bData = {
        "name": "Fat__gonadal_fat_pad",
        "shortName": "Fat_gonadal",
        "tissue": "Fat",
        "organismPart": "gonadal fat pad",
        "pivotAFields": ["3_38_F", "3_39_F", "3_56_F", "3_8_M", "3_9_M", "3_10_M", "3_11_M"],
        "pivotBFields": [["3_38_F", "3_56_F"], "3_39_F", ["3_8_M", "3_9_M"], ["3_10_M", "3_11_M"]]
    }
    BULK_DATA_LIST.append(bData)    

    bData = {
        "name": "Fat__interscapular_brown_adipose_tissue",
        "shortName": "Fat_interscapular",
        "tissue": "Fat",
        "organismPart": "interscapular brown adipose tissue",
        "pivotAFields": ["3_38_F", "3_39_F", "3_56_F", "3_8_M", "3_9_M", "3_10_M", "3_11_M"],
        "pivotBFields": [["3_38_F", "3_56_F"], "3_39_F", ["3_8_M", "3_9_M"], ["3_10_M", "3_11_M"]]
    }
    BULK_DATA_LIST.append(bData)  

    bData = {
        "name": "Large_Intestine__ascending_colon",
        "shortName": "Intestine",
        "tissue": "Large_Intestine",
        "organismPart": "ascending colon",
        "pivotAFields": ["3_38_F", "3_39_F", "3_56_F", "3_8_M", "3_9_M", "3_10_M", "3_11_M"],
        "pivotBFields": [["3_38_F", "3_56_F"], "3_39_F", ["3_8_M", "3_9_M"], ["3_10_M", "3_11_M"]]
    }
    BULK_DATA_LIST.append(bData)       

    RUN_OPTIONS = []
    optionA = {
        "optionName": "A",
        "replicateGroups": [
            {
                "groupName": "Female",
                "replicates": [
                    {
                        "replicateName": "3_38_F",
                        "sampleKey": 412627
                    },
                    {
                        "replicateName": "3_39_F",
                        "sampleKey": 415212
                    },
                    {
                        "replicateName": "3_56_F",
                        "sampleKey": 417666
                    }
                ]
            },
            {
                "groupName": "Male",
                "replicates": [
                    {
                        "replicateName": "3_8_M",
                        "sampleKey": 406044
                    },
                    {
                        "replicateName": "3_9_M",
                        "sampleKey": 426746
                    },
                    {
                        "replicateName": "3_10_M",
                        "sampleKey": 436182
                    }
                ]
            }
        ]
    }
    RUN_OPTIONS.append(optionA)

    optionB = {
        "optionName": "B",
        "replicateGroups": [
            {
                "groupName": "Female",
                "replicates": [
                    {
                        "replicateName": "3_38_F_and_3_56_F",
                        "sampleKey": 406717
                    },
                    {
                        "replicateName": "3_39_F",
                        "sampleKey": 419441
                    }
                ]
            },
            {
                "groupName": "Male",
                "replicates": [
                    {
                        "replicateName": "3_8_M_and_3_9_M",
                        "sampleKey": 415397
                    },
                    {
                        "replicateName": "3_10_M",
                        "sampleKey": 432617
                    }
                ]
            }
        ]
    }
    RUN_OPTIONS.append(optionB)

    sex_2_group_query ='''
        WITH data AS (
            SELECT DISTINCT ON (individual)
                sex,
                _sample_key,
                individual,
                sex AS replicateKey
            FROM gxd_htsample
            WHERE _Experiment_key = 99680 AND individual NOT in('3_57_F', '3_11_M')
        )
        SELECT replicateKey,
            string_agg(DISTINCT _sample_key::text, ',') _sample_keys, count(*) num_of_samples,
            string_agg(DISTINCT (_sample_key::text || '|' || individual), ',') sample_key_names
        FROM data
        GROUP BY replicateKey
    '''

    sex_4_group_query ='''
        WITH data AS (
            SELECT DISTINCT ON (bioreplicate_name)
                sex,
                _sample_key,
                CASE
                    WHEN individual IN ('3_38_F', '3_56_F') THEN '3_38_F_and_3_56_F'
                    WHEN individual IN ('3_8_M', '3_9_M') THEN '3_8_M_and_3_9_M'
                    ELSE individual
                END AS bioreplicate_name,
				sex AS replicateKey
            FROM gxd_htsample
            WHERE _Experiment_key = 99680 AND individual NOT in('3_57_F', '3_11_M')
        )
        SELECT replicateKey,
            string_agg(DISTINCT _sample_key::text, ',') _sample_keys, count(*) num_of_samples,
            string_agg(DISTINCT (_sample_key::text || '|' || bioreplicate_name), ',') sample_key_names
        FROM data
        GROUP BY replicateKey
    '''    