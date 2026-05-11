##########################################################################
#
# Purpose: Create pre-prcoessing files for each Experiment in Baseline set
#
# Usage: rnaseqPPBaselineload.py
#
# Env Vars:
#	 LOGDIR 
#	 BASELINEINPUTDIR - intermediate files generated from RAW input files
#	 BASELINERAW_INPUTDIR - files downloaded from source
#	 BASELINEOUTPUTDIR - rnaseq bcp files
#	 INSTALLDIR - for path to run_join script	 
#	 LOG_CUR - curation log
#	 LOG_DIAG - diagnostic log
#	 EAE_TPMS_PP_FILE_TEMPLATE - path and template for processed eae files
#    EAE_GROUP_PP_FILE_TEMPLATE - path and template for processed eae files
#    AES_SDRF_PP_FILE_TEMPLATE - path and template for processed eae files
#
# Inputs:
#	1. Database: Baseline RNASeq Experiment set
#	2. ArrayExpress files by experiment
#	3. Expression Atlas files by experiment
#	4. Configuration (see rnaseqload.config)
#
# Outputs:
#	 1. preprocessed file for each experiment aes and eae
# 
###########################################################################

import os	 	# for system to execute bcp, getenv
import sys 	 	# to flush stdout
import xml.etree.ElementTree as ET
import db
import mgi_utils 	# for log start/end timestamp

# Expression Atlas Experiment file Template - name of file stored locally
tpmsTemplate = '%s' % os.getenv('EAE_TPMS_LOCAL_FILE_TEMPLATE')
tpmsPPTemplate = '%s' % os.getenv('EAE_TPMS_PP_FILE_TEMPLATE')
groupTemplate = '%s' % os.getenv('EAE_GROUP_LOCAL_FILE_TEMPLATE')
groupPPTemplate = '%s' % os.getenv('EAE_GROUP_PP_FILE_TEMPLATE')
aesTemplate = '%s' % os.getenv('AES_SDRF_LOCAL_FILE_TEMPLATE')

# unique set of raw samples
rawRunList = []

# which samples to this run belong to
# run -> sample
runToSampleDict = {}

#
# Purpose: init lookups
# Returns: 0
# Assumes: Nothing
# Effects: opens a database connection, queries a database
# Throws: Nothing
#
def init():
    global rnaSeqSetResults

    # create the result set of ids to load
    rnaSeqSetResults = db.sql('''
        select a.accid
        from ACC_Accession a, MGI_Set s, MGI_SetMember sm
        where s.name = 'Baseline RNASeq Load Experiment'
        and s._Set_key = sm._Set_key
        and sm._Object_key = a._Object_key
        and a._MGIType_key = 42 --GXD_HTExperiment
        and a._LogicalDB_key = 189
        and a.preferred = 1
        ''', 'auto')

    return 0

# end init() -------------------------------------------------

#
# Purpose: creates EAE_TPMS_PP_FILE_TEMPLATE file in BASELINEINPUTDIR folder for given expID
# Returns:
# Assumes: Nothing
# Effects: creates file in filesystem
# Throws: Nothing
#
# format:
#   ensembm ID
#   marker key
#   marker symbol
#   g1 = 3rd value
#   g2 = 3rd value
#

def ppEAETpmsFile(expID):

    print('in ppEAETpmsFile(expID): %s' % expID)

    #  read the input file
    eaeFile = tpmsTemplate % expID
    print('eaeFile: %s' % eaeFile)
    try:
        fpEae = open(eaeFile, 'r')
    except:
        return 1 # file does not exist

    #  create the output file
    ppFile = tpmsPPTemplate % expID
    try:
        fpPP = open(ppFile, 'w')
    except:
        return 1 # file does not exist

    # iterate thru the fpEae input file
    for line in fpEae.readlines():

        tokens = str.split(line, '\t')

        if tokens[0] == 'GeneID':
            continue

        ensemblID = str.strip(tokens[0])
        markerSymbol = str.strip(tokens[1])
        g1All = str.strip(tokens[2])
        g2All = str.strip(tokens[3])
        
        #print(tokens)
        tokensG1 = g1All.split(',')
        tokensG2 = g2All.split(',')
        g1 = tokensG1[2]
        g2 = tokensG2[2]

        # write to the fpPP file
        fpPP.write('%s\t%s\t%s\t%s\n' % (ensemblID, markerSymbol, g1, g2))

    fpEae.close();
    fpPP.close();

    return 0

# end ppEAETpmsFile ()--------------------------------------------

#
# Purpose: creates EAE_GROUP_PP_FILE_TEMPLATE file in BASELINEINPUTDIR folder for given expID
# Returns:
# Assumes: Nothing
# Effects: creates file in filesystem
# Throws: Nothing
#
# format:
# 	group ID 
# 	label 
# 	Run IDs
#

def ppEAEGroupFile(expID):

    print('in ppEAEGroupFile(expID): %s' % expID)

    #  read the input file
    eaeFile = groupTemplate % expID
    print('eaeFile: %s' % eaeFile)

    #  create the output file
    ppFile = groupPPTemplate % expID
    try:
        fpPP = open(ppFile, 'w')
    except:
        return 1 # file does not exist

    # iterate thru the eaeFile xml file
    #        <assay_group id="g1" label="brown adipose tissue">
    #            <assay>ERR4193656</assay>
    #            <assay>ERR4193654</assay>
    #            <assay>ERR4193655</assay>
    #        </assay_group>

    tree = ET.parse(eaeFile)
    root = tree.getroot()
    assay_groups = root.findall('.//assay_group')
    for ag in assay_groups:
        id = ag.get('id')
        label = ag.get('label')
        runids = []
        for child in ag:
            #print(child.tag, child.text)
            # write to the fpPP output file
            runID = child.text
            sampleID = 'missing'
            if runID in runToSampleDict:
                sampleID = runToSampleDict[runID][0]
            fpPP.write('%s\t%s\t%s\t%s\n' % (id, label, runID, sampleID))

    fpPP.close();

    return 0

# end ppEAEGroupFile ()--------------------------------------------
 
#
# Purpose: creates the run to sample relationship : runToSampleDict
# Returns:
# Assumes: Nothing
# Effects: creates file in filesystem
# Throws: Nothing
#

def ppAESSdrfFile(expID):

    global rawRunList, runToSampleDict

    print('in ppAESSdrfFile(expID): %s' % expID)

    #  read the input file
    aesFile = aesTemplate % expID
    print('aesFile: %s' % aesFile)
    try:
        fpAes = open(aesFile, 'r')
    except:
        return 1 # file does not exist

    # process the header line
    #
    headerList = str.split(fpAes.readline(), '\t')
    if headerList == ['']: # means file is empty
        print ('skipping: missing header: %s\n' % (expID))
        return 1

    # find the idx of the columns we want - they are not ordered
    #
    enaSampleIDX = None
    enaSampleIDX = None
    enaRunIDX = None
    for idx, colName in enumerate(headerList):
        colName = str.strip(colName)
        if str.find(colName, 'Source Name') != -1:
            sourceSampleIDX = idx
        elif str.find(colName, 'ENA_SAMPLE') != -1:
            enaSampleIDX = idx
        elif str.find(colName, 'ENA_RUN') != -1:
            enaRunIDX = idx
    
    # iterate thru the fpEae input file
    for line in fpAes.readlines():

        tokens = str.split(line, '\t')

        if tokens[0] == 'Source Name':
            continue

        #print(tokens)

        sourceSample = None
        enaSample = None
        enaRun = None

        if sourceSampleIDX != None:
            sourceSample = str.strip(tokens[sourceSampleIDX])

        if enaSampleIDX != None:
            enaSample = str.strip(tokens[enaSampleIDX])

        if enaRunIDX == None:
            print('skipping: could not find ENA_RUN in header: %s\n' % (expID))
            return 1

        if str.find(sourceSample, 'ERS') == -1:
            sourceSample = enaSample

        if enaRunIDX != None:
            enaRun = str.strip(tokens[enaRunIDX])

        if sourceSample == None or enaRun == None:
            print('skipping : missing sourceSample : %s, enaRun : %s\n' % (sourceSample, enaRun))
            return 1
            
        # skip of this source is a duplicate
        if enaRun in rawRunList:
            #print('skipping : enaRun already processed : %s,%s' % (expID, enaRun))
            continue

        rawRunList.append(enaRun)

        print('Source/Run: ', sourceSample, enaRun)
        key = enaRun
        value = sourceSample
        if key not in runToSampleDict:
            runToSampleDict[key] = []
        runToSampleDict[key].append(value)

    fpAes.close();

    return 0

# end ppAESSdrfFile ()--------------------------------------------

# Purpose: main processing function. 
#   Gets the list of experiment IDs from the 'Baseline RNASeq Load Experiment' set 
#   Processes experiments one by one
# Returns: 0 if successful; handles return codes from all functions it calls 
# Assumes: Nothing
# Effects:
# Throws: Nothing
#
def process():
    global rawRunList

    #
    # for each expID in the MGI_Set:
    #
    for r in rnaSeqSetResults:

        expID = str.strip(r['accid'])

        # preprocess the eae/tpms file for this expID
        #rc = ppEAETpmsFile(expID)
        #if rc != 0:
        #    print('preprocessing EAE tpms file returned rc %s, skipping file for %s' % (rc, expID))
        #    continue

        # preprocess the aes/sdrf file for this expID to create the runToSampleDict
        rc = ppAESSdrfFile(expID)
        if rc != 0:
            print('preprocessing AES sdrf file returned rc %s, skipping file for %s' % (rc, expID))
            continue

        # preprocess the eae/group file for this expID
        rc = ppEAEGroupFile(expID)
        if rc != 0:
            print('preprocessing EAE group file returned rc %s, skipping file for %s' % (rc, expID))
            continue

    #print(runToSampleDict)
    return 0

# end process ()--------------------------------------------

#
# Main
#

print('Start time: %s' %  mgi_utils.date())
if init() != 0:
     exit(1, 'Error in  init \n')

if process() != 0:
     exit(1, 'Error in  process \n')

print('End time: %s' %  mgi_utils.date())

