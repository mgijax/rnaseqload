##########################################################################
#
# Purpose: Create pre-prcoessing files for each Experiment in Baseline set
#
# Usage: rnaseqPPBaselineload.py
#
# Env Vars:
#	 BASELINEINPUTDIR - post processing files generated from RAW input files
#	 BASELINERAW_INPUTDIR - files downloaded from source
#	 BASELINEOUTPUTDIR - rnaseq bcp files
#	 INSTALLDIR - for path to run_join script	 
#	 EAE_TPMS_PP_FILE_TEMPLATE - path and template for processed eae files
#    EAE_GROUP_PP_FILE_TEMPLATE - path and template for processed eae files
#    AES_SDRF_LOCAL_FILE_TEMPLATE - path and template for processed eae files
#
# Inputs:
#	1. Database: Baseline RNASeq Experiments set
#	2. ArrayExpress files by experiment
#	3. Expression Atlas files by experiment
#	4. Configuration (see rnaseqload.config)
#
# Outputs:
#   BASELINEINPUTDIR
#       xxx.group.txt
#       xxx.tpms.txt
# 
###########################################################################

import os
import sys
import xml.etree.ElementTree as ET
import db
import mgi_utils

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
# Purpose: loads a lookup of samples in the db for the given experiment
# Returns: 0
# Assumes: Nothing
# Effects: queries a database
# Throws: Nothing
#
def loadSamples(expID):
    global sampleInMGI
    sampleInMGI = []

    results = db.sql('''
        select hts.name, hts._Sample_key
        from GXD_HTSample hts, ACC_Accession a
        where hts._Experiment_key = a._Object_key
        and a._MGIType_key = 42 -- experiment
        and a._LogicalDB_key = 189 --ArrayExpress
        and a.preferred = 1
        and a.accID = '%s' 
        ''' % expID, 'auto')
    for r in results:
        key = str.strip(r['name'])
        sampleInMGI.append(key)

    return 0

# end loadSamples()

#
# Purpose: creates EAE_TPMS_PP_FILE_TEMPLATE file in BASELINEINPUTDIR folder for given expID
# Returns: 0 if successful, 1 if unsuccessful
# Assumes: Nothing
# Effects: creates file in filesystem
# Throws: Nothing
#
# format:
#   ensembm ID
#   marker key
#   marker symbol
#   each group (g1, g2, etc. value = 3rd value avg QN TPM)
#

def ppEAETpmsFile(expID):

    print('in ppEAETpmsFile(expID): %s' % expID)

    #  read the input file
    eaeFile = tpmsTemplate % expID
    print('eaeFile: %s' % eaeFile)
    try:
        fpEae = open(eaeFile, 'r')
    except:
        print('skipping: missing -tpms.tsv file: %s' % (expID))
        return 1 # file does not exist

    #  create the output file
    ppFile = tpmsPPTemplate % expID
    try:
        fpPP = open(ppFile, 'w')
    except:
        return 1 # file does not exist

    ensemblDict = {}
    results = db.sql('select accid, _object_key from acc_accession where _logicaldb_key = 60 and _mgitype_key = 2 and preferred = 1', 'auto')
    for r in results:
        key = r['accid']
        value = r['_object_key']
        if key not in ensemblDict:
            ensemblDict[key] = []
        ensemblDict[key].append(value)

    # read the header from fpEae and create header for fpPP
    headerList = str.split(fpEae.readline(), '\t')
    groupSet = []
    for h in headerList[2:]:
        groupSet.append(str.strip(h))
    fpPP.write('ensembl_id\t_marker_key\tsymbol\t' + '\t'.join(groupSet) + '\n')

    # iterate thru the fpEae input file
    for line in fpEae.readlines():

        tokens = str.split(line, '\t')
        ensemblID = str.strip(tokens[0])

        # if ensemblID is not in MGI, then set markerKey = 0
        # will handle this later during TPMS processing
        if ensemblID in ensemblDict:
            markerKey = ensemblDict[ensemblID][0]
        else:
            markerKey = 0
        markerSymbol = str.strip(tokens[1])

        fpPP.write('%s\t%s\t%s' % (ensemblID, markerKey, markerSymbol))

        for g in range(len(groupSet)):
            values = tokens[g+2].split(',')
            fpPP.write('\t' + values[2])

        fpPP.write('\n')

    fpEae.close();
    fpPP.close();

    return 0

# end ppEAETpmsFile()

#
# Purpose: Reads the aesTemplate(AES_SDRF_LOCAL_FILE_TEMPLATE) and processes the runID -> Sample ID match
# Returns: 0 if successful, 1 if unsuccessful
# Assumes: Nothing
# Effects: creates the runToSampleDict file (runID -> Sample ID)
# Throws: Nothing
#
def ppAESSdrfFile(expID):

    global rawRunList, runToSampleDict

    print('in ppAESSdrfFile(expID): %s' % expID)

    #  read the input file
    aesFile = aesTemplate % expID
    try:
        fpAes = open(aesFile, 'r')
    except:
        print('skiiping: missing .sdrf.txt file: %s' % (expID))
        return 1 # file does not exist

    # process the header line
    #
    headerList = str.split(fpAes.readline(), '\t')
    if headerList == ['']: # means file is empty
        print ('skipping: missing header: %s' % (expID))
        return 1

    # load sampleInMGI() for expID
    loadSamples(expID)

    # find the idx of the columns we want - they are not ordered
    sourceSampleIDX = None
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
    if sourceSampleIDX == None:
        print('skipping: missing Source Name column: %s' % (expID))
        return 1
    if enaRunIDX == None:
        print('skipping: missing ENA_RUN column: %s' % (expID))
        return 1

    # iterate thru the fpEae input file
    for line in fpAes.readlines():

        tokens = str.split(line, '\t')

        sourceSample = None
        enaSample = None
        enaRun = None

        if sourceSampleIDX != None:
            sourceSample = str.strip(tokens[sourceSampleIDX])

        # not every sample file contais enaSampleIDX column
        if enaSampleIDX != None:
            enaSample = str.strip(tokens[enaSampleIDX])

        if enaRunIDX != None:
            enaRun = str.strip(tokens[enaRunIDX])

        # skip if this source is a duplicate; but don't report
        if enaRun in rawRunList:
            #print('skipping: enaRun already processed: %s,%s' % (expID, enaRun))
            continue

        if sourceSample not in sampleInMGI:
            if enaSample != None and enaSample in sampleInMGI:
                sourceSample = enaSample
            else:
                print('skipping: sample is not in MGI: %s, sourceSample = %s, enaSample = %s' % (expID, sourceSample, str(enaSample)))
                continue

        rawRunList.append(enaRun)

        key = enaRun
        value = sourceSample
        if key not in runToSampleDict:
            runToSampleDict[key] = []
        runToSampleDict[key].append(value)

    fpAes.close();

    return 0

# end ppAESSdrfFile()

#
# Purpose: Reads groupTermplate (EAE_GROUP_LOCAL_FILE_TEMPLATE) and runToSampleDict()
# Returns: 0 if successful, 1 if unsuccessful
# Assumes: Nothing
# Effects: creates EAE_GROUP_PP_FILE_TEMPLATE file in BASELINEINPUTDIR folder for given expID
# Throws: Nothing
#
# format:
# 	group ID 
# 	label 
# 	Run IDs
#   Sample IDs
#

def ppEAEGroupFile(expID):

    print('in ppEAEGroupFile(expID): %s' % expID)

    #  read the input file
    try:
        eaeFile = groupTemplate % expID
    except:
        print('skipping: missing -configuration.xml: %s' % (expID))
        return 1 # file does not exist

    #  create the output file
    ppFile = groupPPTemplate % expID
    try:
        fpPP = open(ppFile, 'w')
    except:
        return 1 # file does not exist

    #
    # eaeFile is in XML format
    #
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
            runID = child.text
            sampleID = 'missing'
            if runID in runToSampleDict:
                sampleID = runToSampleDict[runID][0]
            fpPP.write('%s\t%s\t%s\t%s\n' % (id, label, runID, sampleID))

    fpPP.close();

    return 0

# end ppEAEGroupFile()
 
#
# Purpose: pre processing
#   read EAE tpms, group and AES sdrf files
#   generate tpms output file, group output file
# Returns: 0 if successful; handles return codes from all functions it calls 
# Assumes: Nothing
# Effects:
# Throws: Nothing
#
def process():
    global rawRunList

    results = db.sql('''
        select a.accid
        from ACC_Accession a, MGI_Set s, MGI_SetMember sm
        where s.name = 'Baseline RNASeq Load Experiments'
        and s._Set_key = sm._Set_key
        and sm._Object_key = a._Object_key
        and a._MGIType_key = 42 --GXD_HTExperiment
        and a._LogicalDB_key = 189
        and a.preferred = 1
        ''', 'auto')

    #
    # for each expID in the MGI_Set:
    #
    for r in results:

        expID = str.strip(r['accid'])

        # process the eae/tpms file for this expID
        rc = ppEAETpmsFile(expID)
        if rc != 0:
            print('processing EAE tpms file returned rc %s, skipping file for %s' % (rc, expID))
            continue

        # process the aes/sdrf file for this expID to create the runToSampleDict
        rc = ppAESSdrfFile(expID)
        if rc != 0:
            print('processing AES sdrf file returned rc %s, skipping file for %s' % (rc, expID))
            continue

        # process the eae/group file for this expID
        rc = ppEAEGroupFile(expID)
        if rc != 0:
            print('processing EAE group file returned rc %s, skipping file for %s' % (rc, expID))
            continue

    return 0

# end process()

#
# Main
#

print('start time: %s' %  mgi_utils.date())
if process() != 0:
     exit(1, 'Error in process()\n')
print('end time: %s' %  mgi_utils.date())
