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
#	 AES_PP_FILE_TEMPLATE - preprocessed to just runID, sampleID
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

# paths to input and two output files
logDir = os.getenv('BASELINELOG_PP')
inputDir = os.getenv('BASELINEINPUTDIR')
rawInputDir = os.getenv('BASELINERAW_INPUTDIR')
outputDir = os.getenv('BASELINEOUTPUTDIR')

# Expression Atlas Experiment file Template - name of file stored locally
tpmsTemplate = '%s' % os.getenv('EAE_TPMS_LOCAL_FILE_TEMPLATE')
tpmsPPTemplate = '%s' % os.getenv('EAE_TPMS_PP_FILE_TEMPLATE')
groupTemplate = '%s' % os.getenv('EAE_GROUP_LOCAL_FILE_TEMPLATE')
groupPPTemplate = '%s' % os.getenv('EAE_GROUP_PP_FILE_TEMPLATE')
aesTemplate = '%s' % os.getenv('AES_SDRF_LOCAL_FILE_TEMPLATE')
aesTemplate = '%s' % os.getenv('AES_SDRF_PP_FILE_TEMPLATE')

#
# Purpose:  create db connection, init lookups
# Returns: 0
# Assumes: Nothing
# Effects: opens a database connection, queries a database
# Throws: Nothing
#
def init():
    global fpLog, rnaSeqSetResults

    fpLog = open (logDir, 'w')

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

    fpLog.write('in ppEAETpmsFile(expID): %s' % expID)

    #  read the input file
    eaeFile = tpmsTemplate % expID
    fpLog.write('eaeFile: %s' % eaeFile)
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
        
        #fpLog.write(tokens)
        tokensG1 = g1All.split(',')
        tokensG2 = g2All.split(',')
        g1 = tokensG1[2]
        g2 = tokensG2[2]
        markerKey = ensemblMarkerDict[ensemblID][0]['_object_key']
        #fpLog.write(markerKey, g1, g2)

        # write to the fpPP file
        fpPP.write('%s\t%s\t%s\t%s\t%s\n' % (ensemblID, markerKey, markerSymbol, g1, g2))

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

    fpLog.write('in ppEAEGroupFile(expID): %s' % expID)

    #  read the input file
    eaeFile = groupTemplate % expID
    fpLog.write('eaeFile: %s' % eaeFile)

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
            fpLog.write(child.tag, child.text)
            runids.append(child.text)

        # write to the fpPP output file
        #fpLog.write('%s\t%s\t%s\n' % (id, label, runids))
        fpPP.write('%s\t%s\t%s\n' % (id, label, '|'.join(runids)))

    fpPP.close();

    return 0

# end ppEAEGroupFile ()--------------------------------------------
 
#
# Purpose: creates AES_SDRF_PP_FILE_TEMPLATE file in BASELINEINPUTDIR folder for given expID
# Returns:
# Assumes: Nothing
# Effects: creates file in filesystem
# Throws: Nothing
#
# format:
# 1 Source Name     
# 2 Comment[ENA_SAMPLE]     
# 34 Comment[ENA_RUN]        
#

def ppAESSdrfFile(expID):

    fpLog.write('in ppAESSdrfFile(expID): %s' % expID)

    #  read the input file
    eaeFile = aesTemplate % expID
    fpLog.write('eaeFile: %s' % eaeFile)
    try:
        fpEae = open(eaeFile, 'r')
    except:
        return 1 # file does not exist

    #  create the output file
    ppFile = aesPPTemplate % expID
    try:
        fpPP = open(ppFile, 'w')
    except:
        return 1 # file does not exist

    # iterate thru the fpEae input file
    for line in fpEae.readlines():

        tokens = str.split(line, '\t')

        if tokens[0] == 'Source Name':
            continue

        sourceName = tokens[0]
        enaSample = tokens[1]
        enaRun = tokens[33]
        fpLog.write(sourceName, enaSample, enaRun)

    fpAes.close();
    fpPP.close();

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
    global experimentNotInDbList

    #
    # for each expID in the MGI_Set:
    #
    for r in rnaSeqSetResults:
        expID = str.strip(r['accid'])

        # preprocess the eae/tpms file for this expID
        sys.stdout.flush()
        rc = ppEAETpmsFile(expID)
        if rc != 0:
            fpLog.write('''preprocessing EAE tpms file returned rc %s, skipping file for %s''' % (rc, expID))
            continue

        # preprocess the eae/group file for this expID
        sys.stdout.flush()
        rc = ppEAEGroupFile(expID)
        if rc != 0:
            fpLog.write('''preprocessing EAE group file returned rc %s, skipping file for %s''' % (rc, expID))
            continue

        # preprocess the aes/sdrf file for this expID
        sys.stdout.flush()
        rc = ppAESSdrfFile(expID)
        if rc != 0:
            fpLog.write('''preprocessing AES sdrf file returned rc %s, skipping file for %s''' % (rc, expID))
            continue

    return 0

# end process ()--------------------------------------------

#
# Main
#

fpLog.write('Start time: %s' %  mgi_utils.date())
if init() != 0:
     exit(1, 'Error in  init \n' )

if process() != 0:
     exit(1, 'Error in  process \n' )

fpLog.write('End time: %s' %  mgi_utils.date())
fpLog.close()

