##########################################################################
#
# Purpose: For each experiment ID in the 'RNA Seq Load Experiment' MGI_Set
#	Create intermediate files for EAE and AES files (pulls out just the data 
#	    needed)
#	Join the two sets of intermediate files
#	Determine biological replicates
#	Create a pandas DataFrame from the joined file and the biological
#	    replicate info
#	Run Quantile Normalization
#	Create bcp files
#	Execute bcp
#  	Write QC reports
#
# Usage: rnaseqload.py
# Env Vars:
#	 1. LOGDIR 
#	 2a. INPUTDIR - intermediate files generated from RAW input files
#	 2b. RAW_INPUTDIR - files downloaded from source
#	 3. OUTPUTDIR - rnaseq bcp files
#	 4. INSTALLDIR - for path to run_join script	 
#	 5. LOG_CUR - curation log
#	 6. LOG_DIAG - diagnostic log
#	 7a. RNASEQ_BCP - bcp filename suffix for rnaseq table - 
#			we will append expID
#	 7b. COMBINED_BCP - bcp filename suffix for combined table - 
#			we will append expID
#	 8. PG_DBUTILS - for path to bcpin.csh script
#	 9. AES_LOCAL_FILE_TEMPLATE - path and template for downloaded aes files
#	 10. AES_PP_FILE_TEMPLATE - preprocessed to just runID, sampleID
#	 11. EAE_LOCAL_FILE_TEMPLATE - path and template for downloaded eae files
#	 12. EAE_PP_FILE_TEMPLATE - preprocessed to just geneID, runID, TPM
#	 13. JOINED_PP_FILE_TEMPLATE - path & template for joined aes/eae files
# Inputs:
#	1. Database: RNA Seq Experiment set
#	2. ArrayExpress files by experiment
#	3. Expression Atlas files by experiment
#	4. Configuration (see rnaseqload.config)
#
# Outputs:
#	 1. aes and eae files from source
#	 2. preprocessed file for each experiment aes and eae
#	 3. joined aes and eae file for each experiment
#	 4. bcp files - one per experiment
#	 5. curator and diagnostic log
# 
# Exit Codes:
#
#      0:  Successful completion
#      1:  An exception occurred
#
#  Assumes:  Nothing
#
#  Notes:  None
#
# HISTORY
# 
#       sc - fixed bug in determining non-relevant samples
#         - https://mgi-jira.atlassian.net/browse/WTS2-514
#         - more RNA seq data (plus more EA links) (TR13314)
###########################################################################

import os	 	# for system to execute bcp, getenv
import mgi_utils 	# for log start/end timestamp
import loadlib	 	# for bcp friendly date/timestamp
import string
import db
import sys 	 	# to flush stdout
import time	 	# used for its time.time() function (for timestamps)
import subprocess
import numpy as np	# used for stddev	
import pandas as pd	# used for QN
import quantileNormalize # module within this product

# paths to input and two output files
logDir =  os.getenv('LOGDIR')
inputDir =  os.getenv('INPUTDIR')
rawInputDir = os.getenv('RAW_INPUTDIR')
outputDir = os.getenv('OUTPUTDIR')
binDir = '%s/bin' % os.getenv('INSTALLDIR')

# curation and diagnostic logs
fpCur = open (os.getenv('LOG_CUR'), 'a')
fpDiag = open (os.getenv('LOG_DIAG'), 'a')

# bcp stuff
rnaSeqBcpFile = os.getenv('RNASEQ_BCP')
fpRnaSeqBcp = None # this will get assigned later, once per experiment
combinedBcpFile = os.getenv('COMBINED_BCP')
fpCombinedBcp = None  # this will get assigned later, once per experiment

bcpCommand = os.getenv('PG_DBUTILS') + '/bin/bcpin.csh'
bcpCommandList = []
rnaSeqTable = 'GXD_HTSample_RNASeq'
combinedTable = 'GXD_HTSample_RNASeqCombined'

# status stuff
stdDevCutoff = float(os.getenv('STDDEV_CUTOFF'))
highAveStdDevList = []

# {expID|sampleID:[list of all tpm aveSD for sample], ...}
sampleAveSDDict = {}

# special report for Richard
fpStudentRpt = open('%s/studentRnaSeq.txt' % outputDir, 'w')

# constants
TAB = '\t'
CRT = '\n'
PIPE = '|'

# load date
loaddate = loadlib.loaddate

# rnaseqload MGI_User
createdByKey = 1613

# ArrayExpress Sample File Template  - name of file stored locally
aesTemplate = '%s' % os.getenv('AES_LOCAL_FILE_TEMPLATE')
# ArrayExpress Sample File Preprocessed Template
aesPPTemplate = '%s' % os.getenv('AES_PP_FILE_TEMPLATE')

# Expression Atlas Experiment file Template - name of file stored locally
eaeTemplate  =  '%s' % os.getenv('EAE_LOCAL_FILE_TEMPLATE')
# Expression Atlas Experiment file Preprocessed Template
eaePPTemplate  =  '%s' % os.getenv('EAE_PP_FILE_TEMPLATE')

# the joined PP files template
joinedPPTemplate = '%s' % os.getenv('JOINED_PP_FILE_TEMPLATE')

# GXT HT Experiment IDs in the database
#experimentInDbSet = set()
experimentInDbDict = {}

# samples with strain J:DO; genotypeKey = 90560 in the db
JDOSampleSet = set()

# samples flagged as relevant in the db
relevantSampleSet = set()

# ensembl IDs assoc w/markers
# used to see if the ensembl ID is in the database
# Sets are faster; lists are sequential scan
ensemblMarkerSet = set()

# ensembl IDs assoc w/markers, to get the marker Key
# if multi markers, the last marker will be stored
# ens assoc/multi markers will be filtered out
# and will not be looked up in this dictionary
ensemblMarkerDict = {}

# ensembl IDs assoc w/sequences
ensemblSequenceSet = set()

# ensembl IDs assoc w/ >1 marker
multiMarkerEnsemblDict = {}

# > 1 Ens IDs for a marker
symbolToMultiEnsIdDict = {} # {symbol:[list of ensIDs], ...}

# above set by ensID
# {ensID:symbol, ...}
multiEnsemblMarkerDict = {}

#
# Level bins
#
HIGH = 50430889
MED = 50430890
LOW = 50430891
BELOW_CUTOFF = 50430892
 
#
# QC data structures
#

# number of experiments in Connie's input file
exptCt = 0

# expIDs in File not in database
experimentNotInDbList = []

# runID not in ArrayExpress sampleID/runID file
runIdNotInAEList = []

# runIDs not in Expression Atlas geneID/runID/TPM file
runIdNotInEAList = []

# sampleID from ArrayExpress sampleID/runID file  not in database
sampleIdNotInDbList = []

# run ID, sample ID columns not in the AES file
noRunOrSampleColList = []

# J:DO strain samples skipped
JDOSkippedList = []

# non-relevant samples skipped
nonRelSkippedList = []

# empty TPM values set to 0.0
noTpmValueList = []

# ensembl IDs assoc w/multi markers
multiMarkerEnsemblList = []

# markers assoc w/multi ensIDs
multiEnsemblMarkerList = []

# ensembl ID is not in MGI at all
ensemblNotInMGIList = []

# ensembl ID is only associated with a sequence (not a marker)
ensemblIsOrphanList = []

# ambiguous expID/sampleID in the database
ambiguousSampleInDbList = []

# next available rnaSeq Key - set by init()
rnaSeqKey = None

# next available combined key - set by init()
combinedKey = None

#
# Purpose:  create db connection, init lookups
# Returns: 0
# Assumes: Nothing
# Effects: opens a database connection, queries a database
# Throws: Nothing
#
def init():
    global experimentInDbDict, sampleInDbList, JDOSampleSet
    global relevantSampleSet, ensemblMarkerSet, ensemblSequenceSet
    global multiEnsemblMarkerDict, multiMarkerEnsemblDict, symbolToMultiEnsIdDict
    global rnaSeqKey, combinedKey, ensemblMarkerDict, rnaSeqSetResults

    db.useOneConnection(1)

    # We have truncated the tables
    rnaSeqKey = 1
    combinedKey = 1

    results = db.sql('''select accid, _object_key
        from ACC_Accession
        where _MGIType_key = 42
        and preferred = 1''', 'auto')
    for r in results:
        experimentInDbDict[r['accid']] = r['_object_key']

    results =  db.sql('''select name
        from GXD_HTSample
        where _Genotype_key = 90560''', 'auto')
    for r in results:
        JDOSampleSet.add(str.strip(r['name']))

    results =  db.sql('''select a.accid as exptID, s.name
        from GXD_HTSample s, ACC_Accession a
        where s._Relevance_key = 20475450
        and s._experiment_key = a._object_key
        and a._mgitype_key = 42
        order by a.accid, s.name''', 'auto')
    for r in results:
        relevantSampleSet.add('%s|%s' % (str.strip(r['exptID']), str.strip(r['name'])))
    results = db.sql('''select accid, _Object_key
        from ACC_Accession
        where _LogicalDB_key =  60
        and _MGIType_key = 2
        and preferred = 1 ''', 'auto')
    for r in results:
        accid = r['accid']
        ensemblMarkerSet.add(accid)
        ensemblMarkerDict[accid] = r['_Object_key']

    db.sql('''select accid
        into temporary table multiEns
        from ACC_Accession
        where _LogicalDB_key =  60
        and _MGIType_key = 2
        and preferred = 1
        group by accid having count(*) > 1 ''', None)

    db.sql('''create index idx1 on multiEns(accid)''')

    results = db.sql('''
        select a.accid, m.symbol
        from multiEns me, ACC_Accession a, MRK_Marker m
        where me.accid = a.accid
        and a._LogicalDB_key =  60
        and a._MGIType_key = 2
        and a.preferred = 1
        and a._Object_key = m._Marker_key
        and m._MarkerStatus_key = 1
        ''', 'auto')

    for r in results:
        ensID = r['accid']
        symbol = r['symbol']
        if ensID not in multiMarkerEnsemblDict:
            multiMarkerEnsemblDict[ensID] = []
        multiMarkerEnsemblDict[ensID].append(symbol)

    db.sql('''select _Object_key
        into temporary table multi
        from ACC_Accession
        where _MGIType_key = 2
        and _LogicalDB_key = 60
        group by _Object_key 
        having count(*) > 1''', None)
    
    db.sql('create index idx2 on multi(_Object_key)')
    
    results = db.sql('''
        select a.accid, m.symbol
        from ACC_Accession a, multi mm, MRK_Marker m
        where a._MGIType_key = 2
        and a._LogicalDB_key = 60
        and a._Object_key = mm._Object_key
        and a._Object_key = m._Marker_key
        and m._MarkerStatus_key = 1
        order by m.symbol
        ''', 'auto')

    for r in results:
        ensID = r['accid']
        symbol = r['symbol']
        if symbol not in symbolToMultiEnsIdDict:
            symbolToMultiEnsIdDict[symbol] = []
        symbolToMultiEnsIdDict[symbol].append(ensID)
        multiEnsemblMarkerDict[ensID] = symbol

    results = db.sql('''select accid
        from ACC_Accession
        where _LogicalDB_key =  60
        and _MGIType_key = 19
        and preferred = 1 ''', 'auto')
    for r in results:
        ensemblSequenceSet.add(r['accid'])

    # create the result set of ids to load
    rnaSeqSetResults = db.sql('''select a.accid
    from ACC_Accession a, MGI_Set s, MGI_SetMember sm
    where s.name = 'RNASeq Load Experiments'
    and s._Set_key = sm._Set_key
    and sm._Object_key = a._Object_key
    and a._MGIType_key = 42 --GXD_HTExperiment
    and a._LogicalDB_key = 189
    and a.preferred = 1''', 'auto')

    return 0

# end init() -------------------------------------------------

#
# Purpose: loads a lookup of samples in the db for the given experiment
#       we did this because sample names are not uniq across experiments
# Returns: 1 if environment variable not set
# Assumes: Nothing
# Effects: queries a database
# Throws: Nothing
#
def loadSampleInDbDict(expID):
    # ArrayExpress (excludes GEO) GXD HT Sample names for expID in the database
    # {sampleID:[sampleKey1, sampleKeyN], ...}
    global sampleInDbDict
    sampleInDbDict = {}

    results = db.sql('''select hts.name, hts._Sample_key
        from GXD_HTSample hts, ACC_Accession a
        where hts._Experiment_key = a._Object_key
        and a._MGIType_key = 42 -- experiment
        and a._LogicalDB_key = 189 --ArrayExpress
        and a.preferred = 1
        and a.accID = '%s' ''' % expID, 'auto')
    for r in results:
        sampleID = str.strip(r['name'])
        if sampleID not in sampleInDbDict:
            sampleInDbDict[sampleID] = []
        sampleInDbDict[sampleID].append(r['_Sample_key'])

    return 0

# end loadSampleInDbDict() -------------------------------------------------

#
# Purpose: Initializes the current GXD_HTSample_RNASeq bcp file descriptor
#	There is one file per experiment
# Returns: 0
# Assumes: Nothing
# Effects: creates file in filesystem
# Throws: Nothing
#
def initRnaSeqBcp(bcpFile):
    global fpRnaSeqBcp

    fpRnaSeqBcp = open(bcpFile, 'w')
    
    return 0

# end initRnaSeqBcp() -------------------------------------------------

#
# Purpose: Initializes the current GXD_HTSample_RNASeqCombined bcp file descriptor. 
#       There is one file per experiment
# Returns: 0
# Assumes: Nothing
# Effects: creates file in filesystem
# Throws: Nothing
#
def initCombinedBcp(bcpFile):
    global fpCombinedBcp

    fpCombinedBcp = open(bcpFile, 'w')

    return 0

# end initCombinedBcp() -------------------------------------------------

#
# Purpose: Closes the current GXD_HTSample_RNASeq bcp file descriptor.
# 
def closeRnaSeqBcp():

    fpRnaSeqBcp.close()

    return 0

# end closeRnaSeqBcp() -------------------------------------------------

#
# Purpose: Closes the current GXD_HTSample_RNASeqCombined bcp file descriptor.
# 
def closeCombinedBcp():

    fpCombinedBcp.close()

    return 0

# end closeCombinedBcp() -------------------------------------------------

#
# Purpose: calculates the average of numList
# Returns: the calculated average
# 
def calcAve(numList, precision):
    tSum = sum(numList)
    tLen =  max(len(numList), 1)
    ave = tSum / tLen
    fAve = float(ave)

    return round(fAve, precision)

# end calcAve () -------------------------------------------------

#
# Purpose: creates file descriptor and reads the ArrayExpress file pulling out 
#	the RunID/SamplID mapping. Writes to the AES preprocessed intermediate file
# Returns: 0 if successful, 1 if AES file does not exist, 2 if AES file is empty
#	3 if No run column in file, 4 if no sample column in file, 5 if ambiguous
#	expID/sampleID in database
# Assumes: Nothing
# Effects: creates file in filesystem
# Throws: Nothing
#

def ppAESFile(expID):
    global noRunOrSampleColList, sampleIdNotInDbList, currentAesPPFile
    global aesRunIdSet, ambiguousSampleInDbList

    #print 'in ppAESFile(expID): %s' % expID
    # run IDs in the AES set
    aesRunIdSet = set()

    # QC expID/sampleID pairs - if more than one sample record in database same sampleName/ExpID
    # report and skip experiment
    # {sampleID:dbCt, ...}
    sampleDbCtDict = {}

    # 
    # create aes input file descriptor for this experiment
    #
    aesFile = aesTemplate % expID
    #print 'aesFile: %s' % aesFile
    try:
        fpAes = open(aesFile, 'r')
    except:
        return 1 # file does not exist

    # path the aes preprocessed file for this experiment
    #
    currentAesPPFile = aesPPTemplate %  expID
    fpCurrentPP = open(currentAesPPFile, 'w')
    ppList = [] # so we can wee out dupes before writing to the PP file
    
    # process the header line
    #
    headerList = str.split(fpAes.readline(), TAB)
    if headerList == ['']: # means file is empty
        return 2
    # find the idx of the columns we want - they are not ordered
    #
    sourceSampleIDX = None
    enaSampleIDX = None
    runIDX = None
    for idx, colName in enumerate(headerList):
        colName = str.strip(colName)
        if str.find(colName, 'Source Name') != -1:
            sourceSampleIDX = idx
        elif str.find(colName, 'ENA_SAMPLE') != -1:
            enaSampleIDX = idx
        elif str.find(colName, 'ENA_RUN') != -1:
            runIDX = idx
    
    # now iterate through each gene/run/tpm in the file
    #
    for line in fpAes.readlines():
        if line == '\n':
            continue 
        tokens = str.split(line, TAB)
        sampleID = '' # The sampleID we load to the lookup
        ssID = '' # ID from 'Source Name' column
        enaID = '' # ID from ENA_SAMPLE column 
        runID = '' # ID from ENA_RUN column
        if runIDX == None:
            # report that the READ column cannot be found
            #
            noRunOrSampleColList.append('Exp: %s file has no Run column' % expID)
            return 3

        # get the sampleID from one or both of the fields
        #
        if sourceSampleIDX != None:
            ssID = str.strip(tokens[sourceSampleIDX])

        if enaSampleIDX != None:
            enaID = str.strip(tokens[enaSampleIDX])

        # report and return if we did not find either sample column
        #
        if not (ssID or enaID): 
            noRunOrSampleColList.append('Exp: %s file has neither Sample column' % expID)
            return 4
        # choose the ID to use: not empty and in database. Checking the database
        # is the only way we know if we have a sampleID (and not something else)
        #
        inDb = 0 # sample ID in the database? 0/1
        if ssID and ssID in sampleInDbDict:
            inDb = 1
            sampleID = ssID
        elif enaID and enaID in sampleInDbDict:
            inDb = 1
            sampleID = enaID

        if inDb == 0:
            # report that the ID is not in the database
            sampleIdNotInDbList.append('Exp: %s Sample Name ID %s and Sample ENA ID %s not in the database' % (expID, ssID, enaID))
        elif len(sampleInDbDict[sampleID]) > 1: 
            # report ambiguous expID/sampleID in database
            ambiguousSampleInDbList.append('Exp: %s Sample ID: %s' % (expID, sampleID))
            return 5
        else:
            runID = str.strip(tokens[runIDX])
            aesRunIdSet.add(runID)
            line = '%s%s%s%s' % (runID, TAB, sampleID, CRT)
            if line not in ppList:
                ppList.append(line)
    # we write all samples except those not in the database to the file. 
    # We will exclude non-relevant and J:DO strain samples later so we may 
    # differentiate btwn a)non-relevant b)J:DO c)not in the aes file 
    # (excluding those not in the database
    #
    for line in ppList:
        fpCurrentPP.write(line)
    fpCurrentPP.close();

    return 0

# end ppAESFile ()---------------------------------------------------

#
# Purpose: Removes the first token of a set of tokens. Used to strip the
#	Gene ID column from the header and the gene lines so the Run IDs 
#	can be more easily processed
# Returns: incoming tokens with first one removed
# 

def splitOffGene(tokens):
    del tokens[0] # remove 'Gene ID' 

    return tokens

# end splitOffGene ()--------------------------------------------

#
# Purpose: writes all QC in various data structures to the curation log with proper
#       header and total count
# Returns: 0
# Assumes: Nothing
# Effects: writes to a file in filesystem
# Throws: Nothing
#

def writeQC():
    fpCur.write("%sExperimentIDs From Experiment Input File Not in the Database%s" % (CRT, CRT))
    fpCur.write('----------------------------------------------------%s' % CRT)
    for e in experimentNotInDbList:
        fpCur.write('%s%s' % (e, CRT))
    fpCur.write('Total: %s%s' % (len(experimentNotInDbList), CRT))

    fpCur.write('%sExperiments missing Run Column or Sample Columns, or both%s' % (CRT, CRT))
    fpCur.write('--------------------------------------------%s' % CRT)
    for e in noRunOrSampleColList:
        fpCur.write('%s%s' % (e, CRT))
    fpCur.write('Total: %s%s' % (len(noRunOrSampleColList), CRT))

    fpCur.write('%sSample IDs from AE file not in the Database%s' % (CRT, CRT))
    fpCur.write('---------------------------------------------%s' % CRT)
    for e in sampleIdNotInDbList:
        fpCur.write('%s%s' % (e, CRT))
    fpCur.write('Total: %s%s' % (len(sampleIdNotInDbList), CRT))

    fpCur.write('%sEnsembl ID in EAE file is associated with > 1 marker in MGI%s' % (CRT, CRT))
    fpCur.write('-------------------------------------------------%s' % CRT)
    for e in multiMarkerEnsemblList:
        fpCur.write('%s%s' % (e, CRT))
    fpCur.write('Total: %s%s' % (len(multiMarkerEnsemblList), CRT))

    fpCur.write('%sMultiple Ensembl IDs in EAE file are associated with the same marker in MGI%s' % (CRT, CRT))
    fpCur.write('-------------------------------------------------%s' % CRT)
    for e in multiEnsemblMarkerList:
        fpCur.write('%s%s' % (e, CRT))
    fpCur.write('Total: %s%s' % (len(multiEnsemblMarkerList), CRT))

    fpCur.write('%sEnsembl ID in EAE file is not in MGI%s' % (CRT, CRT))
    fpCur.write('-------------------------------------------------%s' % CRT)
    for e in ensemblNotInMGIList:
        fpCur.write('%s%s' % (e, CRT))
    fpCur.write('Total: %s%s' % (len(ensemblNotInMGIList), CRT))

    fpCur.write('%sEnsembl ID in EAE file is an orphan in MGI, only has sequence association%s' % (CRT, CRT))
    fpCur.write('-------------------------------------------------%s' % CRT)
    for e in ensemblIsOrphanList:
        fpCur.write('%s%s' % (e, CRT))
    fpCur.write('Total: %s%s' % (len(ensemblIsOrphanList), CRT))

    fpCur.write('%sAmbiguous expID/sampleID in database%s' % (CRT, CRT))
    fpCur.write('-------------------------------------------------%s' % CRT)
    for e in ambiguousSampleInDbList:
        fpCur.write('%s%s' % (e, CRT))
    fpCur.write('Total: %s%s' % (len(ambiguousSampleInDbList), CRT))
    
    fpCur.write('%sRuns with empty TPM value set to 0.0%s' % (CRT, CRT))
    fpCur.write('-------------------------------------------------%s' % CRT)
    for e in noTpmValueList:
        fpCur.write('%s%s' % (e, CRT))
    fpCur.write('Total: %s%s' % (len(noTpmValueList), CRT))

    fpCur.write('%sRun IDs not in ArrayExpress File%s' % (CRT, CRT))
    fpCur.write('-------------------------%s' % CRT)
    for e in runIdNotInAEList:
        fpCur.write('%s%s' % (e, CRT))
    fpCur.write('Total: %s%s' % (len(runIdNotInAEList), CRT))

    fpCur.write('%sRun IDs not in Expression Atlas File%s' % (CRT, CRT))
    fpCur.write('-------------------------%s' % CRT)
    for e in runIdNotInEAList:
        fpCur.write('%s%s' % (e, CRT))
    fpCur.write('Total: %s%s' % (len(runIdNotInEAList), CRT))

    fpCur.write('%sSamples not flagged as Relevant in the Database%s' % (CRT, CRT))
    fpCur.write('-------------------------------------------------%s' % CRT)
    for e in nonRelSkippedList:
        fpCur.write('%s%s' % (e, CRT))
    fpCur.write('Total: %s%s' % (len(nonRelSkippedList), CRT))

    fpCur.write('%sSamples with strain: J:DO%s' % (CRT, CRT))
    fpCur.write('-------------------------------------------------%s' % CRT)
    for e in JDOSkippedList:
        fpCur.write('%s%s' % (e, CRT))
    fpCur.write('Total: %s%s' % (len(JDOSkippedList), CRT))

    fpCur.write('%sSamples where average of the average TPM SD across all genes is > %s%s' % (CRT, stdDevCutoff, CRT))
    fpCur.write('-------------------------------------------------%s' % CRT)
    for e in highAveStdDevList:
        fpCur.write('%s%s' % (e, CRT))
    fpCur.write('Total: %s%s' % (len(highAveStdDevList), CRT))

    return 0

# end writeQC ()--------------------------------------------

#
# Purpose: creates file descriptors and reads the Expression Atlas file pulling out
#       the Gene ID, Run ID and TPM. Writes to the EAE preprocessed intermediate file
# Returns: 0 if successful, 1 if EAE file does not exist, 2 if AES file is empty
#       3 if and additional header is found within the file. This indicates a bad download
#	Once we decouple and augment the download of these files, this error should not
#	happen
# Assumes: Nothing
# Effects: creates file in filesystem
# Throws: Nothing
#

def ppEAEFile(expID):

    global noTpmValueList, multiMarkerEnsemblList, currentEaePPFile
    global  ensemblNotInMGIList, ensemblIsOrphanList, eaeRunIdSet
    global multiEnsemblMarkerList

    #print 'in ppEAEFile(expID): %s' % expID
    start_time = time.time()

    # we create a set here, because eaeRunIdList will have NL char
    # after last run ID, need this later for set difference
    #
    eaeRunIdSet = set()

    #  create eae input file descriptor for this experiment
    #
    eaeFile = eaeTemplate % expID
    #print 'eaeFile: %s' % eaeFile
    try:
        fpEae = open(eaeFile, 'r')
    except:
        return 1 # file does not exist

    # path the aes preprocessed file for this experiment
    #
    currentEaePPFile = eaePPTemplate % expID
    fpCurrentPP = open(currentEaePPFile, 'w')

    # get the runIDs from the header, we will use the index of te`e tpms
    # in each gene/tpm line to get the runID
    #
    headerList = str.split(fpEae.readline(), TAB)
    #print 'headerList: %s' % headerList
    if headerList == ['']: # means file is empty
        return 2
    eaeRunIdList = splitOffGene(headerList)
    #print 'Num runIDs for %s: %s %s' % (expID, len(eaeRunIdList), eaeRunIdList)
    sys.stdout.flush()

    # process each gene and it's sample tpm's
    #
    for line in fpEae.readlines():
        tokens = str.split(line, TAB)
        geneID = str.strip(tokens[0])
        
        # multi marker per ensembl
        if geneID in multiMarkerEnsemblDict:
            symbols = ', '.join(multiMarkerEnsemblDict[geneID])
            msg = '%s: %s' % (geneID, symbols)
            if msg not in multiMarkerEnsemblList:
                multiMarkerEnsemblList.append(msg)
            continue
        
        # multi ensembl per marker
        elif geneID in multiEnsemblMarkerDict:
            symbol = multiEnsemblMarkerDict[geneID]
            #if len(symbolToMultiEnsIdDict[symbol]) > 1:
            ensIDs = ', '.join(symbolToMultiEnsIdDict[symbol])
            msg = '%s: %s' % (symbol, ensIDs)
            if msg not in multiEnsemblMarkerList:
                multiEnsemblMarkerList.append(msg)
            continue

        # not in MGI
        elif geneID not in ensemblMarkerSet and geneID not in ensemblSequenceSet:
            msg = '%s' % (geneID)
            if msg not in ensemblNotInMGIList:
                ensemblNotInMGIList.append(msg)
            continue

        # assoc only with a sequence
        elif geneID not in ensemblMarkerSet and geneID in ensemblSequenceSet:
            msg = '%s' % (geneID)
            if msg not in ensemblIsOrphanList:
                ensemblIsOrphanList.append(msg)
            continue

        # list of runs (tpm) for geneID
        tpmList = splitOffGene(tokens)

        for idx, tpm in enumerate(tpmList):
            # get the runID that goes with this tpm value
            try:
                runID = str.strip(eaeRunIdList[idx])
            except:
                # this shouldn't happen once we finish the file download script
                print('Multi header record. idx: %s tpm: %s' % (idx, tpm))
                return 3
            eaeRunIdSet.add(runID)
            tpm = str.strip(tpm)

            # report blank tpm values and set them to 0.0
            if tpm == '':
                tpm = 0.0
                noTpmValueList.append('ExpID: %s geneID: %s runID: %s' % \
                    (expID, geneID, runID))

            fpCurrentPP.write('%s%s%s%s%s%s' % (geneID, TAB, runID, TAB, tpm, CRT))

    fpCurrentPP.close();
    elapsed_time = time.time() - start_time
    print('%sTIME: Processing Runs from the EAE file %s %s%s' % (CRT, expID, time.strftime("%H:%M:%S", time.gmtime(elapsed_time)), CRT)) 

    return 0

# end ppEAEFile ()--------------------------------------------

#
# Purpose: joins the preprocessed ArrayExpress and Expression Atlas file 
#	on runID. This is how we get the sampleID for each run
#	calls the product script 'join_files' which uses unix join
# Returns: 0 if successful
# Assumes: Nothing
# Effects: creates file in filesystem
# Throws: Nothing
#

def createJoinedFile(joinedFile):

    start_time = time.time()
    cmd = "%s/run_join %s %s %s" % (binDir, currentEaePPFile, currentAesPPFile, joinedFile)
    #rc = os.system(cmd)
    #if rc != 0:
    #    msg = 'join cmd did not succeed: %s%s' % (cmd, CRT)
    #    fpDiag.write(msg)
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    stdout = result.stdout
    stderr = result.stderr
    statusCode = result.returncode
    if statusCode != 0:
        msg = 'Joining %s and %s failed statusCode: %s stderr: %s%s' % (currentEaePPFile, currentAesPPFile, statusCode, stderr, CRT)
        fpLog.write(msg)
        return statusCode
    else:
        elapsed_time = time.time() - start_time
        print('%sTIME: Creating the join file %s %s%s' % (CRT, joinedFile, time.strftime("%H:%M:%S", time.gmtime(elapsed_time)), CRT))
        sys.stdout.flush()

        return 0

# end createJoinedFile ()--------------------------------------------

#
# Purpose: processes the joined file, doing QC and creating the geneDict 
# dictionary for those records that pass QC
# 'geneDict' maps the ensemblGeneID to a dictionary of sampleIDs mapped to
#	their raw TPM values
# Returns: 0 if successful
# Assumes: Nothing
# Effects: reads a file in the file system
# Throws: Nothing
#

def processJoinedFile(expID, joinedFile):
    global runIdNotInAEList, nonRelSkippedList, JDOSkippedList
    print('joinedFile: %s' % joinedFile)
    start_time = time.time()
    # {geneID: {sampleID:[tpm1, ...], ...}, ...}
    geneDict = {}

    fpJoined = open(joinedFile, 'r')
    line = fpJoined.readline()
    
    while line:
        tokens = str.split(line, TAB)
        geneID = tokens[0]
        runID = tokens[1]
        tpm = tokens[2]
        sampleID = ''

        # some join files have no column 4, do a try/except
        try:
            sampleID = str.strip(tokens[3])
        except:
            sampleID = ''
        # some join files have empty column 4
        if sampleID == '':
            msg = '%s: %s' % (expID, runID)
            if msg not in  runIdNotInAEList:
                runIdNotInAEList.append(msg)
            line = fpJoined.readline()
            continue
        # is sampleID flagged as relevant in the db? If not report/skip
        # we do this here rather than excluding when create the aes pp file
        # so we may differentiate between sample id 1) not in database and
        # 2) a) not relevant b) in J:DO strain set
        lookupVal = '%s|%s' % (expID, sampleID)
        if lookupVal not in relevantSampleSet:
            msg = '%s: %s' % (expID, sampleID)
            if msg not in  nonRelSkippedList:
                nonRelSkippedList.append(msg)
            line = fpJoined.readline()
            continue
        elif sampleID in JDOSampleSet:
            # Report/skip if the sampleID is in the J:DO strain set
            msg = '%s: %s' % (expID, sampleID)
            if msg not in JDOSkippedList:
                JDOSkippedList.append(msg)
            line = fpJoined.readline()
            continue

        if geneID not in geneDict:
            geneDict[geneID] = {}
        if sampleID not in geneDict[geneID]:
            geneDict[geneID][sampleID] = []
        geneDict[geneID][sampleID].append(float(tpm))
        line = fpJoined.readline()

    elapsed_time = time.time() - start_time
    print('%sTIME: Processing the join file %s %s%s' % (CRT, joinedFile, time.strftime("%H:%M:%S", time.gmtime(elapsed_time)), CRT))

    return geneDict

# end processJoinedFile -----------------------------------------------------------

#
# Purpose: calculate the level for the average QN TPM
# Returns: proper level key for aveQnTpm
# Assumes: Nothing
# Effects: 
# Throws: Nothing
#
def calcLevel(aveQnTpm):
    level = None
    if aveQnTpm < 0.5:	
        level = BELOW_CUTOFF 
    elif aveQnTpm >=  0.5 and aveQnTpm <= 10:
        level = LOW
    elif aveQnTpm >=  11 and aveQnTpm <= 1000:
        level = MED
    else:  # aveQnTpm > 1000:
        level = HIGH

    return  level
# end calcLevel ------------------------------------------------------------------
        
#
# Purpose: creates RNASeq and Combined bcp files from a list of matrices (2-dimensional
#       arrays) created by another function
# Returns: 0 if successful
# Assumes: Nothing
# Effects: creates files in filesystem
# Throws: Nothing
#
def writeBCP(expID, matrixList):
    global rnaSeqKey, combinedKey, bcpCommandList

    #print 'creating bcp files for expID: %s' % expID
    start_time =  time.time()

    combinedFile = '%s.%s' % (combinedBcpFile, expID)
    combinedFilePath = '%s/%s' % (outputDir, combinedFile)
    initCombinedBcp(combinedFilePath)
 
    rnaSeqFile = '%s.%s' % (rnaSeqBcpFile, expID)
    rnaSeqFilePath = '%s/%s' % (outputDir, rnaSeqFile)
    initRnaSeqBcp(rnaSeqFilePath)

    # counts in the bcp file
    combinedLineCt = 0
    rnaSeqLineCt = 0 

    for matrix in matrixList:
        for row in matrix:
            #print row
            rowLen = len(row)
            # position 0 is always markerKey
            # position rowLen -1 is # Biological Replicates
            # position rowLen - 2 is ave QN TPM (all samples)
            markerKey =  row[0] 
            numBioRepl = row[-1]
            aveQnTpm = row[-2]
            #print 'avQnTpm: %s' % aveQnTpm
            levelKey = calcLevel(aveQnTpm)
            #print 'levelKey: %s' % levelKey
            line = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (combinedKey, TAB, markerKey, TAB, levelKey, TAB, numBioRepl, TAB, aveQnTpm, TAB, createdByKey, TAB, createdByKey, TAB, loaddate, TAB, loaddate, CRT)
            #print 'combined: %s' %line
            fpCombinedBcp.write(line)
            combinedLineCt += 1

            # positions 1 to rowLen-2 are samples: sampleKey|aveTPM|qnTpm
            # to get samples we iterate over rowLen-3 starting with i+1 
            for i in range(rowLen-3):
                #print 'next sample: %s' % row[i+1]
                sampleKey, aveTpm, qnTpm = str.split(row[i+1], PIPE)
                line = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (rnaSeqKey, TAB, sampleKey, TAB, combinedKey, TAB, markerKey, TAB, aveTpm, TAB, qnTpm, TAB, createdByKey, TAB, createdByKey, TAB, loaddate, TAB, loaddate, CRT)
                #print 'rnaseq: %s' % line
                fpRnaSeqBcp.write(line)
                rnaSeqLineCt += 1
                rnaSeqKey += 1
            combinedKey +=1

    closeCombinedBcp()
    closeRnaSeqBcp()

    # gather bcp commands in a List so we can execute them together at the
    # end of processing
    #
    cmd = '%s %s %s %s %s %s "\\t" "\\n" mgd' % \
    (bcpCommand, db.get_sqlServer(), db.get_sqlDatabase(), combinedTable, outputDir, combinedFile)
    bcpCommandList.append(cmd)
    #print 'Num combined bcp lines for expID %s: %s' % (expID, combinedLineCt)

    cmd = '%s %s %s %s %s %s "\\t" "\\n" mgd' % \
    (bcpCommand, db.get_sqlServer(), db.get_sqlDatabase(), rnaSeqTable, outputDir, rnaSeqFile)
    bcpCommandList.append(cmd)
    #print 'Num rnaSeq bcp lines for expID %s: %s' % (expID, rnaSeqLineCt)

    elapsed_time = time.time() - start_time
    elapsed_time = time.time() - start_time
    print('%sTIME: Creating bcp files %s %s%s' % (CRT, expID, time.strftime("%H:%M:%S", time.gmtime(elapsed_time)), CRT))

    return 0

# end writeBCP ----------------------------------------------------------------

#
# Purpose: calculates the aveTPM and Std Dev and ave Std Dev for all samples of the gene 
# Returns: dictionary of average TPM values (by gene) for each sample of the gene
# This is the format needed for input to the quantile normalize function
# {sampleKey:{markerKey:aveTPM, ...}, sampleKey2:{markerKey:aveTPM, ...}, ...}
# Assumes: Nothing
# Effects: 
# Throws: Nothing
#

def calcTPMAveSD(expID, geneDict):
    global sampleAveSDDict, highAveStdDevList

    start_time =  time.time()

    # {sampleKey:{markerKey:aveTPM, ...}, sampleKey2:{markerKey:aveTPM, ...}, ...}
    # input to QN
    aveTPMDict = {}

    #print 'calculating TPM AVE and SD for expID: %s' % expID
    for geneID in geneDict:
        sampleDict = geneDict[geneID]
        for sampleID in sampleDict:
            # dupe samples for current expID have already been QC, if exist this expID is
            # skipped. we can assume only one sampleID in the list
            sampleKey = sampleInDbDict[sampleID][0]
            tpmList = sampleDict[sampleID]
            techReplCt = len(tpmList)       # number of runs for the sample (1-n)
            aveTpm = calcAve(tpmList, 2)       # calculate the average tpm
            stdDev = round(np.std(tpmList), 2)  # calc the SD of the tpms
            # calculate the SD average
            stdDevAve = 0.0
            if aveTpm != 0.0:
                stdDevAve = round(stdDev / aveTpm, 2)

            # get the marker key
            markerKey = ensemblMarkerDict[geneID]

            # load the QN input dictionary 
            if sampleKey not in aveTPMDict:
                aveTPMDict[sampleKey] = {}
                #print 'calcTPMAveSD sampleID; %s sampleKey: %s' % (sampleID, sampleKey)
            aveTPMDict[sampleKey][markerKey] = aveTpm

            # collect the stdDevAve of the technical replicates - we will want to
            # take the average of all aveSD across all genes of a sample to
            # determine samples to include/exclude.
            key = '%s%s%s%s%s' % (expID, PIPE, sampleID, PIPE, sampleKey)
            if key not in sampleAveSDDict:
                sampleAveSDDict[key] = []
            sampleAveSDDict[key].append(stdDevAve)

            # this reports the run tpm's, the tpm average, the stdDev of the run
            # tpms and the stdDevAve (stdDev/aveTpm) as well as the count of
            # technical replicates (number of runs per sample)
            #
            fpStudentRpt.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (expID, TAB, geneID, TAB, sampleID, TAB, ', '.join(list(map(str, tpmList))),  TAB, aveTpm, TAB, stdDev, TAB, stdDevAve, TAB, techReplCt, CRT))

    # Now iterate through the sampleAveSDDict and calculate the average of
    # each aveSD. Report any above a threshold weeding out dupes.
    #
    keyList = sorted(sampleAveSDDict.keys())

    for key in keyList:
        eID, sampleID, sampleKey = str.split(key, PIPE)
        aveAveSdAllGenes = calcAve(sampleAveSDDict[key], 2)

        # if average of the average SD across genes in > stdDevCutoff (config)
        # report it and remove all genes for that sample from the aveTPMDict
        if aveAveSdAllGenes > stdDevCutoff:
            del aveTPMDict[sampleKey]
            reportLine = '%s%s%s%s%s%s' % \
                (eID, TAB, sampleID, TAB, aveAveSdAllGenes, CRT)
            if reportLine not in highAveStdDevList:
                highAveStdDevList.append(reportLine)

    elapsed_time = time.time() - start_time
    print('%sTIME: calcTPMAveSD %s %s%s' % (CRT, expID, time.strftime("%H:%M:%S", time.gmtime(elapsed_time)), CRT))

    return  aveTPMDict

# end calcTPMAveSD -------------------------------------------------------------

#
# Purpose: executes all bcp commands for RNASeq tables and sets their primary 
#	key sequence
# Returns: 0 if successful
# Assumes: Nothing
# Effects: creates records in the database, writes to the diagnostic log
# Throws: Nothing
#

def execBCP ():

    # execute all the bcp files
    for bcpCmd in bcpCommandList:
        fpDiag.write('%s\n' % bcpCmd)
        os.system(bcpCmd)

    # reset the rnaseq primary key sequence
    db.sql('''select setval('gxd_htsample_rnaseq_seq', (select max(_rnaseq_key) from gxd_htsample_rnaseq))''', None)

    # reset the rnaseq combined primary key sequence
    db.sql('''select setval('gxd_htsample_rnaseqcombined_seq', (select max(_rnaseqcombined_key) from gxd_htsample_rnaseqcombined))''', None)

    db.commit()

    return 0

# end execBCP ------------------------------------------------------------------

#
# Purpose: Query the RNASeq set/setmembers and find biological replicates
# 	and their list of samples
# Returns: dictionary mapping biological replicate key to its set of samples
# Assumes: Nothing
# Effects:
# Throws: Nothing
#

def getBioReplicates(expID):
    expKey = experimentInDbDict[expID]
    results = db.sql('''select s.*, sm._Sample_key 
        from GXD_HTSample_RNASeqSet s,
        GXD_HTSample_RNASeqSetMember sm
        where s._experiment_key = %s
        and s._rnaSeqSet_key = sm._rnaSeqSet_key''' % expKey, 'auto')
    #print 'len bio replicates results: %s' % len(results)

    # {attributeKey:set(sampleKeys), ...}
    replisetDict = {}

    # iterate through results and map the replicate to its set of sample keys
    for r in results:
        sampleKey = r['_sample_key']
        expKey = r['_experiment_key']
        age = r['age']
        orgKey = r['_organism_key']
        sexKey = r['_sex_key']
        emapaKey = r['_emapa_key']
        stageKey = r['_stage_key']
        genotypeKey = r['_genotype_key']
        note = r['note']
        if note == None:
            note = ''

        key = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (expKey, PIPE, age, PIPE, orgKey, PIPE, sexKey, PIPE, emapaKey, PIPE, stageKey, PIPE, genotypeKey, PIPE, note)
        if key not in replisetDict:
            replisetDict[key] = set()
        replisetDict[key].add(sampleKey)
    #print replisetDict 
    return replisetDict
# end getBioReplicates ------------------------------------------------------------------

#
# Purpose: main processing function. Gets the list of experiment IDs from the
#	'RNA Seq Load Experiment' set and processes experiments one by one. Does some QC
# Returns: 0 if successful; handles return codes from all functions it calls 
# Assumes: Nothing
# Effects:
# Throws: Nothing
#

def process():
    global  experimentNotInDbList, highAveStdDevList

    start_time = time.time()
    
    # write the header to the student report up front
    fpStudentRpt.write('expID%sgeneID%ssampleID%stechRepl%saveTpm%sstdDev%sstdDevAve%stechRepCt%s' % (TAB, TAB, TAB, TAB, TAB, TAB, TAB, CRT))

    #
    # for each expID in  the 'RNA Seq Load Experiment' MGI_Set:
    #
    for r in rnaSeqSetResults:
        expID = str.strip(r['accid'])

        # load the list of sampleIDS and their keys from the database 
        # for current experiment there can be > 1
        loadSampleInDbDict(expID)

        # the list of matrices for this experiment, one per repliset
        matrixList = []

        # report if expID not in the database and skip
        #
        if expID not in experimentInDbDict:
            experimentNotInDbList.append(expID)
            continue
        # preprocess the aes file for this expID (run, sample)
        #
        sys.stdout.flush()
        rc = ppAESFile(expID)
        if rc != 0:
            print('''preprocessing AES file returned rc %s, 
                        skipping file for %s''' % (rc, expID))
            continue

        # preprocess the eae file for this expID (gene, run, tpm)
        #
        sys.stdout.flush()
        rc = ppEAEFile(expID)
        if rc != 0:
            print('''preprocessing EAE file returned rc %s, 
                    skipping file for %s''' % (rc, expID))
            continue

        # check for AES runIDs not in the EAE file
        #
        diff = aesRunIdSet.difference(eaeRunIdSet)
        if len(diff):
            diffString = ', '.join(diff)
            runIdNotInEAList.append('%s: %s' % (expID, diffString))

        # Create the joined file from the two preprocessed (pp) files
        #
        joinedFile =  joinedPPTemplate % expID
        cjf_rc = createJoinedFile(joinedFile)
        if cjf_rc != 0:
            return 1 

        # process the joined file, creating dict by geneID
        # {geneID: {sampleID:[tpm1, ...], ...}, ...}
        #
        geneDict = processJoinedFile(expID, joinedFile)
        
        # {sampleKey:{markerKey:aveTPM, ...}, 
        #	sampleKey2:{markerKey:aveTPM, ...}, ...}
        aveTPMDict = calcTPMAveSD(expID, geneDict)
        #print 'aveTPMDict keys: %s' % aveTPMDict.keys()
        if not aveTPMDict: # nothing to normalize
            continue

        # {attributeKey:set(sampleKeys), ...}
        replisetDict = getBioReplicates(expID)
        # number of genes
        numRows = 0 # set this later when we know!

        # Preprocess the replisetDict:
        # a) remove samples in the db not in input
        # b) # samples is  used to determine # columns in the matrix
        for key in replisetDict:
            sampleSet = replisetDict[key]
            newSet = set()	
            for sampleKey in sampleSet:
                # if we find samples in DB not in aveTPMDict we need to 
                # remove from the set
                # remember - QC can kick out some samples
                # qnOutputDict is produced from aveTPMDict so don't need
                # to check it too
                if sampleKey not in aveTPMDict:
                    #print 'sampleKey: %s for expID: %s not in aveTPMDict' % \
                    #    (sampleKey, expID)
                    #sampleSet.remove(sampleKey) # can't do this get 
                        # RuntimeError: Set changed size during iteration
                        # Create a new set
                    #print 'removing %s from sampleSet' % sampleKey
                    continue
                newSet.add(sampleKey)
                # set this only once, all samples have same number of genes
                if numRows == 0: 
                    numRows = len(aveTPMDict[sampleKey])
            # in case we changed the sampleSet, reset the new set in the Dict
            #print 'Adding newSet %s to replisetDict key: %s' % (newSet, key)
            replisetDict[key] = newSet

        # now we have removed any kicked out samples from replisetDict
        # recall 'key' is our compound key of all attributes
        for key in replisetDict:
            #print 'process replisetDict key: %s' % key
            # gather all the numBioReplicates by rowNum
            #  {rowNum:numBioReplicates, ...}
            numBioReplDict = {}
           
            sampleSet = replisetDict[key]

            # if the sample set is empty then we skip this repliset
            if not sampleSet:
                #print 'no samples for this repliset'
                continue

            #print 'process replisetDict sampleSet: %s' % sampleSet
            totalSamples = len(sampleSet)

            #
            # get the aveTPMs for this repliset so we may QN them
            #
            replisetAveTpmDict = {}
            for sampleKey in sampleSet:
                replisetAveTpmDict[sampleKey] = aveTPMDict[sampleKey]
            # Now QN them
            qnInput =  pd.DataFrame(replisetAveTpmDict)
            qnOutput = quantileNormalize.qn(qnInput, list(replisetAveTpmDict.keys()))

            # {sampleKey:{markerKey:qnAveTPM, ...},
            #       sampleKey2:{markerKey:qnAveTPM, ...}, ...}
            qnOutputDict = qnOutput.to_dict()

            #print 'process replisetDict key: %s samples: %s numBioReplicates: %s' % (key, sampleSet, len(sampleSet))
            
            # gather all the QN TPMs for averaging across all samples of a gene
            # rowNum = gene, we use number so we can sort the dict keys to 
            # assign the proper aveQNTPM to the proper gene (note python 3.* has
            # 'OrderedDict'
            # {rowNum:[list of qnTPM for the gene], ...}
            aveQNTpmDict = {}

            # 1  + number of samples  + 2
            # gene + total samples + qn tpm across all samples +
            # num bio replicates
            numColumns = 1 + totalSamples + 2

            # now our sampleSet, numColumns, numRows are set and correct
            #  create an empty matrix with str.data type
            matrix = np.empty((numRows, numColumns), dtype='object')

            # now we pull everything into a 2-D matrix from aveTPMDict 
            # and qnOutputDict for each biological replicate set for easy bcp 
            # file creation
            currentColumnNum = 1    # The first sample column
            currentRowNum = 0       # the first row
            
            for sampleKey in sampleSet:
                #print 'next sampleKey: %s' % sampleKey
                try:
                    # {markerKey:aveTPM, ...}
                    aveTpmByGeneDict = aveTPMDict[sampleKey]
                except:
                    # this shouldn't happend because we synced things up
                    # above
                    'sampleKey %s not in aveTPMDict' % sampleKey

                # {markerKey:qnAveTPM, ...}
                qnTpmByGeneDict  = qnOutputDict[sampleKey]
        
                for mKey in aveTpmByGeneDict:
                    # add number of replicates to dict by row number (gene )
                    numBioReplDict[currentRowNum] = totalSamples

                    aveTPM = aveTpmByGeneDict[mKey]
                    qnTPM = round(qnTpmByGeneDict[mKey], 2)
                    sampleValues = '%s%s%s%s%s' % (sampleKey, PIPE, aveTPM, PIPE, qnTPM)
                    # create/add to the list of qnTPM for the gene that will
                    # will be processed later after we have built up the matrix
                    # with its genes and sample aveTPM and qnTPM
                    if currentRowNum not in aveQNTpmDict:
                        aveQNTpmDict[currentRowNum] = []
                    aveQNTpmDict[currentRowNum].append(qnTPM)
                    matrix[currentRowNum][0] = mKey
                    matrix[currentRowNum][currentColumnNum] = sampleValues
                    #print 'matrix[%s][%s]: %s' % (currentRowNum, currentColumnNum, matrix[currentRowNum][currentColumnNum])
                    currentRowNum += 1
                currentRowNum = 0
                currentColumnNum += 1

            # add the ave QN to the matrix for each gene
            for rowNum in aveQNTpmDict:
                ave = calcAve(aveQNTpmDict[rowNum], 1)
                if ave >= 1:
                    ave = int(round(ave, 0)) # if >1 round then truncate decimal
                matrix[rowNum][numColumns - 2] = ave

            # add the number of bio  replicates to matrix for each gene
            for rowNum in numBioReplDict:
                matrix[rowNum][numColumns - 1] = numBioReplDict[rowNum]
            # Now matrix is complete
            matrixList.append(matrix)

            #print 'Row 1 of matrix:'
            #print  matrix[0]

            #print '\nRow 2 of matrix'
            #print matrix[1]
            
            #print '\nRow 3 of matrix'
            #print matrix[2]

        # write out bcp for all replisets of this experiment
        writeBCP(expID, matrixList)

    # we've processed all experiments, now execute bcp
    execBCP()
    return 0

# end process ()--------------------------------------------

#
# Purpose: closes all file descriptors that remain open 
# 	across all experiments
#

def closefiles():
    fpCur.close()
    fpDiag.close()
    fpStudentRpt.close()
    db.useOneConnection(0)

    return 0

# end closefiles ()--------------------------------------------

#
# Main
#

# -------------------------------------------------------------
START_TIME = time.time()

print('Start time: %s' %  mgi_utils.date())
sys.stdout.flush()
if init() != 0:
     exit(1, 'Error in  init \n' )
elapsed_time = time.time() - START_TIME
print('TIME to run init function %s' %  time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
sys.stdout.flush()

# -------------------------------------------------------------
TIME = time.time()

if process() != 0:
     exit(1, 'Error in  process \n' )
elapsed_time = time.time() - TIME
print('TIME to run process function %s' % time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
sys.stdout.flush()

# -------------------------------------------------------------
TIME = time.time()

if writeQC() != 0:
     exit(1, 'Error in writeQC \n' )
elapsed_time = time.time() - TIME
print('TIME to run writeQC function %s' % time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
sys.stdout.flush()

# -------------------------------------------------------------

closefiles()

elapsed_time = time.time() - START_TIME

print('End time: %s'  % mgi_utils.date())

print('Total run time: %s' %  time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
