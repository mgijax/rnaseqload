##########################################################################
#
# Purpose: For each experiment ID in the 'Baseline RNASeq Load Experiment' MGI_Set
#	Create intermediate files for EAE and AES files (pulls out just the data needed)
#	Join the two sets of intermediate files
#	Determine biological replicates
#	Run Quantile Normalization
#	Create bcp files
#	Execute bcp
#  	Write QC reports
#
# Usage: rnaseqload.py
#
# Env Vars:
#	 1.  LOGDIR 
#	 2a. BASELINEINPUTDIR - intermediate files generated from RAW input files
#	 2b. BASELINERAW_INPUTDIR - files downloaded from source
#	 3.  BASELINEOUTPUTDIR - rnaseq bcp files
#	 4.  INSTALLDIR - for path to run_join script	 
#	 5.  LOG_CUR - curation log
#	 6.  LOG_DIAG - diagnostic log
#	 7a. RNASEQ_BCP - bcp filename suffix for rnaseq table - we will append expID
#	 7b. COMBINED_BCP - bcp filename suffix for combined table - we will append expID
#	 8.  PG_DBUTILS - for path to bcpin.csh script
#	 9.  AES_LOCAL_FILE_TEMPLATE - path and template for processed aes files
#	 10. AES_PP_FILE_TEMPLATE - preprocessed to just runID, sampleID
#	 11. EAE_TPMS_PP_FILE_TEMPLATE - path and template for processed eae files
#    12. EAE_GROUP_PP_FILE_TEMPLATE - path and template for processed eae files
#    13. AES_SDRF_PP_FILE_TEMPLATE - path and template for processed eae files
#
# Inputs:
#	1. Database: Baseline RNASeq Experiment set
#	2. ArrayExpress files by experiment
#	3. Expression Atlas files by experiment
#	4. Configuration (see rnaseqload.config)
#
# Outputs:
#	 1. aes and eae files from source
#	 2. preprocessed file for each experiment aes and eae
#	 3. bcp files - one per experiment
#	 4. curator and diagnostic log
# 
###########################################################################

import os	 	# for system to execute bcp, getenv
import sys 	 	# to flush stdout
import time	 	# used for its time.time() function (for timestamps)
import xml.etree.ElementTree as ET
import db
import mgi_utils 	# for log start/end timestamp
import loadlib	 	# for bcp friendly date/timestamp

# paths to input and two output files
logDir =  os.getenv('LOGDIR')
inputDir =  os.getenv('BASELINEINPUTDIR')
rawInputDir = os.getenv('BASELINERAW_INPUTDIR')
outputDir = os.getenv('BASELINEOUTPUTDIR')

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

# Expression Atlas Experiment file Template - name of file stored locally
tpmsTemplate = '%s' % os.getenv('EAE_TPMS_LOCAL_FILE_TEMPLATE')
tpmsPPTemplate = '%s' % os.getenv('EAE_TPMS_PP_FILE_TEMPLATE')
groupTemplate = '%s' % os.getenv('EAE_GROUP_LOCAL_FILE_TEMPLATE')
groupPPTemplate = '%s' % os.getenv('EAE_GROUP_PP_FILE_TEMPLATE')

# the experients that will be processed 
expIDToProcess = set()

# GXT HT Experiment IDs in the database
experimentInDbDict = {}

# samples with strain J:DO; genotypeKey = 90560 in the db
JDOSampleSet = set()

# samples flagged as relevant in the db
relevantSampleSet = set()

# ensembl IDs assoc w/ =1 marker
ensemblMarkerDict = {}

# marker IDs assoc w/ >1 ensembl
multiMarkerEnsemblDict = {}

# ensembl IDs assoc w/ >1 marker
multiEnsemblMarkerDict = {}

#
# QC data structures
#

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
    global relevantSampleSet
    global ensemblMarkerDict, multiMarkerEnsemblDict, multiEnsemblMarkerDict
    global rnaSeqKey, combinedKey, rnaSeqSetResults

    db.useOneConnection(1)

    # We have truncated the tables
    rnaSeqKey = 1
    combinedKey = 1

    results = db.sql('''
        select accid, _object_key
        from ACC_Accession
        where _MGIType_key = 42
        and preferred = 1
        ''', 'auto')
    for r in results:
        experimentInDbDict[r['accid']] = r['_object_key']

    results =  db.sql('''
        select name
        from GXD_HTSample
        where _Genotype_key = 90560
        ''', 'auto')
    for r in results:
        JDOSampleSet.add(str.strip(r['name']))

    results =  db.sql('''
        select a.accid as exptID, s.name
        from GXD_HTSample s, ACC_Accession a
        where s._Relevance_key = 20475450
        and s._experiment_key = a._object_key
        and a._mgitype_key = 42
        order by a.accid, s.name
        ''', 'auto')
    for r in results:
        relevantSampleSet.add('%s|%s' % (str.strip(r['exptID']), str.strip(r['name'])))

    results = db.sql('''
        WITH ensembls AS (
        select accid
        from ACC_Accession
        where _LogicalDB_key =  60
        and _MGIType_key = 2
        and preferred = 1
        group by accid having count(*) = 1
        )
        select a.accid, a._Object_key, m.symbol
        from ensembls e, ACC_Accession a, MRK_Marker m
        where e.accid = a.accid
        and a._MGIType_key = 2
        and a.preferred = 1
        and a._Object_key = m._Marker_key
        ''', 'auto')
    for r in results:
        key = r['accid']
        value = r
        if key not in ensemblMarkerDict:
            ensemblMarkerDict[key] = []
        ensemblMarkerDict[key].append(value)

    results = db.sql('''
        WITH ensembls AS (
        select _Object_key
        from ACC_Accession
        where _LogicalDB_key =  60
        and _MGIType_key = 2
        and preferred = 1
        group by _Object_key having count(*) > 1
        )
        select a.accid, a._Object_key, m.symbol
        from ensembls e, ACC_Accession a, MRK_Marker m
        where e._Object_key = a._Object_key
        and a._LogicalDB_key =  60
        and a._MGIType_key = 2
        and a.preferred = 1
        and a._Object_key = m._Marker_key
        ''', 'auto')
    for r in results:
        key = r['accid']
        value = r['symbol']
        if key not in multiMarkerEnsemblDict:
            multiMarkerEnsemblDict[key] = []
        multiMarkerEnsemblDict[key].append(value)

    results = db.sql('''
        WITH ensembls AS (
        select accid
        from ACC_Accession
        where _LogicalDB_key =  60
        and _MGIType_key = 2
        and preferred = 1
        group by accid having count(*) > 1
        )
        select a.accid, a._Object_key, m.symbol
        from ensembls e, ACC_Accession a, MRK_Marker m
        where e.accid = a.accid
        and a._MGIType_key = 2
        and a.preferred = 1
        and a._Object_key = m._Marker_key
        ''', 'auto')
    for r in results:
        key = r['accid']
        value = r['symbol']
        if key not in multiEnsemblMarkerDict:
            multiEnsemblMarkerDict[key] = []
        multiEnsemblMarkerDict[key].append(value)

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
# Purpose: creates EAE_GROUP_PP_FILE_TEMPLATE file in BASELINEINPUTDIR folder for given expID
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

    global multiMarkerEnsemblList
    global ensemblNotInMGIList, ensemblIsOrphanList
    global multiEnsemblMarkerList

    print('in ppEAETpmsFile(expID): %s' % expID)
    start_time = time.time()

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

        tokens = str.split(line, TAB)

        if tokens[0] == 'GeneID':
            continue

        ensemblID = str.strip(tokens[0])
        markerSymbol = str.strip(tokens[1])
        g1All = str.strip(tokens[2])
        g2All = str.strip(tokens[3])
        
        # multi marker per ensembl
        if ensemblID in multiMarkerEnsemblDict:
            symbols = multiMarkerEnsemblDict[ensemblID]
            msg = '%s: %s' % (ensemblID, symbols)
            if msg not in multiMarkerEnsemblList:
                multiMarkerEnsemblList.append(msg)
            continue
        
        # multi ensembl per marker
        elif ensemblID in multiEnsemblMarkerDict:
            symbols = multiEnsemblMarkerDict[ensemblID]
            msg = '%s: %s' % (ensemblID, symbols)
            if msg not in multiEnsemblMarkerList:
                multiEnsemblMarkerList.append(msg)
            continue

        # ensembl/marker 1:1 does not exist (should be in one of the other sets above)
        elif ensemblID not in ensemblMarkerDict:
            msg = '%s' % (ensemblID)
            if msg not in ensemblNotInMGIList:
                ensemblNotInMGIList.append(msg)
            continue

        #print(tokens)
        tokensG1 = g1All.split(',')
        tokensG2 = g2All.split(',')
        g1 = tokensG1[2]
        g2 = tokensG2[2]
        markerKey = ensemblMarkerDict[ensemblID][0]['_object_key']
        #print(markerKey, g1, g2)

        # write to the fpPP file
        fpPP.write('%s\t%s\t%s\t%s\t%s\n' % (ensemblID, markerKey, markerSymbol, g1, g2))
        expIDToProcess.add(expID)

    fpEae.close();
    fpPP.close();
    elapsed_time = time.time() - start_time
    print('%sTIME: Processing Runs from the EAE file %s %s%s' % (CRT, expID, time.strftime("%H:%M:%S", time.gmtime(elapsed_time)), CRT)) 

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
    start_time = time.time()

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
            print(child.tag, child.text)
            runids.append(child.text)

        # write to the fpPP output file
        #print('%s\t%s\t%s\n' % (id, label, runids))
        fpPP.write('%s\t%s\t%s\n' % (id, label, '|'.join(runids)))

    fpPP.close();
    elapsed_time = time.time() - start_time
    print('%sTIME: Processing Runs from the EAE file %s %s%s' % (CRT, expID, time.strftime("%H:%M:%S", time.gmtime(elapsed_time)), CRT)) 

    return 0

# end ppEAEGroupFile ()--------------------------------------------
 
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

    print('creating bcp files for expID: %s' % expID)
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

    # gather bcp commands in a List so we can execute them together at the
    # end of processing
    #
    cmd = '%s %s %s %s %s %s "\\t" "\\n" mgd' % \
    (bcpCommand, db.get_sqlServer(), db.get_sqlDatabase(), combinedTable, outputDir, combinedFile)
    bcpCommandList.append(cmd)
    #print('Num combined bcp lines for expID %s: %s' % (expID, combinedLineCt))

    cmd = '%s %s %s %s %s %s "\\t" "\\n" mgd' % \
    (bcpCommand, db.get_sqlServer(), db.get_sqlDatabase(), rnaSeqTable, outputDir, rnaSeqFile)
    bcpCommandList.append(cmd)
    #print('Num rnaSeq bcp lines for expID %s: %s' % (expID, rnaSeqLineCt))

    elapsed_time = time.time() - start_time
    elapsed_time = time.time() - start_time
    print('%sTIME: Creating bcp files %s %s%s' % (CRT, expID, time.strftime("%H:%M:%S", time.gmtime(elapsed_time)), CRT))

    return 0

# end writeBCP ----------------------------------------------------------------

#
# Purpose: executes all bcp commands for RNASeq tables and sets their primary key sequence
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

# Purpose: main processing function. 
#   Gets the list of experiment IDs from the 'Baseline RNASeq Load Experiment' set 
#   Processes experiments one by one
#   Does some QC
# Returns: 0 if successful; handles return codes from all functions it calls 
# Assumes: Nothing
# Effects:
# Throws: Nothing
#

def process():
    global experimentNotInDbList, highAveStdDevList

    start_time = time.time()
    
    # write the header to the student report up front
    fpStudentRpt.write('expID%sgeneID%ssampleID%stechRepl%saveTpm%sstdDev%sstdDevAve%stechRepCt%s' % (TAB, TAB, TAB, TAB, TAB, TAB, TAB, CRT))

    #
    # for each expID in the MGI_Set:
    #
    for r in rnaSeqSetResults:
        expID = str.strip(r['accid'])

        # if expID is not in the acc_accession table (via the gxdhtload) then skip
        if expID not in experimentInDbDict:
            experimentNotInDbList.append(expID)
            continue

        # preprocess the eae/tpms file for this expID
        #sys.stdout.flush()
        #rc = ppEAETpmsFile(expID)
        #if rc != 0:
        #    print('''preprocessing EAE tpms file returned rc %s, skipping file for %s''' % (rc, expID))
        #    continue

        #if expID not in expIDToProcess:
        #    print('''ppEAETpmsFile had issues, skipping ppEAEGroupFile file processing for %s''' % (expID))
        #    continue
            
        # preprocess the eae/group file for this expID
        sys.stdout.flush()
        rc = ppEAEGroupFile(expID)
        if rc != 0:
            print('''preprocessing EAE group file returned rc %s, skipping file for %s''' % (rc, expID))
            continue

        if expID not in expIDToProcess:
            print('''ppEAEGroupFile had issues, skipping XXXX file processing for %s''' % (expID))
            continue
            
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
