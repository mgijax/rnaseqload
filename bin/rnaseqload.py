#!/usr/local/bin/python
# /usr/bin/python is 2.6.* or 2.7.* depending on server. It is required for downloadFiles.py
##########################################################################
#
# Purpose:
#       
#
# Usage: rnaseqload.py
# Env Vars:
#	 1. INPUT_FILE - Connie's file of experiment IDs
#	 2. LOGDIR 
#	 3. INPUTDIR - files are downloaded to this directory
#		and generated intermediate files too
#	 4. OUTPUTDIR - rnaseq bcp files
#	 5. INSTALLDIR - for path to run_join script	 
#	 6. AES_LOCAL_FILE_TEMPLATE - path and template for downloaded aes files
#	 7. AES_PP_FILE_TEMPLATE - preprocessed to just runID, sampleID
#	 8. EAE_LOCAL_FILE_TEMPLATE - path and template for downloaded eae files
#	 9. EAE_PP_FILE_TEMPLATE - preprocessed to just geneID, runID, TPM
#	 10. JOINED_PP_FILE_TEMPLATE - path & template for joined aes/eae files
# Inputs:
#	1. INPUTFILE - Connie's file of experiment IDs
#	2. ArrayExpress files by experiment
#	3. Expression Atlas files by experiment
#	4. Configuration (see rnaseqload.config)
#
# Outputs:
#	 1.  preprocessed file for each experiment aes and eae
#	 2.  joined aes and eae file for each experiment
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
###########################################################################

import os
import mgi_utils
import string
import db
import sys
import time             # used for its time.time() function (for timestamps)

# paths to input and two output files
inFilePath= os.getenv('INPUT_FILE')
fpInfile = open(inFilePath, 'r')
logDir =  os.getenv('LOGDIR')
inputDir =  os.getenv('INPUTDIR')
outputDir = os.getenv('OUTPUTDIR')
binDir = '%s/bin' % os.getenv('INSTALLDIR')

# curation and diagnostic logs
fpCur = open (os.environ['LOG_CUR'], 'a')
fpDiag = open (os.environ['LOG_DIAG'], 'a')

# bcp file
bcpFilePath = os.getenv('RNASEQ_BCP')
fpRnaSeqBcp = open(bcpFilePath, 'w')

# constants
TAB = '\t'
CRT = '\n'

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
experimentInDbSet = set()

# GXD HT Sample name in the database
# {sampleID:sampleKey, ...}
sampleInDbDict = {}

# samples with strain J:DO; genotypeKey = 90560 in the db
JDOSampleSet = set()

# samples flagged as relevant in the db
relevantSampleSet = set()

# ensembl IDs assoc w/markers
ensemblMarkerSet = set()

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

def init():
    global experimentInDbSet, sampleInDbList, JDOSampleSet
    global relevantSampleSet, ensemblMarkerSet, ensemblSequenceSet
    global multiEnsemblMarkerDict, multiMarkerEnsemblDict, symbolToMultiEnsIdDict

    db.useOneConnection(1)

    results = db.sql('''select accid 
        from ACC_Accession
	where _MGIType_key = 42
	and preferred = 1''', 'auto')
    for r in results:
	experimentInDbSet.add(r['accid'])

    results = db.sql('''select name, _Sample_key
	from GXD_HTSample''', 'auto')
    for r in results:
	sampleID = string.strip(r['name'])
	sampleInDbDict[sampleID] = r['_Sample_key']

    results =  db.sql('''select name
        from GXD_HTSample
	where _Genotype_key = 90560''', 'auto')
    for r in results:
        JDOSampleSet.add(string.strip(r['name']))

    results =  db.sql('''select name
        from GXD_HTSample
        where _Relevance_key = 20475450''', 'auto')
    for r in results:
        relevantSampleSet.add(string.strip(r['name']))

    results = db.sql('''select accid
	from ACC_Accession
	where _LogicalDB_key =  60
	and _MGIType_key = 2
	and preferred = 1 ''', 'auto')
    for r in results:
	ensemblMarkerSet.add(r['accid'])

    db.sql('''select accid
 	into temporary table multiEns
        from ACC_Accession
        where _LogicalDB_key =  60
        and _MGIType_key = 2
        and preferred = 1
	group by accid having count(*) > 1 ''', None)

    db.sql('''create index idx1 on multiEns(accid)''')

    results = db.sql('''select a.accid, m.symbol
	from multiEns me, ACC_Accession a, MRK_Marker m
	where me.accid = a.accid
	and a._LogicalDB_key =  60
		and a._MGIType_key = 2
		and a.preferred = 1
	and a._Object_key = m._Marker_key''', 'auto')

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
    
    results = db.sql('''select a.accid, m.symbol
	from ACC_Accession a, multi mm, MRK_Marker m
	where a._MGIType_key = 2
	and a._LogicalDB_key = 60
	and a._Object_key = mm._Object_key
	and a._Object_key = m._Marker_key
	order by m.symbol''', 'auto')

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

    return 0

# end init() -------------------------------------------------

def mean(numList):
    raw =  float(sum(numList)) / max(len(numList), 1)

    return round(raw, 2)

# end mean () -------------------------------------------------

def ppAESFile(expID):
    global noRunOrSampleColList, sampleIdNotInDbList, currentAesPPFile
    global aesRunIdSet

    # run IDs in the AES set
    aesRunIdSet = set()

    # path to the aes input file for this experiment
    aesFile = aesTemplate % expID
    try:
	fpAes = open(aesFile, 'r')
    except:
	return 1 # file does not exist
    # path the aes preprocessed file for this experiment
    currentAesPPFile = aesPPTemplate %  expID
    fpCurrentPP = open(currentAesPPFile, 'w')
    ppList = [] # so we can wee out dupes before writing to the PP file
    #
    # process the header line
    #
    header = fpAes.readline()
    headerList = string.split(header, TAB)

    # find the idx of the columns we want - they are not ordered
    sourceSampleIDX = None
    enaSampleIDX = None
    runIDX = None
    for idx, colName in enumerate(headerList):
        colName = string.strip(colName)
	if string.find(colName, 'Source Name') != -1:
	    sourceSampleIDX = idx
	elif string.find(colName, 'ENA_SAMPLE') != -1:
	    enaSampleIDX = idx
	elif string.find(colName, 'ENA_RUN') != -1:
	    runIDX = idx
    #
    # now iterate through each gene/run/tpm in the file
    #
    for line in fpAes.readlines():
	if line == '\n':
	    continue 
        tokens = string.split(line, TAB)
	sampleID = '' # The sampleID we load to the lookup
	ssID = '' # ID from 'Source Name' column
	enaID = '' # ID from ENA_SAMPLE column 
	runID = '' # ID from ENA_RUN column
	if runIDX == None:
	    # report that the READ column cannot be found
	    noRunOrSampleColList.append('Exp: %s file has no Run column' % expID)
	    return 2

	# get the sampleID from one or both of the fields
  	if sourceSampleIDX != None:
	    ssID = string.strip(tokens[sourceSampleIDX])

	if enaSampleIDX != None:
	    enaID = string.strip(tokens[enaSampleIDX])

	# report and return if we did not find either sample column
	if not (ssID or enaID): 
	    noRunOrSampleColList.append('Exp: %s file has neither Sample column' % expID)
	    return 3
	# choose the ID to use: not empty and in database. Checking the database
	# is the only way we know if we have a sampleID (and not something else)
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
	else:
	    runID = string.strip(tokens[runIDX])
	    aesRunIdSet.add(runID)
	    line = '%s%s%s%s' % (runID, TAB, sampleID, CRT)
	    if line not in ppList:
		ppList.append(line)
	# we write all samples except those not in the database to the file. 
	# We will exclude non-relevant and J:DO strain samples later so we may 
	# differentiate btwn a)non-relevant b)J:DO c)not in the aes file 
	# (excluding those not in the database
    for line in ppList:
	fpCurrentPP.write(line)
    fpCurrentPP.close();

    return 0

# end ppAESFile ()---------------------------------------------------

def splitOffGene(tokens):
    del tokens[0] # remove 'Gene ID' 

    return tokens

# end splitOffGene ()--------------------------------------------

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


    return 0

# end writeQC ()--------------------------------------------

def ppEAEFile(expID):

    global noTpmValueList, multiMarkerEnsemblList, currentEaePPFile
    global  ensemblNotInMGIList, ensemblIsOrphanList, eaeRunIdSet
    global multiEnsemblMarkerList
    # we create a set here, because eaeRunIdList will have NL char
    # after last run ID, need this later for set difference
    eaeRunIdSet = set()

    start_time = time.time()

    # create file name and open file descriptor
    #print 'expID: %s' % expID
    eaeFile = eaeTemplate % expID
    try:
	fpEae = open(eaeFile, 'r')
    except:
	return 1 # file does not exist
    # path the aes preprocessed file for this experiment
    currentEaePPFile = eaePPTemplate % expID
    fpCurrentPP = open(currentEaePPFile, 'w')

    # get the runIDs from the header, we will use the index of te`e tpms
    # in each gene/tpm line to get the runID
    header = string.split(fpEae.readline(), TAB)
    eaeRunIdList = splitOffGene(header)
    print 'Num runIDs for %s: %s %s' % (expID, len(eaeRunIdList), eaeRunIdList)
    sys.stdout.flush()

    # process each gene and it's sample tpm's
    for line in fpEae.readlines():
	#print 'line: %s' % line
	tokens = string.split(line, TAB)
	geneID = string.strip(tokens[0])
	#print geneID
	# multi marker per ensembl
	if geneID in multiMarkerEnsemblDict:
	    symbols = string.join(multiMarkerEnsemblDict[geneID], ', ')
	    msg = '%s: %s' % (geneID, symbols)
	    if msg not in multiMarkerEnsemblList:
		multiMarkerEnsemblList.append(msg)
	    continue
	#multiEnsemblMarkerDict
	#symbolToMultiEnsIdDict
	# multi ensembl per marker
	elif geneID in multiEnsemblMarkerDict:
	    #if geneID in ('ENSMUSG00000113702', 'ENSMUSG00000113662', 'ENSMUSG00000113781'): # ('Gm35558', 'Gm2464'):
		#print 'DEBUG multiEnsemblMarker: geneID: %s symbol: %s' % (geneID, multiEnsemblMarkerDict[geneID])
	    symbol = multiEnsemblMarkerDict[geneID]
	    #if len(symbolToMultiEnsIdDict[symbol]) > 1:
	    ensIDs = string.join(symbolToMultiEnsIdDict[symbol], ', ')
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

	# list of runs for geneID
	tpmList = splitOffGene(tokens)

	for idx, tpm in enumerate(tpmList):
	    #print 'idx: %s tpm: %s' % (idx, tpm)
	    # get the runID that goes with this tpm value
            runID = string.strip(eaeRunIdList[idx])
	    eaeRunIdSet.add(runID)
	    tpm = string.strip(tpm)
	    # report blank tpm values and set them to 0.0
	    if tpm == '':
		tpm = 0.0
		noTpmValueList.append('ExpID: %s geneID: %s runID: %s' % \
		    (expID, geneID, runID))

	    fpCurrentPP.write('%s%s%s%s%s%s' % (geneID, TAB, runID, TAB, tpm, CRT))

    fpCurrentPP.close();
    elapsed_time = time.time() - start_time
    print '%sTIME: Processing Runs from the EAE file %s %s%s' % (CRT, expID, time.strftime("%H:%M:%S", time.gmtime(elapsed_time)), CRT) 
    return 0

# end ppEAEFile ()--------------------------------------------

def process():
    global  experimentNotInDbList, runIdNotInAEList
    global nonRelSkippedList, JDOSkippedList
    
    # for each expID in Connie's file
    for line in fpInfile.readlines():
        expID = string.strip(string.split(line)[0])
	
        # report if expID not in the database and skip
        if expID not in experimentInDbSet:
            experimentNotInDbList.append(expID)
            continue

        # preprocess the aes file for this expID (run, sample)
        sys.stdout.flush()
        rc = ppAESFile(expID)
        if rc != 0:
            print '''preprocessing AES file failed with rc%s, 
			skipping file for %s''' % (rc, expID)
            continue

        # preprocess the eae file for this expID (gene, run, tpm)
        sys.stdout.flush()
        rc = ppEAEFile(expID)
        if rc != 0:
            print '''preprocessing EAE file failed with rc %s, 
		    skipping file for %s''' % (rc, expID)
            continue

	# check for AES runIDs not in the EAE file
	#print 'aesRunIdSet: %s' % aesRunIdSet
	#print 'eaeRunIdSet: %s' % eaeRunIdSet
	diff = aesRunIdSet.difference(eaeRunIdSet)
	if len(diff):
	    diffString = string.join(diff, ', ')
	    runIdNotInEAList.append('%s: %s' % (expID, diffString))
        # start of the join
	start_time = time.time()

        joinedFile =  joinedPPTemplate % expID
	cmd = "%s/run_join %s %s %s" % (binDir, currentEaePPFile, currentAesPPFile, joinedFile)
	rc = os.system(cmd)
        if rc != 0:
            msg = 'join cmd failed: %s%s' % (cmd, CRT)
            fpDiag.write(msg)

	# end of the join
 	elapsed_time = time.time() - start_time
	print '%sTIME: Creating the join file %s %s%s' % (CRT, expID, time.strftime("%H:%M:%S", time.gmtime(elapsed_time)), CRT)
	sys.stdout.flush()

        # now open the joined file and process
	start_time = time.time()

	geneDict = {}
	fpJoined = open(joinedFile, 'r') 
	line = fpJoined.readline()
	while line:
	    tokens = string.split(line, TAB)
	    geneID = tokens[0]
	    runID = tokens[1]
	    tpm = tokens[2]
	    sampleID = ''
	    # some join files have no column 4, do a try/except
	    try:
		sampleID = string.strip(tokens[3])
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
            if sampleID not in relevantSampleSet:
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
	    geneDict[geneID][sampleID].append(tpm)
	    line = fpJoined.readline()

	# for Maria
	#for geneID in geneDict:
	#    sampleDict = geneDict[geneID]
	#    for s in sampleDict:
	#	fpRnaSeqBcp.write('%s%s%s%s%s%s%s%s' % (expID, TAB, geneID, TAB, s, TAB, string.join(sampleDict[s], ', '),  CRT))
	#calcTpm - grsd-73
	#writeRnaSeq(geneDict)  - grsd-37

	elapsed_time = time.time() - start_time
        print '%sTIME: Iterating through join file %s %s%s' % (CRT, expID, time.strftime("%H:%M:%S", time.gmtime(elapsed_time)), CRT)
    return 0

# end process ()--------------------------------------------

def closefiles():
    fpCur.close()
    fpDiag.close()
    fpInfile.close()
    fpRnaSeqBcp.close()
    db.useOneConnection(0)

    return 0

# end closeFiles ()--------------------------------------------

#
# Main
#

# -------------------------------------------------------------
START_TIME = time.time()

print 'Start time: %s' %  mgi_utils.date()
sys.stdout.flush()
init()
elapsed_time = time.time() - START_TIME
print 'TIME to run init function %s' %  time.strftime("%H:%M:%S", time.gmtime(elapsed_time))
sys.stdout.flush()

# -------------------------------------------------------------
TIME = time.time()

process()
elapsed_time = time.time() - TIME
print 'TIME to run process function %s' % time.strftime("%H:%M:%S", time.gmtime(elapsed_time))
sys.stdout.flush()

# -------------------------------------------------------------
TIME = time.time()

writeQC()
elapsed_time = time.time() - TIME
print 'TIME to run writeQC function %s' % time.strftime("%H:%M:%S", time.gmtime(elapsed_time))
sys.stdout.flush()

# -------------------------------------------------------------

closefiles()

elapsed_time = time.time() - START_TIME

print 'End time: %s'  % mgi_utils.date()

print 'Total run time: %s' %  time.strftime("%H:%M:%S", time.gmtime(elapsed_time))

