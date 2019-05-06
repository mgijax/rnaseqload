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
#	 4. AES_URL_TEMPLATE - url template for Array Express
#	 5. EAE_URL_TEMPLATE - url template for Expression Atlas
#
# Inputs:
#	1. INPUTFILE - Connie's file of experiment IDs
#	2. ArrayExpress files by experiment
#	3. Expression Atlas files by experiment
#	4. Configuration (see rnaseqload.config)
#
# Outputs:
#	 1.  bcp file for gxd_htsample_rnaseq
#	 2.  Expression Atlas file for each experiment
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

# paths to input and two output files
inFilePath= os.getenv('INPUT_FILE')
fpInfile = open(inFilePath, 'r')
logDir =  os.getenv('LOGDIR')
inputDir =  os.getenv('INPUTDIR')

# curation log
fpCur = open (os.environ['LOG_CUR'], 'a')
#fpDiag = open (os.environ['LOG_DIAG'], 'a')
# constants
TAB = '\t'
CRT = '\n'

# ArrayExpress Sample File URL  Template
aesTemplate = '%s' % os.getenv('AES_LOCAL_FILE_TEMPLATE')

# Mapping of runID to sampleID
runIdToSampleIdDict = {}

# Expression Atlas Experiment file URL Templage
eaeTemplate  =  '%s' % os.getenv('EAE_LOCAL_FILE_TEMPLATE')

# GXT HT Experiment IDs in the database
experimentInDbList = []

# GXD HT Sample name in the database
# {sampleID:sampleKey, ...}
sampleInDbDict = {}

# Churchill samples with genotypeKey = 90560 in the db`:w
churchillSampleList = []

# samples flagged as relevant in the db
relevantSampleList = []

#
# QC data structures
#

# expIDs in File not in database
experimentNotInDbList = []

# runID not in ArrayExpress sampleID/runID file
runIdNotInAEList = []

# sampleID from ArrayExpress sampleID/runID file  not in database
sampleIdNotInDbList = []

# run ID, sample ID columns not in the AES file
noRunOrSampleColList = []

# churchill samples skipped
churchillSkippedList = []

# non-relevant samples skipped
nonRelSkippedList = []

def init():
    global experimentInDbList, sampleInDbList, churchillSampleList
    global relevantSampleList

    db.useOneConnection(1)

    results = db.sql('''select accid 
        from ACC_Accession
	where _MGIType_key = 42
	and preferred = 1''', 'auto')
    for r in results:
	experimentInDbList.append(r['accid'])

    results = db.sql('''select name, _Sample_key
	from GXD_HTSample''', 'auto')
    for r in results:
	sampleID = string.strip(r['name'])
	sampleInDbDict[sampleID] = r['_Sample_key']

    results =  db.sql('''select name
        from GXD_HTSample
	where _Genotype_key = 90560''', 'auto')
    for r in results:
        churchillSampleList.append(string.strip(r['name']))

    results =  db.sql('''select name
        from GXD_HTSample
        where _Relevance_key = 20475450''', 'auto')
    for r in results:
        relevantSampleList.append(string.strip(r['name']))

    return 0

def mean(numList):
    raw =  float(sum(numList)) / max(len(numList), 1)

    return round(raw, 2)

def processAESFile(expID):
    global runIdToSampleIdDict
    
    runIdToSampleIdDict = {}

    aesFile = aesTemplate % expID
    fpAes = open(aesFile, 'r')
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
    if sourceSampleIDX == None:
	print 'Source Name not found for %s' % expID
    if enaSampleIDX == None:
	print 'ENA_SAMPLE not found for %s' % expID
    if runIDX == None:
	print 'ENA_RUN not found for %s' % expID

    # now iterate through each gene/run/tpm in the file
    for line in fpAes.readlines():
        tokens = string.split(line, TAB)
	sampleID = '' # The sampleID we load to the lookup
	ssID = '' # ID from 'Source Name' column
	enaID = '' # ID from ENA_SAMPLE column 
	runID = '' # ID from ENA_RUN column
	if runIDX == None:
	    # report that the READ column cannot be found
	    noRunOrSampleColList.append('Exp: %s file has no Run column' % expID)
	    return 1

	# get the sampleID from one or both of the fields
  	if sourceSampleIDX != None:
	    ssID = string.strip(tokens[sourceSampleIDX])

	if enaSampleIDX != None:
	    enaID = string.strip(tokens[enaSampleIDX])

	# report and return if we did not find either sample column
	if not (ssID or enaID): 
	    noRunOrSampleColList.append('Exp: %s file has neither Sample column' % expID)
	    return 2

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
	    runIdToSampleIdDict[runID] = sampleID
	# we load all samples except those not in the database. We will exclude
	# non-relevant and churchill samples later so we may differentiate btwn
	# a)non-relevant b)churchill c)not in the aes file (excluding those not
	# in the database

    return 0

def splitOffGene(tokens):
    del tokens[0] # remove 'Gene ID' 

    return tokens

def writeQC():
    if len(experimentNotInDbList):
	fpCur.write('%sExperimentIDs From Connie File Not in the Database%s' % (CRT, CRT))
	fpCur.write('----------------------------------------------------%s' % CRT)
	for e in experimentNotInDbList:
	    fpCur.write('%s%s' % (e, CRT))
	fpCur.write('Total: %s' % len(experimentNotInDbList))

    if len(runIdNotInAEList):
        fpCur.write('%sRun IDs not in AE File%s' % (CRT, CRT))
        fpCur.write('-------------------------%s' % CRT)
        for e in runIdNotInAEList:
            fpCur.write('%s%s' % (e, CRT))
        fpCur.write('Total: %s' % len(runIdNotInAEList))

    if len(sampleIdNotInDbList):
        fpCur.write('%sSample IDs from AE file not in the Database%s' % (CRT, CRT))
        fpCur.write('---------------------------------------------%s' % CRT)
        for e in sampleIdNotInDbList:
            fpCur.write('%s%s' % (e, CRT))
        fpCur.write('Total: %s' % len(sampleIdNotInDbList))

    if len(noRunOrSampleColList):
        fpCur.write('%sExperiments with no Run or Sample Columns%s' % (CRT, CRT))
        fpCur.write('--------------------------------------------%s' % CRT)
        for e in noRunOrSampleColList:
            fpCur.write('%s%s' % (e, CRT))
        fpCur.write('Total: %s' % len(noRunOrSampleColList))

    if len(nonRelSkippedList):
        fpCur.write('%sSamples not flagged as Relevant in the Database%s' % (CRT, CRT))
        fpCur.write('-------------------------------------------------%s' % CRT)
        for e in nonRelSkippedList:
            fpCur.write('%s%s' % (e, CRT))
        fpCur.write('Total: %s' % len(nonRelSkippedList))

    if len(churchillSkippedList):
        fpCur.write('%sSamples in the Churchill set in the Database%s' % (CRT, CRT))
        fpCur.write('-------------------------------------------------%s' % CRT)
        for e in churchillSkippedList:
            fpCur.write('%s%s' % (e, CRT))
        fpCur.write('Total: %s' % len(churchillSkippedList))

    return 0

def process():
    # for each expID in Connie's file
    for line in fpInfile.readlines():
	expID = string.strip(string.split(line)[0])

	# report if expID not in the database and skip
	if expID not in experimentInDbList:
	    experimentNotInDbList.append(expID)
	    continue

	# load run/sample mapping from the aes file for this expID
	print '%sTIME: Calling processAESFile %s %s%s' % (CRT, expID, mgi_utils.date(), CRT)
	sys.stdout.flush()
	rc = processAESFile(expID)
	if rc != 0:
	    print 'Processing AES file failed with rc%s, skipping file for %s' % (rc, expID)
	    continue

	# process the EAE file
	# create file name and open file descriptor
        print 'expID: %s' % expID
	eaeFile = eaeTemplate % expID
	fpEae = open(eaeFile, 'r')

	# get the runIDs from the header
	header = string.split(fpEae.readline(), TAB)
	runIdList = splitOffGene(header)
	print 'Num runIDs: %s %s' % (len(runIdList), runIdList)
	
	# process each gene and it's sample tpm's
	print '%sTIME: Processing Runs from the EAE file %s %s%s' % (CRT, expID, mgi_utils.date(), CRT)
	sys.stdout.flush()
	sys.stdout.flush()
	for line in fpEae.readlines():
	    tokens = string.split(line, TAB)
	    geneID = tokens[0]
	    ###print 'geneID: %s' % geneID

	    # list of runs for geneID
	    tpmList = splitOffGene(tokens)

	    # The samples for this gene with their one or more TPM values
	    sampleDict = {} # {sampleID: [tmp1, ...], ...}

	    for idx, tpm in enumerate(tpmList):
		tpm = string.strip(tpm)

		# get the runID that goes with this tpm value
		runID = string.strip(runIdList[idx])

		# can the runID be translated to a sampleID? If not report
		if runID not in runIdToSampleIdDict:
		    msg = '%s: %s' % (expID, runID)
		    if msg not in  runIdNotInAEList:
			runIdNotInAEList.append(msg)
		    continue

		# runID can be translated, get the sampleID
		sampleID = string.strip(runIdToSampleIdDict[runID])

		# is sampleID flagged as relevant in the db? If not report/skip
		if sampleID not in relevantSampleList: 
		    nonRelSkippedList.append(sampleID)
		    continue
		elif sampleID in churchillSampleList:
		    # Report/skip if the sampleID is in the churchill set
		    churchillSkippedList.append(sampleID)
		    continue
		# Add to the sample dictionary for this gene
		# later we will average all tpm/sample for this gene
		# added runID to the tpm just for debugging
		if sampleID not in sampleDict:
		    sampleDict[sampleID] = []
		sampleDict[sampleID].append('%s|%s' % (runID,tpm))
	    ###print 'Printing samples for gene ID: %s' % geneID
	    #for sampleID in sampleDict:
	#	if len(sampleDict[sampleID]) > 2:
	#	    print 'sample with runs > 2 %s has %s runs' % (sampleID, len(sampleDict[sampleID]))
		###print 'sampleID:%s Runs: %s' % (sampleID, sampleDict[sampleID])

def closefiles():
    fpCur.close()
    #fpDiag.close()
    fpInfile.close()
    db.useOneConnection(0)

    return 0

#
# Main
#

print 'TIME: Calling init %s' % mgi_utils.date()
sys.stdout.flush()
init()

print 'TIME: Calling process %s' % mgi_utils.date()
sys.stdout.flush()
process()

print 'TIME Calling writeQC %s' % mgi_utils.date()
sys.stdout.flush()
writeQC()

print 'TIME: Calling closefiles %s' % mgi_utils.date()
sys.stdout.flush()
closefiles()

print '%s' % mgi_utils.date()

