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
outputDir = os.getenv('OUTPUTDIR')

# curation log
fpCur = open (os.environ['LOG_CUR'], 'a')
fpDiag = open (os.environ['LOG_DIAG'], 'a')
# constants
TAB = '\t'
CRT = '\n'

# ArrayExpress Sample File URL  Template
aesTemplate = '%s' % os.getenv('AES_LOCAL_FILE_TEMPLATE')

# Mapping of runID to sampleID
#runIdToSampleIdDict = {}

# Expression Atlas Experiment file URL Template
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

# ensembl IDs assoc w/markers
ensemblMarkerList = []

# ensembl IDs assoc w/sequences
ensemblSequenceList = []

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

# empty TPM values set to 0.0
noTpmValueList = []

def init():
    global experimentInDbList, sampleInDbList, churchillSampleList
    global relevantSampleList, ensemblMarkerList, ensemblSequenceList

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

    #results = db.sql('''select accid
#	from ACC_Accession
#	where _LogicalDB_key =  60
#	and _MGIType_key = 2
#	and preferred = 1 ''', 'auto')
#    for r in results:
#	ensemblMarkerList.append(r['accid'])

#    results = db.sql('''select accid
#        from ACC_Accession
#        where _LogicalDB_key =  60
#        and _MGIType_key = 19
#        and preferred = 1 ''', 'auto')
#    for r in results:
#        ensemblSequenceList.append(r['accid'])

    return 0

# end init() -------------------------------------------------

def mean(numList):
    raw =  float(sum(numList)) / max(len(numList), 1)

    return round(raw, 2)

# end mean () -------------------------------------------------

def createAESTable():
    # aesTable is global set by process function
    fpDiag.write('creating table %s%s' % (aesTable, CRT))
    db.sql('''create table radar.%s (
	runID text not null,
	sampleID text not null)''' % aesTable, None)

    return 0

# end createAESTable () ---------------------------------------

def bcpAESTable(bcpFile):
    # aesTable is global set by process function
    fpDiag.write('bcping %s into %s%s' % (bcpFile, aesTable, CRT))
    db.bcp(bcpFile,  aesTable, schema='radar')
    results = db.sql('''select count(*) from %s''' % aesTable)
    fpDiag.write('results: %s%s' % (results, CRT))

    return 0

# end bcpAESTable ()--------------------------------------------

def createAESIndex():
    # aesTable is global set by process function
    fpDiag.write('creating runID index on table %s%s' % (aesTable, CRT))
    db.sql('''create index idxAES on %s(runID)''' % aesTable, None)

    return 0

# end createAESIndex ()--------------------------------------------

def dropAESTable():
    # aesTable is global set by process function
    fpDiag.write('dropping table %s%s' % (aesTable, CRT))
    db.sql('''drop table %s''' % aesTable, None)

    return 0

# end dropAESTable ()--------------------------------------------

def bcpAESFile(expID):
    global noRunOrSampleColList, sampleIdNotInDbList
       
    # path the aes bcp file for this experiment 
    currentBcp = '%s/%s.aes.bcp' % (outputDir, expID)
    fpCurrentBcp = open(currentBcp, 'w')

    # create the table for the bcp file
    createAESTable()

    # path to the aes input file for this experiment
    aesFile = aesTemplate % expID
    fpAes = open(aesFile, 'r')

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
    if sourceSampleIDX == None:
	print 'Source Name not found for %s' % expID
    if enaSampleIDX == None:
	print 'ENA_SAMPLE not found for %s' % expID
    if runIDX == None:
	print 'ENA_RUN not found for %s' % expID

    #
    # now iterate through each gene/run/tpm in the file
    #
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
	    fpCurrentBcp.write('%s%s%s%s' % (runID, TAB, sampleID, CRT))
	# we load all samples except those not in the database. We will exclude
	# non-relevant and churchill samples later so we may differentiate btwn
	# a)non-relevant b)churchill c)not in the aes file (excluding those not
	# in the database
    fpCurrentBcp.close();
    bcpAESTable(currentBcp)
    createAESIndex()

    return 0

# end bcpAESFile ()---------------------------------------------------

def splitOffGene(tokens):
    del tokens[0] # remove 'Gene ID' 

    return tokens

# end splitOffGene ()--------------------------------------------

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

    if len(noTpmValueList):
        fpCur.write('%sRuns with empty TPM value set to 0.0%s' % (CRT, CRT))
        fpCur.write('-------------------------------------------------%s' % CRT)
        for e in noTpmValueList:
            fpCur.write('%s%s' % (e, CRT))
        fpCur.write('Total: %s' % len(noTpmValueList))

    return 0

# end writeQC ()--------------------------------------------

def createEAETable():
    fpDiag.write('creating table %s%s' % (eaeTable, CRT))
    db.sql('''create table radar.%s (
	geneID text not null,
        runID text not null,
        tpm text not null)''' % eaeTable, None)

    return 0

# end createEAETable ()--------------------------------------------

def bcpEAETable(bcpFile):
    fpDiag.write('bcping %s into %s%s' % (bcpFile, eaeTable, CRT))
    db.bcp(bcpFile,  eaeTable, schema='radar')
    results = db.sql('''select count(*) from %s''' % eaeTable)
    fpDiag.write('results: %s%s' % (results, CRT))

    return 0

# end bcpEAETable ()--------------------------------------------

def createEAEIndex():
    fpDiag.write('creating runID index on table %s%s' % (eaeTable, CRT))
    db.sql('''create index idxEAE on %s(runID)''' % eaeTable, None)

    return 0

# end createEAEIndex ()--------------------------------------------

def dropEAETable():
    fpDiag.write('dropping table %s%s' % (eaeTable, CRT))
    db.sql('''drop table %s''' % eaeTable, None)

    return 0

# end dropEAETable ()--------------------------------------------

def bcpEAEFile(expID):

    global noTpmValueList

    # path the aes bcp file for this experiment
    currentBcp = '%s/%s.eae.bcp' % (outputDir, expID)
    fpCurrentBcp = open(currentBcp, 'w')

    # create the table for the bcp file
    createEAETable()

    # create file name and open file descriptor
    print 'expID: %s' % expID
    eaeFile = eaeTemplate % expID
    fpEae = open(eaeFile, 'r')

    # get the runIDs from the header, we will use the index of te`e tpms
    # in each gene/tpm line to get the runID
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
	#if geneID not in ensemblMarkers:
	#    report
	#    continue
	# elif geneID not in ensemblSequences:
        #    report
        #    continue

	# list of runs for geneID
	tpmList = splitOffGene(tokens)

	for idx, tpm in enumerate(tpmList):

	    # get the runID that goes with this tpm value
            runID = string.strip(runIdList[idx])
	    tpm = string.strip(tpm)

	    # report blank tpm values and set them to 0.0
	    if tpm == '':
		tpm = 0.0
		noTpmValueList.append('ExpID: %s geneID: %s runID: %s' % \
		    (expID, geneID, runID))

	    fpCurrentBcp.write('%s%s%s%s%s%s' % (geneID, TAB, runID, TAB, tpm, CRT))

    fpCurrentBcp.close();
    bcpEAETable(currentBcp)
    createEAEIndex()

    return 0

# end bcpEAEFile ()--------------------------------------------

def process():
    global aesTable, eaeTable, experimentNotInDbList, runIdNotInAEList
    global nonRelSkippedList, churchillSkippedList

    # for each expID in Connie's file
    for line in fpInfile.readlines():
        expID = string.strip(string.split(line)[0])
	
        # report if expID not in the database and skip
        if expID not in experimentInDbList:
            experimentNotInDbList.append(expID)
            continue

	# table names can't contain '-'
        aesTable = '%s_aes' % expID
	aesTable = aesTable.replace('-', '_')
        eaeTable = '%s_eae' % expID
        eaeTable = eaeTable.replace('-', '_')

        # load run/sample mapping from the aes file into db for this expID
        print '%sTIME: Calling bcpAESFile %s %s%s' % \
	    (CRT, expID, mgi_utils.date(), CRT)
        sys.stdout.flush()
        rc = bcpAESFile(expID)
        if rc != 0:
            print 'bcp AES file failed with rc%s, skipping file for %s' % \
		(rc, expID)
            continue

        #  load gene/run/tpm mapping from eae file into db for this expID
        print '%sTIME: Calling bcpEAEFile %s %s%s' % \
	    (CRT, expID, mgi_utils.date(), CRT)
        sys.stdout.flush()
        rc = bcpEAEFile(expID)
        if rc != 0:
            print 'bcp EAE file failed with rc %s, skipping file for %s' % \
		(rc, expID)
            continue

	# using outer join in case there are runIDs missing from the aes file
	# query the two tables for the data
	# changed to union/not exists to help performance
	print '%sTIME: Querying to join file data %s %s%s' % \
            (CRT, expID, mgi_utils.date(), CRT)
        sys.stdout.flush()
        results = db.sql('''select distinct eae.geneID, eae.tpm, 
		aes.sampleID, eae.runID
		from %s eae, %s aes
		where eae.runID = aes.runID
	        union
		select distinct eae.geneID, eae.tpm, null,  eae.runID
		from %s eae
		where not exists( select 1
		from %s aes where eae.runID = aes.runID)''' % (eaeTable, aesTable, eaeTable, aesTable), 'auto')
	geneDict = {}
	print '%sTIME: DONE querying to join file data %s %s%s' % \
            (CRT, expID, mgi_utils.date(), CRT)
        sys.stdout.flush()
	for r in results:
	    geneID = r['geneID']
	    runID = r['runID']
	    sampleID = r['sampleID']
	    tpm = r['tpm']
	    # if the runID from the eae file is not in the aes file the sampleID
	    # will be None since we did a left outer join, report this
	    #print '%s %s %s' % (geneID, runID, sampleID, tpm)
	    if sampleID == None:
	  	#print '%s %s %s %s' % (geneID, runID, sampleID, tpm)
                msg = '%s: %s' % (expID, runID)
                if msg not in  runIdNotInAEList:
		    runIdNotInAEList.append(msg)
		continue
            # is sampleID flagged as relevant in the db? If not report/skip
	    # we do this here rather than excluding when create the aes bcp file
	    # so we may differentiate between sample id 1) not in database and 
	    # 2) a) not relevant b) in churchill set
            elif sampleID not in relevantSampleList:
                nonRelSkippedList.append(sampleID)
                continue
            elif sampleID in churchillSampleList:
                # Report/skip if the sampleID is in the churchill set
                churchillSkippedList.append(sampleID)
                continue

	    if geneID not in geneDict:
		geneDict[geneID] = {}
	    if sampleID not in geneDict[geneID]:
		geneDict[geneID][sampleID] = []
	    geneDict[geneID][sampleID].append(tpm)
	#calcTpm - grsd-73
	#writeRnaSeq(geneDict)  - grsd-37


	dropAESTable()
	dropEAETable()

    return 0

# end process ()--------------------------------------------

def closefiles():
    fpCur.close()
    fpDiag.close()
    fpInfile.close()
    db.useOneConnection(0)

    return 0

# end closeFiles ()--------------------------------------------

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

