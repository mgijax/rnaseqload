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

# paths to input and two output files
inFilePath= os.getenv('INPUT_FILE')
fpInfile = open(inFilePath, 'r')
logDir =  os.getenv('LOGDIR')
inputDir =  os.getenv('INPUTDIR')
outputDir = os.getenv('OUTPUTDIR')
binDir = '%s/bin' % os.getenv('INSTALLDIR')

# curation log
fpCur = open (os.environ['LOG_CUR'], 'a')
fpDiag = open (os.environ['LOG_DIAG'], 'a')
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

def ppAESFile(expID):
    global noRunOrSampleColList, sampleIdNotInDbList, currentAesPPFile
       
    # path to the aes input file for this experiment
    aesFile = aesTemplate % expID
    fpAes = open(aesFile, 'r')

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
	    line = '%s%s%s%s' % (runID, TAB, sampleID, CRT)
	    if line not in ppList:
		ppList.append(line)
	# we write all samples except those not in the database to the file. 
	# We will exclude non-relevant and churchill samples later so we may 
	# differentiate btwn a)non-relevant b)churchill c)not in the aes file 
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

def ppEAEFile(expID):

    global noTpmValueList, currentEaePPFile

    # path the aes preprocessed file for this experiment
    currentEaePPFile = eaePPTemplate % expID
    fpCurrentPP = open(currentEaePPFile, 'w')

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

	    fpCurrentPP.write('%s%s%s%s%s%s' % (geneID, TAB, runID, TAB, tpm, CRT))

    fpCurrentPP.close();

    return 0

# end ppEAEFile ()--------------------------------------------

def process():
    global  experimentNotInDbList, runIdNotInAEList
    global nonRelSkippedList, churchillSkippedList

    # for each expID in Connie's file
    for line in fpInfile.readlines():
        expID = string.strip(string.split(line)[0])
	
        # report if expID not in the database and skip
        if expID not in experimentInDbList:
            experimentNotInDbList.append(expID)
            continue

        # preprocess the aes file for this expID (run, sample)
        print '%sTIME: Calling ppAESFile %s %s%s' % \
	    (CRT, expID, mgi_utils.date(), CRT)
        sys.stdout.flush()
        rc = ppAESFile(expID)
        if rc != 0:
            print '''preprocessing AES file failed with rc%s, 
			skipping file for %s''' % (rc, expID)
            continue

        # preprocess the eae file for this expID (gene, run, tpm)
        print '%sTIME: Calling ppEAEFile %s %s%s' % \
	    (CRT, expID, mgi_utils.date(), CRT)
        sys.stdout.flush()
        rc = ppEAEFile(expID)
        if rc != 0:
            print '''preprocessing EAE file failed with rc %s, 
		    skipping file for %s''' % (rc, expID)
            continue

        joinedFile =  joinedPPTemplate % expID
	print '%sTIME: Calling run_join to create joinedFile: %s %s %s' % (CRT, joinedFile,  mgi_utils.date(), CRT)
	sys.stdout.flush()
	cmd = "%s/run_join %s %s %s" % (binDir, currentEaePPFile, currentAesPPFile, joinedFile)
	fpDiag.write('cmd: %s%s' % (cmd, CRT))
	rc = os.system(cmd)
        if rc != 0:
            msg = 'join cmd failed: %s%s' % (cmd, CRT)
            fpDiag.write(msg)
	print '%sTIME: Done calling run_join %s%s' % (CRT, mgi_utils.date(), CRT)
	sys.stdout.flush()

        # now open the joined file and process
	geneDict = {}
	fpJoined = open(joinedFile, 'r') 
	line = fpJoined.readline()
	print line
	while line:
	    tokens = string.split(line)
	    geneID = tokens[0]
	    runID = tokens[1]
	    tpm = tokens[2]
	    sampleID = ''
	    try:
		sampleID = tokens[3]
	    except:
		# if the runID from the eae file is not in the aes file 
		# the sample column will be missing
		#print '%s %s %s %s' % (geneID, runID, sampleID, tpm)
		if sampleID == '':
		    #print '%s %s %s %s' % (geneID, runID, sampleID, tpm)
		    msg = '%s: %s' % (expID, runID)
		    if msg not in  runIdNotInAEList:
			runIdNotInAEList.append(msg)
		    line = fpJoined.readline()
		    continue
            # is sampleID flagged as relevant in the db? If not report/skip
	    # we do this here rather than excluding when create the aes pp file
	    # so we may differentiate between sample id 1) not in database and 
	    # 2) a) not relevant b) in churchill set
            if sampleID not in relevantSampleList:
                nonRelSkippedList.append(sampleID)
		line = fpJoined.readline()
                continue
            elif sampleID in churchillSampleList:
                # Report/skip if the sampleID is in the churchill set
                churchillSkippedList.append(sampleID)
		line = fpJoined.readline()
                continue

	    #if geneID not in geneDict:
	#	geneDict[geneID] = {}
	#    if sampleID not in geneDict[geneID]:
	#	geneDict[geneID][sampleID] = []
	#    geneDict[geneID][sampleID].append(tpm)
	    line = fpJoined.readline()
	#calcTpm - grsd-73
	#writeRnaSeq(geneDict)  - grsd-37

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

