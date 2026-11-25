##########################################################################
#
# Purpose: Create/Load files for each Experiment in Differential set
#
# processRNASet():
#
# For each Experiment from Differential RNASeq MGI_Set
#   for each GXD_HTSammple
#       GXD_HTSample : name (sample), age, organism, sex, stage, emapa, genotype
#       load into GXD_HTSample_RNASeqSet, GXD_HTSample_RNASeqSetMember
#
# processCombined():
#
# For each Experiment from Differential RNASeq MGI_Set
#   load into GXD_HTSample_RNASeq
#   load into GXD_HTSample_RNASeqCombined
#
# Inputs:
#	MGI_Set = Differential RNASeq Load Experiments
#   Pre-Processed Differential files: DIFFINPUTDIR/xxx.group.txt
#
# Outputs: DIFFOUTPUTDIR
#   GXD_HTSample_RNASeqSet
#   GXD_HTSample_RNASeqSetMember
#   GXD_HTSample_RNASeq
#   GXD_HTSample_RNASeqCombined
# 
###########################################################################

import os
import sys
import loadlib
import db
import numpy as np
import pandas as pd
import quantileNormalize

db.setTrace(True)

# bcp stuff
bcpCommand = os.getenv('PG_DBUTILS') + '/bin/bcpin.csh'
bcpCmdList = []

logDir = os.getenv('LOGDIR')
inputDir = os.getenv('DIFFINPUTDIR')
outputDir = os.getenv('DIFFOUTPUTDIR')

setTable = 'GXD_HTSample_RNASeqSet'
memberTable = 'GXD_HTSample_RNASeqSetMember'
combinedTable = 'GXD_HTSample_RNASeqCombined'
seqTable = 'GXD_HTSample_RNASeq'

setBcp = '%s.bcp' % setTable
memberBcp = '%s.bcp' % memberTable
combinedBcp = '%s.bcp' % combinedTable
seqBcp = '%s.bcp' % seqTable

fpSet = None
fpMember = None
fpCombined = None
fpSeq = None
fpErrorEnsembl = None
fpErrorMarker = None
fpErrorResolved = None
fpErrorUnResolved = None
fpErrorSamples = None

provider = 'Expression Atlas'

# load date
loaddate = loadlib.loaddate

# rnaseqload MGI_User
createdByKey = 1613

setKey = None
memberKey = None
combinedKey = None
seqKey = None

ensemblMarkers = {}
markerEnsembls = {}

stdDevCutoff = float(os.getenv('STDDEV_CUTOFF'))

#
# Level bins
#
HIGH = 50430889
MED = 50430890
LOW = 50430891
BELOW_CUTOFF = 50430892

#
# initialize all
#
def init():
    global ensemblMarkers, markerEnsembls

    #
    # differential experiments
    #
    db.sql('''
        select a._object_key as _experiment_key, a.accid as expID 
        into temporary table experiments
        from MGI_Set s, MGI_SetMember m , ACC_Accession a
        where s.name = 'RNASeq Load Experiments'
        and s._set_key = m._set_key
        and s._mgitype_key = a._mgitype_key
        and m._object_key = a._object_key
        and a._logicaldb_key = 189
        and a.preferred = 1
        --and a._object_key = 15497
        ''', None)

    db.sql('''create index idx1 on experiments (_experiment_key)''', None)

    #
    # differential experiments and samples
    # exclude J:DO genotypes
    # only include relevelance = Yes
    #
    db.sql('''
        select e.*, hts._sample_key, 
            hts.name, hts.age, hts._organism_key, hts._sex_key, hts._stage_key, 
            hts._emapa_key, hts._genotype_key 
        into temporary table samples 
        from experiments e, gxd_htsample hts
        where hts._genotype_key != 90560
        and hts._relevance_key = 20475450
        and hts._experiment_key = e._experiment_key
        ''', None)

    db.sql('''create index idx2 on samples (expID)''', None)
    db.sql('''create index idx3 on samples (name)''', None)

    #
    # ensembl id associated with > 1 marker
    #
    ensemblMarkers = {}
    results = db.sql('''
        WITH ensembls AS (
        select a.accid, count(a._object_key)
        from acc_accession a
        where a._mgitype_key = 2
        and a._logicaldb_key = 60
        and a.preferred = 1
        group by a.accid having count(a._object_key) > 1
        )
        select a.accid, a._object_key, m.symbol
        from ensembls e, acc_accession a, mrk_marker m
        where e.accid = a.accid
        and a._mgitype_key = 2
        and a._logicaldb_key = 60
        and a.preferred = 1
        and a._object_key = m._marker_key
        order by m.symbol
        ''', 'auto')
    for r in results:
        key = r['accid']
        value = r
        if key not in ensemblMarkers:
            ensemblMarkers[key] = []
        ensemblMarkers[key].append(value)
        
    #
    # markers associated with > 1 ensembl id
    #
    markerEnsembls = {}
    results = db.sql('''
        WITH ensembls AS (
        select a._object_key, count(a.accid)
        from acc_accession a
        where a._mgitype_key = 2
        and a._logicaldb_key = 60
        and a.preferred = 1
        group by a._object_key having count(*) > 1
        )
        select a.accid, a._object_key, m.symbol
        from ensembls e, acc_accession a, mrk_marker m
        where e._object_key = a._object_key
        and a._mgitype_key = 2
        and a._logicaldb_key = 60
        and a.preferred = 1
        and a._object_key = m._marker_key
        order by m.symbol
        ''', 'auto')
    for r in results:
        key = r['_object_key']
        value = r
        if key not in markerEnsembls:
            markerEnsembls[key] = []
        markerEnsembls[key].append(value)
        
    return 0

# end init()

#
# initialize Set
#
def initRNASet():
    global fpSet, fpMember, fpErrorResolved, fpErrorUnresolved, fpErrorSamples
    global setKey, memberKey

    fpSet = open('%s/%s' % (outputDir, setBcp), 'w')
    fpMember = open('%s/%s' % (outputDir, memberBcp), 'w')

    # errors that are not in the diagnostic report
    fpErrorResolved = open('%s/mismatchResolvedDifferential.error' % (logDir), 'w')
    fpErrorResolved.write('sample sex pooled\n')
    fpErrorResolved.write('experiment\tEA group\n')
    fpErrorUnresolved = open('%s/mismatchUnresolvedDifferential.error' % (logDir), 'w')
    fpErrorUnresolved.write('experiment\tEA group|sample\n')
    fpErrorSamples = open('%s/samplesDifferential.error' % (logDir), 'w')
    fpErrorSamples.write('missing samples differential\n')
    fpErrorSamples.write('col 1: experiment\n')
    fpErrorSamples.write('col 2: # of samples included in expression atlas grouping\n')
    fpErrorSamples.write('col 3: # of relevant samples in HT index\n')
    fpErrorSamples.write('col 4: sample mismatch\n\n')

    results = db.sql('''select nextval('gxd_htsample_rnaseqset_seq') as maxKey ''', 'auto')
    setKey = results[0]['maxKey']

    results = db.sql('''select nextval('gxd_htsample_rnaseqsetmember_seq') as maxKey ''', 'auto')
    memberKey = results[0]['maxKey']

    return 0

# end initRNASet()

#
# initialize Combined
#
def initCombined():
    global fpCombined, fpSeq, fpErrorEnsembl, fpErrorMarker
    global combinedKey, seqKey

    fpCombined = open('%s/%s' % (outputDir, combinedBcp), 'w')
    fpSeq = open('%s/%s' % (outputDir, seqBcp), 'w')
    fpErrorEnsembl = open('%s/ensemblDifferential.error' % (logDir), 'w')
    fpErrorEnsembl.write('ensemblId associated with > 1 marker OR ensemblId not in MGI\n\n')
    fpErrorMarker = open('%s/markerDifferential.error' % (logDir), 'w')
    fpErrorMarker.write('markers associated with > 1 ensemblId\n\n')

    results = db.sql('''select nextval('gxd_htsample_rnaseqcombined_seq') as maxKey ''', 'auto')
    combinedKey = results[0]['maxKey']

    results = db.sql('''select nextval('gxd_htsample_rnaseq_seq') as maxKey ''', 'auto')
    seqKey = results[0]['maxKey']

    return 0

# end initCombined()

#
# calculate the level for the avgrage QN TPM
# returns proper level key for avgQnTpm
#
def calcLevel(avgQnTpm):
    level = None
    if avgQnTpm < 0.5:	
        level = BELOW_CUTOFF 
    elif avgQnTpm >= 0.5 and avgQnTpm <= 10:
        level = LOW
    elif avgQnTpm >= 11 and avgQnTpm <= 1000:
        level = MED
    else:  # avgQnTpm > 1000:
        level = HIGH

    return  level

# end calcLevel

#
# calculates the average of numList
# returns the calculated average
# 
def calcAve(numList, precision, aveQN=False):
    tSum = sum(numList)
    tLen =  max(len(numList), 1)
    ave = tSum / tLen
    fAve = float(ave)

    if aveQN == True:
    	if fAve < 1:
    		precision = 2
    	else:
    		precision = 1

    return round(fAve, precision)

# end calcAve

#
# calculates the aveTPM and Std Dev and ave Std Dev for all samples of the gene 
# returns dictionary of average TPM values (by gene) for each sample of the gene
#
# input:
# ensemblTPM[ensemblId][sample][list of tpms]
# ENSMUSG00000111744 {
#	'IleumWT3': [0.0, 0.0, 0.0, 0.0], 
#	'IleumWT1': [0.0, 0.0, 0.0, 0.0], 
#	'IleumIL17Ako4': [0.0, 0.0, 0.76, 0.76], 
#	'IleumMR1ko2': [0.0, 0.6, 0.0, 0.59], 
#	'IleumIL17Ako2': [0.94, 0.0, 0.0, 0.0], 
#	'IleumMR1ko1': [0.58, 0.0, 0.0, 0.0], 
#	'ColonWT4': [0.0, 0.0, 0.0, 0.0], 
#	'IleumIL17Ako1': [0.0, 0.0, 0.0, 0.0], 
#	'ColonIL17Ako1': [0.0, 1.01, 0.0, 0.0]
#	}
#
# output:  the input for the quantile normalize function
# {sampleKey:{markerKey:aveTPM, ...}, sampleKey2:{markerKey:aveTPM, ...}, ...}
#
def calcTPMAveSD(expID, ensemblId, markerKey, sampleDict):

    #print('calcTPMAveSD(expID, sampleDict): %s' % expID)
    #print(sampleDict)

    # {sampleKey:{markerKey:aveTPM, ...}, sampleKey2:{markerKey:aveTPM, ...}, ...}
    # input to QN
    aveTPMDict = {}
    sampleAveSDDict = {}

    for sampleId in sampleDict:
        sampleKey = sampleKeyByRun[sampleId]
        tpmList = sampleDict[sampleId]
        techReplCt = len(tpmList)       # number of runs for the sample (1-n)
        aveTpm = calcAve(tpmList, 2)       # calculate the average tpm
        stdDev = round(np.std(tpmList), 2)  # calc the SD of the tpms

        # calculate the SD average
        stdDevAve = 0.0
        if aveTpm != 0.0:
            stdDevAve = round(stdDev / aveTpm, 2)

        # load the QN input dictionary 
        if sampleKey not in aveTPMDict:
            aveTPMDict[sampleKey] = {}
            #print('calcTPMAveSD sampleId; %s sampleKey: %s' % (sampleId, sampleKey))
        aveTPMDict[sampleKey][markerKey] = aveTpm
        #print(aveTPMDict[sampleKey])

        # collect the stdDevAve of the technical replicates - we will want to
        # take the average of all aveSD across all genes of a sample to
        # determine samples to include/exclude.
        key = '%s|%s|%s' % (expID, sampleId, sampleKey)
        if key not in sampleAveSDDict:
            sampleAveSDDict[key] = []
        sampleAveSDDict[key].append(stdDevAve)

    # Now iterate through the sampleAveSDDict and calculate the average of each aveSD. 
    # Report any above a threshold weeding out dupes.
    #print('\naveTPMDict:')
    #print(aveTPMDict)
    #print('\nsampleAveSDDict:')
    #print(sampleAveSDDict)
    keyList = sorted(sampleAveSDDict.keys())
    #print('\nkeyList')
    #print(keyList)

    for key in keyList:
        eID, sampleId, sampleKey = str.split(key, '|')
        aveAveSdAllGenes = calcAve(sampleAveSDDict[key], 2)

        # if average of the average SD across genes in > stdDevCutoff (config)
        # report it and remove all genes for that sample from the aveTPMDict
        if aveAveSdAllGenes > stdDevCutoff:
            del aveTPMDict[int(sampleKey)]
            print('calcTPMAveSD: deleted from aveTPMDict: ', key, str(aveAveSdAllGenes), str(stdDevCutoff))

    return aveTPMDict

# end calcTPMAveSD -------------------------------------------------------------

#
# create BCP files for RNASeqSet, RNASeqSetMember
#
# Read the GXD HT Samples
#	if the GXD HT Sample does not exist in the pre-processed sdrf.txt file, then skip it
#
def processRNASet():
    global fpSet, fpMember, fpErrorResolved, fpErrorUnresolved, fpErrorSamples
    global setKey, memberKey

    db.sql('''
        select n._object_key as _sample_key, n.note 
        into temporary table sampleNotes 
        from mgi_note n
        where n._notetype_key = 1048 
        and n._mgitype_key = 43 
        ''', None)

    db.sql('''create index idx4 on sampleNotes (_sample_key)''', None)

    #
    # for each expID
    #
    results = db.sql(''' select expID from experiments ''', 'auto')
    for r in results:

    	expID = r['expID']

    	# create whereClause from set of sample names that in the Sdrf file
    	sampleInSdrf = []
    	try:
    	    fpSdrf = open('%s/%s.sdrf.txt' % (inputDir, expID), 'r')
    	except:
            print('experiment does not exist in %s/%s.sdrf.txt' % (inputDir, expID))
            continue
    	for line in fpSdrf.readlines():
    	    tokens = str.split(line[:-1], '\t')
    	    sampleInSdrf.append("'" + tokens[0] + "'")
    	fpSdrf.close()
    	whereClause = "and s.name in (" + ",".join(sampleInSdrf) + ")"
    	#print('\nsampleInSdrf')
    	#print(expID, sampleInSdrf)
        
	#
	# for each sample in the experiment
	# create sets, group by:
        #     age, _organism_key, _sex_key, _stage_key, _emapa_key, _genotype_key, n.note 
	# only include samples that exist in the Sdrf file (see whereClause)
	#

    	sampleSet = {}

    	sampleResults = db.sql('''
	    select array_to_string(array_agg(distinct s._sample_key),'|') as sampleSet,
                s._experiment_key, s.age, s._organism_key, s._sex_key, s._stage_key, s._emapa_key, s._genotype_key, n.note 
            from samples s
            left outer join sampleNotes n on (s._sample_key = n._sample_key)
	    where s.expID = '%s'
	    %s
            group by s._experiment_key, s.age, s._organism_key, s._sex_key, s._stage_key, s._emapa_key, s._genotype_key, n.note
            ''' % (expID, whereClause), 'auto')

    	for s in sampleResults:
        	note = s['note']
        	if note == None:
            		note = ''

        	key = s['sampleSet']
        	value = s
        	sampleSet[key] = []
        	sampleSet[key].append(value)

    	for s in sampleSet:

        	sampleKeys = s.split('|')

        	expKey = sampleSet[s][0]['_experiment_key']
        	age = sampleSet[s][0]['age']
        	orgKey = sampleSet[s][0]['_organism_key']
        	sexKey = sampleSet[s][0]['_sex_key']
        	emapaKey = sampleSet[s][0]['_emapa_key']
        	stageKey = sampleSet[s][0]['_stage_key']
        	genotypeKey = sampleSet[s][0]['_genotype_key']

        	fpSet.write('%s\t%s\t%s\t\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (\
                	setKey, expKey, provider, age, orgKey, sexKey, emapaKey, stageKey, genotypeKey, note, \
                	createdByKey, createdByKey, loaddate, loaddate))

        	for sKey in sampleKeys:
        		fpMember.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (\
                		memberKey, setKey, sKey, createdByKey, createdByKey, loaddate, loaddate))
        		memberKey += 1

        	setKey += 1

    fpSet.close()
    fpMember.close()

    return 0

# end processRNASet()

#
# create BCP files for RNASeqCombined and RNASeq
#
# for each experiment
# 1. get gene length in basepairs mgigff3/work/mbl_composite.lengths.tsv
# 2. convert raw counts to tpm (apples -> oranges)
# 3. run calcTPMAveSD
# 4. run quantile normalize
# 5. add rows to fpCombined
# 6. add rows to fpSeq
#

def processCombined():
    global fpCombined, fpSeq, fpErrorEnsembl, fpErrorMarker
    global combinedKey, seqKey
    global sampleKeyByRun

    ensemblError = {}
    markerError = {}

    #
    # for each expID
    #
    #results = db.sql(''' select expID from experiments ''', 'auto')
    #results = db.sql(''' select expID from experiments where expID = 'E-GEOD-21860' ''', 'auto')
    results = db.sql(''' select expID from experiments where expID = 'E-MTAB-6547' ''', 'auto')
    #results = db.sql(''' select expID from experiments where expID in('E-MTAB-698','E-MTAB-6547') ''', 'auto')
    for r in results:

        expID = r['expID']

	# infoBySample['name'] = expID, name, _sample_key, _rnaseqset_key
        infoBySample = {}
        sampleResults = db.sql('''
            select distinct s.expID, s.name, s._sample_key, g._rnaseqset_key
            from samples s, GXD_HTSample_RNASeqSetMember g
	    where s._sample_key = g._sample_key
	    and s.expID = '%s'
            ''' % (expID), 'auto')
        for s in sampleResults:
            key = s['name']
            value = s
            infoBySample[key] = []
            infoBySample[key].append(value)
	# GSM543645 1 = {E-GEOD-21860,GSM543645 1,8385,3843}
        print('\ninfoBySample: ', expID)
        print(infoBySample)

        sampleByRun = {}
        sampleKeyByRun = {}
        try:
            fpSdrf = open('%s/%s.sdrf.txt' % (inputDir, expID), 'r')
        except:
            print('experiment does not exist in %s/%s.sdrf.txt' % (inputDir, expID))
            continue
        for line in fpSdrf.readlines():
            tokens = str.split(line[:-1], '\t')
            sampleName = tokens[0]
            key = sampleName
            value = tokens[1].split(',')
            sampleKeyByRun[key] = []
            sampleKeyByRun[key] = infoBySample[tokens[0]][0]['_sample_key']
            # for each run (value):
            for v in value:
                key = v
                sampleByRun[key] = []
                sampleByRun[key] = sampleName
        fpSdrf.close()

        print('\nsampleByRun')
        print(sampleByRun)
        print('\nsampleKeyByRun')
        print(sampleKeyByRun)
        
	# countByRun
        try:
            fpRawCounts = open('%s/%s.raw-counts.txt' % (inputDir, expID), 'r')
        except:
            print('experiment does not exist in %s/%s.raw-countsgroup.txt' % (inputDir, expID))
            continue
        
        # read the header from fpRawCounts
        # generate a runTpmSet
        headerList = str.split(fpRawCounts.readline(), '\t')
        runTpmSet = []
        for h in headerList[5:]:
            runTpmSet.append(str.strip(h))
        print('\nrunTpmSet')
        print(runTpmSet)
        # read these 2 headers; but do nothing
	# these rows are used for testing/troubleshooting
	# use infoBySample() instead
        headerList2 = str.split(fpRawCounts.readline(), '\t')
        headerList3 = str.split(fpRawCounts.readline(), '\t')

        # ensemblTPM = {}
        # {ensemblId: {sampleId:[float(tokens[ss+5], ...], ...}, ...}
        ensemblTPM = {}

        # set for fpSeq
        fpSeqValues = {}

        # ENSMUSG00000000001      9199    Gnai3   116     96      184     41      119     141
        for line in fpRawCounts.readlines():

            tokens = str.split(line[:-1], '\t')
            ensemblId = tokens[0]
            markerKey = int(tokens[1])
            markerSymbol = tokens[2]
                
            # ensemblId has > 1 marker or ensemblId does not exist in MGI
            if ensemblId in ensemblMarkers or markerKey == 0:
                if ensemblId not in ensemblError:
                    ensemblError[ensemblId] = []
                    if ensemblId not in ensemblMarkers:
                        continue
                    for e in ensemblMarkers[ensemblId]:
                        ensemblError[ensemblId].append(e['symbol'])
                continue

            if markerKey in markerEnsembls:
                if markerSymbol not in markerError:
                    markerError[markerSymbol] = []
                    for m in markerEnsembls[markerKey]:
                        markerError[markerSymbol].append(m['accid'])
                continue

            if ensemblId not in ensemblTPM:
                ensemblTPM[ensemblId] = {}

            for ss in range(len(runTpmSet)):
                runId = runTpmSet[ss]
                sampleId = sampleByRun[runId]
                tpm = float(tokens[ss+5])
                if sampleId not in ensemblTPM[ensemblId]:
                    ensemblTPM[ensemblId][sampleId] = []
                ensemblTPM[ensemblId][sampleId].append(tpm)

            #print('\nensemblTPM')
            #print(ensemblId, ensemblTPM[ensemblId])
            #print(sampleKeyByRun)

            # process each marker/sample
            aveTPMDict = calcTPMAveSD(expID, ensemblId, markerKey, ensemblTPM[ensemblId])

            # process QN
            qnInput = pd.DataFrame(aveTPMDict)
            qnOutput = quantileNormalize.qn(qnInput, list(aveTPMDict.keys()))
            qnOutputDict = qnOutput.to_dict()

            # for each sample in this row
            for s in ensemblTPM[ensemblId]:
               sampleKey = sampleKeyByRun[s]
               if sampleKey not in aveTPMDict:
                   #print('skipping: sampleKey was deleted in calcTPMAveSD step: ', s)
                   continue

               rnaSeqSetKey = infoBySample[s][0]['_rnaseqset_key']
               aveTpm = aveTPMDict[sampleKey][markerKey]
               qnTpm = round(qnOutputDict[sampleKey][markerKey], 2)

               # per marker, per rnaSeqSetKey (bioreplicate)
               #print(expID, ensemblId, markerKey, rnaSeqSetKey, sampleKey, aveTpm, qnTpm)
               if markerKey not in fpSeqValues:
                 fpSeqValues[markerKey] = {}
               if rnaSeqSetKey not in fpSeqValues[markerKey]:
                 fpSeqValues[markerKey][rnaSeqSetKey] = []
               fpSeqValues[markerKey][rnaSeqSetKey].append([sampleKey, aveTpm, qnTpm])

        print('\ngenerating fpSeq, fpCombined')
        for marker in fpSeqValues:
           for rnaseqset in fpSeqValues[marker]:
               #print(fpSeqValues[marker][rnaseqset])
               numberofbiologicalreplicants = len(fpSeqValues[marker][rnaseqset])
               qnList = []
               for biorep in fpSeqValues[marker][rnaseqset]:
                   qnList.append(biorep[2])
               aveQnTpm = calcAve(qnList, 2, True)
               levelKey = calcLevel(aveQnTpm)
               #print(numberofbiologicalreplicants, aveQnTpm, levelKey)

               for biorep in fpSeqValues[marker][rnaseqset]:
                   sampleKey = biorep[0]
                   aveTpm = biorep[1]
                   qnTpm = biorep[2]
                   fpSeq.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (\
                        seqKey, sampleKey, combinedKey, marker, aveTpm, qnTpm, \
                        createdByKey, createdByKey, loaddate, loaddate))
                   seqKey += 1
        
               fpCombined.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (\
                        combinedKey, rnaseqset, marker, levelKey, \
                        numberofbiologicalreplicants, aveQnTpm, \
                        createdByKey, createdByKey, loaddate, loaddate))
               combinedKey += 1
        
        fpRawCounts.close()

    fpCombined.close()
    fpSeq.close()

    for e in sorted(ensemblError):
        fpErrorEnsembl.write(e + '\t' + '\t'.join(ensemblError[e]) + '\n')
    fpErrorEnsembl.close()

    for m in sorted(markerError):
        fpErrorMarker.write(m + '\t' + '\t'.join(markerError[m]) + '\n')
    fpErrorMarker.close()

    return 0

# end processCombined()

#
# load the bcp files into the database for Set
#
def execSetBCP():

    bcpCmd = '%s %s %s %s %s %s "\\t" "\\n" mgd' % \
    (bcpCommand, db.get_sqlServer(), db.get_sqlDatabase(), setTable, outputDir, setBcp)
    print('%s' % bcpCmd)
    os.system(bcpCmd)

    bcpCmd = '%s %s %s %s %s %s "\\t" "\\n" mgd' % \
    (bcpCommand, db.get_sqlServer(), db.get_sqlDatabase(), memberTable, outputDir, memberBcp)
    print('%s' % bcpCmd)
    os.system(bcpCmd)

    # reset the primary key sequence
    db.sql(''' select setval('gxd_htsample_rnaseqset_seq', (select max(_rnaseqset_key) from GXD_HTSample_RNASeqSet)) ''', None)
    db.commit()
    db.sql(''' select setval('gxd_htsample_rnaseqsetmember_seq', (select max(_rnaseqsetmember_key) from GXD_HTSample_RNASeqSetMember));
''', None)
    db.commit()

    return 0

# end execSetBCP()

#
# load the bcp files into the database for Combined
#
def execCombinedBCP():

    bcpCmd = '%s %s %s %s %s %s "\\t" "\\n" mgd' % \
    (bcpCommand, db.get_sqlServer(), db.get_sqlDatabase(), seqTable, outputDir, seqBcp)
    print('%s' % bcpCmd)
    os.system(bcpCmd)

    bcpCmd = '%s %s %s %s %s %s "\\t" "\\n" mgd' % \
    (bcpCommand, db.get_sqlServer(), db.get_sqlDatabase(), combinedTable, outputDir, combinedBcp)
    print('%s' % bcpCmd)
    os.system(bcpCmd)

    db.sql(''' select setval('gxd_htsample_rnaseq_seq', (select max(_rnaseq_key) from GXD_HTSample_RNASeq)); ''', None)
    db.commit()

    db.sql(''' select setval('gxd_htsample_rnaseqcombined_seq', (select max(_rnaseqcombined_key) from GXD_HTSample_RNASeqCombined)); ''', None)
    db.commit()

    return 0

# end execCombinedBCP()

#
# Main
#

init()
initRNASet()
processRNASet()
execSetBCP()
initCombined()
processCombined()
execCombinedBCP()

