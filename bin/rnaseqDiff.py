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
#   Pre-Processed Differential files: DIFFINPUTDIR/xxx.group.txt, xxx.tpms.txt
#
# Outputs: DIFFOUTPUTDIR
#   GXD_HTSample_RNASeqSet
#   GXD_HTSample_RNASeqSetMember
#   GXD_HTSample_RNASeqCombined
# 
###########################################################################

import os
import sys
import loadlib
import db

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

# Constants
TAB = '\t'
CRT = '\n'
PIPE = '|'

# load date
loaddate = loadlib.loaddate

# rnaseqload MGI_User
createdByKey = 1613

setKey = None
memberKey = None
combinedKey = None
seqey = None

ensemblMarkers = {}
markerEnsembls = {}

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
# create BCP files for RNASeqSet, RNASeqSetMember
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

    	sampleResults = db.sql('''
            select distinct s._experiment_key, s.expID, s.name, s._experiment_key, s._sample_key, s.age,
                    s._organism_key, s._sex_key, s._stage_key, s._emapa_key, s._genotype_key, 
                    n.note 
            from samples s
            left outer join sampleNotes n on (s._sample_key = n._sample_key)
	    where s.expID = '%s'
            ''' % (expID), 'auto')

	#
	# for each sample in the experiment
	#
    	for s in sampleResults:

        	expKey = s['_experiment_key']
        	sampleKey = s['_sample_key']
        	sample = s['name']
        	age = s['age']
        	orgKey = s['_organism_key']
        	sexKey = s['_sex_key']
        	emapaKey = s['_emapa_key']
        	stageKey = s['_stage_key']
        	genotypeKey = s['_genotype_key']

        	note = s['note']
        	if note == None:
            		note = ''

        	fpSet.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (\
            		setKey, TAB, expKey, TAB, provider, TAB, TAB, \
            		age, TAB, orgKey, TAB, sexKey, TAB, \
            		emapaKey, TAB, stageKey, TAB, genotypeKey, TAB, note, TAB, \
            		createdByKey, TAB, createdByKey, TAB, loaddate, TAB, loaddate, CRT))

        	fpMember.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (\
                	memberKey, TAB, setKey, TAB, sampleKey, TAB, createdByKey, TAB, createdByKey, TAB, loaddate, TAB, loaddate, CRT))

        	memberKey += 1
        	setKey += 1

    fpSet.close()
    fpMember.close()

    return 0

# end processRNASet()

#
# create BCP files for RNASeqCombined and RNASeq
#
def processCombined():
    global fpCombined, fpSeq, fpErrorEnsembl, fpErrorMarker
    global combinedKey, seqKey

    ensemblError = {}
    markerError = {}

    #
    # for each expID
    #
    results = db.sql(''' select expID from experiments ''', 'auto')
    for r in results:

        expID = r['expID']

        # number of bioreplicates per experiment by groupSet
        replicates = {}
        bioresults = db.sql('''
                select distinct s.expID, rm._rnaseqset_key, rs.groupset, count(rm._rnaseqsetmember_key) as countMember
                from samples s, gxd_htsample_rnaseqsetmember rm, gxd_htsample_rnaseqset rs
                where s._sample_key = rm._sample_key
                and rm._rnaseqset_key = rs._rnaseqset_key
                and s.expID = '%s'
                group by s.expID, rm._rnaseqset_key, rs.groupset
            ''' % (expID), 'auto')
        for b in bioresults:
            key = b['groupset']
            value = b
            replicates[key] = []
            replicates[key].append(value)

        #
        # read the "tpms" file
        #
        try:
            fpTpms = open('%s/%s.tpms.txt' % (inputDir, expID), 'r')
        except:
            print('skipping: experiment does not exist in %s/%s.tpms.txt' % (inputDir, expID))
            continue

        # read the header from fpTpms
        # generate a groupSet
        headerList = str.split(fpTpms.readline(), '\t')
        groupSet = []
        for h in headerList[3:]:
            groupSet.append(str.strip(h))

        for line in fpTpms.readlines():

            tokens = str.split(line[:-1], TAB)
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

            for g in range(len(groupSet)):
                gKey = groupSet[g]
                # if a mismatch fails, then the groupSet will not be in the RNASeqSet
                if gKey not in replicates:
                    continue
                avgQnTpm = float(tokens[g+3])
                levelKey = calcLevel(avgQnTpm)
                countMember = replicates[gKey][0]['countMember']
                rnaSeqSetKey = replicates[gKey][0]['_rnaseqset_key']

                fpCombined.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (\
                            combinedKey, TAB, rnaSeqSetKey, TAB, markerKey, TAB, levelKey, TAB, \
                            countMember, TAB, avgQnTpm, TAB, \
                            createdByKey, TAB, createdByKey, TAB, loaddate, TAB, loaddate, CRT))
    
                combinedKey += 1

        fpTpms.close()

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
#initCombined()
#processCombined()
#execCombinedBCP()

