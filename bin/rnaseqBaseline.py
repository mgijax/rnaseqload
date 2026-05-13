##########################################################################
#
# Purpose:
#
# processSets():
#
# For each Experiments from Baseline RNSeq MGI_Set
#   compare GXD_HTSample to BASELINEINPUTDIR/xxx.group.txt
#       GXD_HTSample : name (sample), age, organism, sex, stage, emapa, genotype
#       xxx.group.txt : group (g1, g2, etc.)
#
#   if no mismatch, then add to RNASeqSet
#   else if only mismatch is with sex, then set sex = Pooled, add to RNASeqSet
#   else, skip
#
#   load into GXD_HTSample_RNASeqSet, GXD_HTSample_RNASeqSetMember
#
# processCombined():
#
# Inputs:
#	MGI_Set = Baseline RNASeq Load Experiment
#   Pre-Processed Baseline files: BASELINEINPUTDIR
#
# Outputs: BASELINEOUTPUTDIR
#   GXD_HTSample_RNASeqSet
#   GXD_HTSample_RNASeqSetMember
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

inputDir = os.getenv('BASELINEINPUTDIR')
outputDir = os.getenv('BASELINEOUTPUTDIR')

setTable = 'GXD_HTSample_RNASeqSet'
memberTable = 'GXD_HTSample_RNASeqSetMember'
combinedTable = 'GXD_HTSample_RNASeqCombined'

setBcp = '%s.bcp' % setTable
memberBcp = '%s.bcp' % memberTable
combinedBcp = '%s.bcp' % combinedTable

fpSet = None	# created in init()
fpMember = None
fpCombined = None

provider = 'Expression Atlas'

# Constants
TAB = '\t'
CRT = '\n'
PIPE = '|'

# load date
loaddate = loadlib.loaddate

# rnaseqload MGI_User
createdByKey = 1673

setKey = None
memberKey = None
combinedKey = None

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

    #
    # baseline experiments
    #
    db.sql('''
        select a._object_key as _experiment_key, a.accid as expID 
        into temporary table experiments
        from MGI_Set s, MGI_SetMember m , ACC_Accession a
        where s.name = 'Baseline RNASeq Load Experiment'
        and s._set_key = m._set_key
        and s._mgitype_key = a._mgitype_key
        and m._object_key = a._object_key
        and a._logicaldb_key = 189
        and a.preferred = 1
        ''', None)

    db.sql('''create index idx1 on experiments (_experiment_key)''', None)

    #
    # baseline experiments and samples
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

    return 0

# end init()

#
# initialize Set
#
def initSet():
    global fpSet, fpMember
    global setKey, memberKey

    fpSet = open('%s/%s' % (outputDir, setBcp), 'w')
    fpMember =  open('%s/%s' % (outputDir, memberBcp), 'w')

    results = db.sql('''select nextval('gxd_htsample_rnaseqset_seq') as maxKey ''', 'auto')
    setKey = results[0]['maxKey']

    results = db.sql('''select nextval('gxd_htsample_rnaseqsetmember_seq') as maxKey ''', 'auto')
    memberKey = results[0]['maxKey']

    return 0

# end initSet()

#
# initialize Combined
#
def initCombined():
    global fpCombined
    global combinedKey

    fpCombined =  open('%s/%s' % (outputDir, combinedBcp), 'w')

    results = db.sql('''select nextval('gxd_htsample_rnaseqcombined_seq') as maxKey ''', 'auto')
    combinedKey = results[0]['maxKey']

    return 0

# end initCombined()

#
# calculate the level for the avgrage QN TPM
# returns proper level key for aveQnTpm
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
# end calcLevel

#
# create BCP files for RNASeqSet, RNASeqSetMember
#
def processSets():
    global fpSet, fpMember
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
    #results = db.sql(''' select distinct s.expID from samples s where expID in ('E-MTAB-7637') ''', 'auto')
    results = db.sql(''' select distinct s.expID from samples s ''', 'auto')
    for r in results:

        #
        # read the "group.txt" file
        # create groupDict by group (g1, g2, etc.)
        # each group contains the set of samples that belong to that group
        #
        groupDict = {}
        prevSample = ''
        expID = r['expID']

        try:
            fpGroup = open('%s/%s.group.txt' % (inputDir, expID), 'r')
            for line in fpGroup.readlines():
                tokens = str.split(line, TAB)
                groupSet = tokens[0]
                sample = "'" + str.strip(tokens[3]) + "'"
                key = groupSet
                value = sample
                if sample != prevSample:
                    if key not in groupDict:
                        groupDict[key] = []
                    groupDict[key].append(value)
                    prevSample = sample
            fpGroup.close()
        except:
            print('experiment does not exist in %s/%s.group.txt' % (inputDir, expID))

        #
        # for each group
        #   select MGI rows for all samples in the group
        #
        for groupSet in groupDict:

            checkAllDict = {}
            checkNoSexDict = {}
            sampleKeySet = []

            byGroup = ','.join(groupDict[groupSet])

            sampleResults = db.sql('''
                select distinct s.expID, s.name, s._experiment_key, s._sample_key, s.age,
                    s._organism_key, s._sex_key, s._stage_key, s._emapa_key, s._genotype_key, 
                    n.note 
                from samples s
                left outer join sampleNotes n on (s._sample_key = n._sample_key)
                where s.expID = '%s' and rtrim(s.name) in (%s)
                ''' % (expID, byGroup), 'auto')

            #
            # compare samples _organism_key, age, _emapa_key, _stage_key, _sex_key, _genotype_key
            #
            for s in sampleResults:

                expKey = s['_experiment_key']
                age = s['age']
                orgKey = s['_organism_key']
                sexKey = s['_sex_key']
                emapaKey = s['_emapa_key']
                stageKey = s['_stage_key']
                genotypeKey = s['_genotype_key']

                note = s['note']
                if note == None:
                    note = ''

                sampleKey = s['_sample_key']
                sampleKeySet.append(sampleKey)

                key = '%s|%s|%s|%s|%s|%s|%s|%s' % (expKey, age, orgKey, sexKey, emapaKey, stageKey, genotypeKey, note)
                if key not in checkAllDict:
                    checkAllDict[key] = []
                checkAllDict[key].append(sampleKey)

                key = '%s|%s|%s|%s|%s|%s|%s' % (expKey, age, orgKey, emapaKey, stageKey, genotypeKey, note)
                if key not in checkNoSexDict:
                    checkNoSexDict[key] = []
                checkNoSexDict[key].append(sampleKey)

            #print(len(checkAllDict))
            #print(checkAllDict)

            # no mismatch
            if len(checkAllDict) == 1:

                #print('no mismatch')
                fpSet.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (\
                    setKey, TAB, expKey, TAB, provider, TAB, groupSet, TAB, \
                    age, TAB, orgKey, TAB, sexKey, TAB, \
                    emapaKey, TAB, stageKey, TAB, genotypeKey, TAB, note, TAB, \
                    createdByKey, TAB, createdByKey, TAB, loaddate, TAB, loaddate, CRT))

                for sKey in sampleKeySet:
                    fpMember.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (\
                        memberKey, TAB, setKey, TAB, sKey, TAB, createdByKey, TAB, createdByKey, TAB, loaddate, TAB, loaddate, CRT))
                    memberKey += 1

                setKey += 1

            # only mismatch is due to Sex
            elif len(checkAllDict) > 1 and len(checkNoSexDict) == 1:

                print('mismatch sex only')
                sexKey = 315166

                fpSet.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (\
                    setKey, TAB, expKey, TAB, provider, TAB, groupSet, TAB, \
                    age, TAB, orgKey, TAB, sexKey, TAB, \
                    emapaKey, TAB, stageKey, TAB, genotypeKey, TAB, note, TAB, \
                    createdByKey, TAB, createdByKey, TAB, loaddate, TAB, loaddate, CRT))

                for sKey in sampleKeySet:
                    fpMember.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (\
                        memberKey, TAB, setKey, TAB, sKey, TAB, createdByKey, TAB, createdByKey, TAB, loaddate, TAB, loaddate, CRT))
                    memberKey += 1

                setKey += 1

            # other mismatch
            else:
                print('skipped due to mismatch : %s\n' % (expID))
                print(checkAllDict)
                print(checkNoSexDict)

    fpSet.close()
    fpMember.close()

    return 0

# end processSets()

#
# create BCP files for RNASeqCombined
#
def processCombined():
    global fpCombined, combinedKey

    #
    # for each expID
    #
    results = db.sql(''' select expID from experiments where expid = 'E-GEOD-55966' ''', 'auto')
    #results = db.sql(''' select expID from experiments ''', 'auto')
    for r in results:

        expID = r['expID']

        # number of bioreplicates per experiment by groupSet
        replicates = {}
        results = db.sql('''
                select distinct s.expID, rm._rnaseqset_key, rs.groupset, count(rm._rnaseqsetmember_key) as replicatesCount
                from samples s, gxd_htsample_rnaseqsetmember rm, gxd_htsample_rnaseqset rs
                where s._sample_key = rm._sample_key
                and rm._rnaseqset_key = rs._rnaseqset_key
                and s.expID = '%s'
                group by s.expID, rm._rnaseqset_key, rs.groupset
            ''' % (expID), 'auto')
        for r in results:
            key = r['groupset']
            value = r['replicatesCount']
            replicates[key] = []
            replicates[key].append(value)

        #
        # read the "tpms" file
        #
        try:
            fpTpms = open('%s/%s.tpms.txt' % (inputDir, expID), 'r')

            # read the header from fpTpms
            headerList = str.split(fpTpms.readline(), '\t')
            groupSet = []
            for h in headerList[3:]:
                groupSet.append(str.strip(h))

            for line in fpTpms.readlines():
                tokens = str.split(line[:-1], TAB)
                ensemblId = tokens[0]

                markerKey = tokens[1]
                markerSymbol = tokens[2]

                if markerKey == '0':
                    print('skipping: ensemblId/marker invalid %s, %s\n' % (expID, ensemblId))
                    continue

                for g in range(len(groupSet)):
                    gKey = groupSet[g]
                    aveQnTpm = float(tokens[g+3])
                    levelKey = calcLevel(aveQnTpm)
                    replicatesCount = replicates[gKey][0]

                    fpCombined.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (\
                            combinedKey, TAB, markerKey, TAB, levelKey, TAB, \
                            replicatesCount, TAB, aveQnTpm, TAB, \
                            createdByKey, TAB, createdByKey, TAB, loaddate, TAB, loaddate, CRT))
    
                    combinedKey += 1

            fpTpms.close()
        except:
            print('skipping: experiment does not exist in %s/%s.tpms.txt' % (inputDir, expID))

    fpCombined.close()

    return 0

# end processCombined()

#
# load the bcp files into the database for Set
#
def execSetBCP():

    bcpCmd = '%s %s %s %s %s %s "\\t" "\\n" mgd' % \
    (bcpCommand, db.get_sqlServer(), db.get_sqlDatabase(), setTable, outputDir, setBcp)
    print('%s\n' % bcpCmd)
    os.system(bcpCmd)

    bcpCmd = '%s %s %s %s %s %s "\\t" "\\n" mgd' % \
    (bcpCommand, db.get_sqlServer(), db.get_sqlDatabase(), memberTable, outputDir, memberBcp)
    print('%s\n' % bcpCmd)
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
    (bcpCommand, db.get_sqlServer(), db.get_sqlDatabase(), combinedTable, outputDir, combinedBcp)
    print('%s\n' % bcpCmd)
    os.system(bcpCmd)

    db.sql(''' select setval('gxd_htsample_rnaseqcombined_seq', (select max(_rnaseqcombined_key) from GXD_HTSample_RNASeqCombined)); ''', None)
    db.commit()

    return 0

# end execCombinedBCP()

#
# Main
#

init()
initSet()
processSets()
execSetBCP()
initCombined()
processCombined()
execCombinedBCP()

