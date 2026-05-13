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

setBcp = '%s.bcp' % setTable
memberBcp = '%s.bcp' % memberTable

fpSet = None	# created in init()
fpMember = None

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

#
# initialize some things
#
def init():
    global setKey, memberKey, fpSet, fpMember

    # existing set deleted from wrapper script

    fpSet = open('%s/%s' % (outputDir, setBcp), 'w')
    fpMember =  open('%s/%s' % (outputDir, memberBcp), 'w')

    results = db.sql('''select nextval('gxd_htsample_rnaseqset_seq') as maxKey ''', 'auto')
    setKey = results[0]['maxKey']

    results = db.sql('''select nextval('gxd_htsample_rnaseqsetmember_seq') as maxKey ''', 'auto')
    memberKey = results[0]['maxKey']

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
# create BCP files for RNASeqSet, RNASeqSetMember
#
def processSets():
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

        print(groupDict)

        #
        # for each group
        #   select MGI rows for all samples in the group
        #
        for groupSet in groupDict:

            checkAllDict = {}
            checkNoSexDict = {}
            sampleKeySet = []

            byGroup = ','.join(groupDict[groupSet])
            print(byGroup)

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

            print(len(checkAllDict))
            print(checkAllDict)

            # no mismatch
            if len(checkAllDict) == 1:

                print('no mismatch')
                fpSet.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (\
                    setKey, TAB, expKey, TAB, provider, TAB, groupSet, TAB, \
                    age, TAB, orgKey, TAB, sexKey, TAB, \
                    emapaKey, TAB, stageKey, TAB, genotypeKey, TAB, note, TAB, createdByKey, TAB, \
                    createdByKey, TAB, loaddate, TAB, loaddate, CRT))

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
                    emapaKey, TAB, stageKey, TAB, genotypeKey, TAB, note, TAB, createdByKey, TAB, \
                    createdByKey, TAB, loaddate, TAB, loaddate, CRT))

                for sKey in sampleKeySet:
                    fpMember.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (\
                        memberKey, TAB, setKey, TAB, sKey, TAB, createdByKey, TAB, createdByKey, TAB, loaddate, TAB, loaddate, CRT))
                    memberKey += 1

                setKey += 1

            # other mismatch
            else:
                print('skipped due to mismatch')
                print(checkAllDict)
                print(checkNoSexDict)

    fpSet.close()
    fpMember.close()

    return 0

# end processSets()

#
#
#
def processCombined():

#
# _rnaseqcombined_key
# _marker_key
# _level_key
# numberofbiologicalreplicates -> from rnaseqsetmember?
# averagequantilenormalizedtpm
#
    #
    # for each expID
    #
    results = db.sql(''' select expID from experiments where expid = 'E-GEOD-55966' ''', 'auto')
    for r in results:

        #
        # read the "tpms" file
        #
        tpmsDict = {}
        expID = r['expID']

        try:
            fpTpms = open('%s/%s.tpms.txt' % (inputDir, expID), 'r')
            for line in fpTpms.readlines():
                tokens = str.split(line, TAB)
                print(tokens)
                markerKey = tokens[1]
                markerSymbol = tokens[2]
                g1 = tokens[3]
                g2 = tokens[4]
            fpTpms.close()
        except:
            print('experiment does not exist in %s/%s.tpms.txt' % (inputDir, expID))

        #print(tpmsDict)

    return 0

# end processCombined()

#
# load the bcp files into the database
#
def execBCP():

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

# end execBCP()

#
# Main
#

init()
processSets()
#processCombined()
execBCP()

