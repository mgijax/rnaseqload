##########################################################################
#
# Purpose:
# Load the GXD_HTSample_RNASeqSet and GXD_HTSample_RNASeqSetMember
# with the biological replicates       
#
# Inputs:
#	MGI_Set = Baseline RNASeq Load Experiment
#
# Outputs:
#	MGI_Set.bcp
#   MGI_SetMember.bcp
# 
###########################################################################

import os
import sys
import mgi_utils
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

def init():
    global setKey, memberKey, fpSet, fpMember

    # existing set deleted from wrapper script

    fpSet = open('%s/%s' % (outputDir, setBcp), 'w')
    fpMember =  open('%s/%s' % (outputDir, memberBcp), 'w')

    results = db.sql('''select nextval('gxd_htsample_rnaseqset_seq') as maxKey ''', 'auto')
    setKey = results[0]['maxKey']

    results = db.sql('''select nextval('gxd_htsample_rnaseqsetmember_seq') as maxKey ''', 'auto')
    memberKey = results[0]['maxKey']

    return 0

# end init()

def process():
    global setKey, memberKey

    db.sql('''
        select a._object_key as _experiment_key, a.accid as expID 
        into temporary table temp1
        from MGI_Set s, MGI_SetMember m , ACC_Accession a
        where s.name = 'Baseline RNASeq Load Experiment'
        and s._set_key = m._set_key
        and s._mgitype_key = a._mgitype_key
        and m._object_key = a._object_key
        and a._logicaldb_key = 189
        and a.preferred = 1
        ''', None)

    db.sql('''create index idx1 on temp1 (_experiment_key)''', None)

    db.sql('''
        select t1.*, hts._sample_key, 
            hts.name, hts.age, hts._organism_key, hts._sex_key, hts._stage_key, 
            hts._emapa_key, hts._genotype_key 
        into temporary table temp2 
        from temp1 t1, gxd_htsample hts
        where hts._genotype_key != 90560 --J:DO
        and hts._relevance_key = 20475450 --Yes
        and hts._experiment_key = t1._experiment_key
        ''', None)

    db.sql('''create index idx2 ON temp2 (_emapa_key)''', None)

    db.sql('''
        select n._object_key as _sample_key, n.note 
        into temporary table temp4 
        from mgi_note n
        where n._notetype_key = 1048 
        and n._mgitype_key = 43 
        ''', None)

    db.sql('''create index idx4 on temp4 (_sample_key)''', None)

    #
    # for each expID
    #
    #results = db.sql(''' select distinct t2.expID from temp2 t2 where expID in ('E-MTAB-7637') ''', 'auto')
    results = db.sql(''' select distinct t2.expID from temp2 t2 ''', 'auto')
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
                gr = tokens[0]
                sample = "'" + str.strip(tokens[3]) + "'"
                key = gr
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
        for g in groupDict:

            checkAllDict = {}
            checkNoSexDict = {}
            sampleKeySet = []

            byGroup = ','.join(groupDict[g])
            print(byGroup)

            sampleResults = db.sql('''
                select distinct t2.expID, t2.name, t2._experiment_key, t2._sample_key, t2.age,
                    t2._organism_key, t2._sex_key, t2._stage_key, t2._emapa_key, t2._genotype_key, 
                    t4.note 
                from temp2 t2 left outer join temp4 t4 on (t2._sample_key = t4._sample_key)
                where t2.expID = '%s' and rtrim(t2.name) in (%s)
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
                fpSet.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (\
                    setKey, TAB, expKey, TAB, age, TAB, orgKey, TAB, sexKey, TAB, \
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

                fpSet.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (\
                    setKey, TAB, expKey, TAB, age, TAB, orgKey, TAB, sexKey, TAB, \
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

# end process()

def execBCP():

    bcpCmd = '%s %s %s %s %s %s "\\t" "\\n" mgd' % \
    (bcpCommand, db.get_sqlServer(), db.get_sqlDatabase(), setTable, outputDir, setBcp)
    print('%s\n' % bcpCmd)
    os.system(bcpCmd)

    bcpCmd = '%s %s %s %s %s %s "\\t" "\\n" mgd' % \
    (bcpCommand, db.get_sqlServer(), db.get_sqlDatabase(), memberTable, outputDir, memberBcp)
    print('%s\n' % bcpCmd)
    os.system(bcpCmd)

    # reset the rnaseq primary key sequence
    db.sql(''' select setval('gxd_htsample_rnaseq_seq', (select max(_rnaseq_key) from gxd_htsample_rnaseq)) ''', None)

    db.commit()

    return 0

# end execBCP()

#
# Main
#

print('Start time: %s' %  mgi_utils.date())
init()
process()
execBCP()
print('End time: %s'  % mgi_utils.date())

