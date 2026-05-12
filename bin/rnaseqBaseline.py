##########################################################################
#
# Purpose:
# Load the GXD_HTSample_RNASeqSet and GXD_HTSample_RNASeqSetMember
# with the biological replicates       
#
# Inputs:
#	1. MGI_Set = Baseline RNASeq Load Experiment
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

# bcp stuff
bcpCommand = os.getenv('PG_DBUTILS') + '/bin/bcpin.csh'
bcpCmdList = []

inputDir = os.getenv('BASELINEINPUTDIR')
outputDir = os.getenv('OUTPUTDIR')

setTable = 'GXD_HTSample_RNASeqSet_baseline'
memberTable = 'GXD_HTSample_RNASeqSetMember_baseline'

setBcp = '%s.bcp' % setTable
memberBcp = '%s.bcp' % memberTable

fpSet = None	# created in init()
fpMember = None

expIDString = None # expIDs to plug into query

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
        select sm._object_key as _experiment_key, a.accid as expID 
        into temporary table temp1
        from ACC_Accession a,  MGI_Set s, MGI_SetMember sm
        where s.name = 'Baseline RNASeq Load Experiment'
        and s._Set_key = sm._Set_key
        and sm._Object_key = a._Object_key
        and a._MGIType_key = 42 --GXD_HTExperiment
        and a._LogicalDB_key = 189
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

    groupDict = {}
    results = db.sql(''' select distinct t2.expID, t2.name from temp2 t2 ''', 'auto')
    for r in results:
        expID = r['expID']
        name = r['name']

        try:
            fpGroup = open('%s/%s.group.txt' % (inputDir, expID), 'r')
            for line in fpAes.readlines():
                tokens = line[:-1].split('\t')
                key = tokens[0]
                value = tokens
                if key not in groupDict:
                    groupDict[key] = []
                groupDict[key].append(value)
            fpGroup.close()
        except:
            print('experiment does not exist in %s/%s.group.txt' % (inputDir, expID))
    print(groupDict)

    results = db.sql('''
        select distinct t2.expID, t2.name, t2._experiment_key, t2._sample_key, t2.age,
            t2._organism_key, t2._sex_key, t2._stage_key, t2._emapa_key, t2._genotype_key, 
            t4.note 
        from temp2 t2 
        left outer join temp4 t4 on (t2._sample_key = t4._sample_key)
        ''', 'auto')
    #print 'len results: %s' % len(results)
    repliconDict = {}
    return 0

    # iterate through results and map the replicate to its set of sample keys	
    for r in results:

        expKey = r['_experiment_key']
        age = r['age']
        orgKey = r['_organism_key']
        sexKey = r['_sex_key']
        emapaKey = r['_emapa_key']
        stageKey = r['_stage_key']
        genotypeKey = r['_genotype_key']
        note = r['note']
        if note == None:
            note = ''

        sampleKey = r['_sample_key']

        key = '%s|%s|%s|%s|%s|%s|%s|%s' % (expKey, age, orgKey, sexKey, emapaKey, stageKey, genotypeKey, note)

        if key not in repliconDict:
            repliconDict[key] = set()
        repliconDict[key].add(sampleKey)

    # get the set of sample keys for each replicate
    # create bcp line for the Set and the SetMembers
    for rep in repliconDict:

        expKey, age, orgKey, sexKey, emapaKey, stageKey, genotypeKey, note = str.split(rep, PIPE)
        sampleKeySet = repliconDict[rep]

        fpSet.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (setKey, TAB, expKey, TAB, age, TAB, orgKey, TAB, sexKey, TAB, emapaKey, TAB, stageKey, TAB, genotypeKey, TAB, note, TAB, createdByKey, TAB, createdByKey, TAB, loaddate, TAB, loaddate, CRT))

        for sKey in sampleKeySet:
            fpMember.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (memberKey, TAB, setKey, TAB, sKey,  TAB, createdByKey, TAB, createdByKey, TAB, loaddate, TAB, loaddate, CRT))
            memberKey += 1

        setKey += 1

    fpSet.close()
    fpMember.close()

    return 0

# end process()

def execBCP():

    bcpCmd = '%s %s %s %s %s %s "\\t" "\\n" mgd' % \
    (bcpCommand, db.get_sqlServer(), db.get_sqlDatabase(), setTable, outputDir, setBcp)
    fpDiag.write('%s\n' % bcpCmd)
    os.system(bcpCmd)

    bcpCmd = '%s %s %s %s %s %s "\\t" "\\n" mgd' % \
    (bcpCommand, db.get_sqlServer(), db.get_sqlDatabase(), memberTable, outputDir, memberBcp)
    fpDiag.write('%s\n' % bcpCmd)
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
#execBCP()
print('End time: %s'  % mgi_utils.date())

