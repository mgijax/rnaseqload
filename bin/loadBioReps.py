##########################################################################
#
# Purpose:
# Load the GXD_HTSample_RNASeqSet and GXD_HTSample_RNASeqSetMember
# with the biological replicates       
#
# Usage: 
# Env Vars:
#      1. INPUT_FILE_DEFAULT - Connie's file of experiment IDs
#      2. LOGDIR
#      3. INPUTDIR
#      4. OUTPUTDIR
#      4. LOG_CUR
#      5. LOG_DIAG
#      6. PG_DBUTILS
#
# Inputs:
#	1. INPUTFILE - Connie's file of experiment IDs
#	2. Configuration (see rnaseqload.config)
#
# Outputs:
#	 1. Records in the database
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

import os	 # for system to execute bcp, getenv
import mgi_utils # for log start/end timestamp
import loadlib	 # for bcp friendly date/timestamp
import string
import db
import sys 	 # to flush stdout
import time	 # used for its time.time() function (for timestamps)

# paths to input and two output files
logDir =  os.getenv('LOGDIR')
inputDir =  os.getenv('INPUTDIR')
outputDir = os.getenv('OUTPUTDIR')

# curation and diagnostic logs
fpCur = open (os.getenv('LOG_CUR'), 'a')
fpDiag = open (os.getenv('LOG_DIAG'), 'a')

# bcp stuff
bcpCommand = os.getenv('PG_DBUTILS') + '/bin/bcpin.csh'
bcpCmdList = []

setTable = 'GXD_HTSample_RNASeqSet'
memberTable = 'GXD_HTSample_RNASeqSetMember'

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

    fpSet = open('%s/%s' % (outputDir, setBcp), 'w')
    fpMember =  open('%s/%s' % (outputDir, memberBcp), 'w')

    db.useOneConnection(1)

    # Don't need to do this as drop and reload, keep for if we decide to do updates
    #results = db.sql('''select nextval('gxd_htsample_rnaseqset_seq') as maxKey ''', 'auto')
    #setKey = results[0]['maxKey']
    #print 'setKey: %s' % setKey

    #results = db.sql('''select nextval('gxd_htsample_rnaseqsetmember_seq') as maxKey ''', 'auto')
    #memberKey = results[0]['maxKey']
    #print 'memberKey: %s' % memberKey

    setKey = 1
    memberKey = 1

    return 0

# end init() -------------------------------------------------


def process():
    global setKey, memberKey

    db.sql('''select sm._object_key as _experiment_key, a.accid as exp_id 
        into temporary table temp1
        from ACC_Accession a,  MGI_Set s, MGI_SetMember sm
        where s.name = 'RNASeq Load Experiments'
    and s._Set_key = sm._Set_key
    and sm._Object_key = a._Object_key
    and a._MGIType_key = 42 --GXD_HTExperiment
    and a._LogicalDB_key = 189
    and a.preferred = 1''', None)
    db.sql('''create index idx1 on temp1 (_experiment_key)''', None)
    db.sql('''select t1.*, hts._sample_key, 
            hts.name, hts.age, hts._organism_key, hts._sex_key, hts._stage_key, 
            hts._emapa_key, hts._genotype_key 
        into temporary table temp2 
        from temp1 t1, gxd_htsample hts
        where hts._genotype_key != 90560 --J:DO
        and hts._relevance_key = 20475450 --Yes
        and hts._experiment_key = t1._experiment_key''', None)
    db.sql('''create index idx2 ON temp2 (_emapa_key)''', None)
    db.sql('''select n._object_key as _sample_key, n.note 
        into temporary table temp4 
        from mgi_note n
        where n._notetype_key = 1048 
        and n._mgitype_key = 43 
        ''', None)
    db.sql('''create index idx4 on temp4 (_sample_key)''', None)
    results = db.sql('''select distinct  t2._experiment_key, t2._sample_key, t2.age,
            t2._organism_key, t2._sex_key, t2._stage_key, t2._emapa_key, t2._genotype_key, 
            t4.note 
        from temp2 t2 
        left outer join temp4 t4 on (t2._sample_key = t4._sample_key)''', 'auto')
    #print 'len results: %s' % len(results)
    repliconDict = {}

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
    # creat bcp line for the Set and the SetMembers
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

# end process() -------------------------------------------------

def execBCP ():

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

    db.useOneConnection(0)

    return 0

# end execBCP --------------------------------------------------------------------


#
# Main
#

START_TIME = time.time()

print('Start time: %s' %  mgi_utils.date())
sys.stdout.flush()

# ----------------------------------------------------------------------------------
TIME = time.time()
init()
elapsed_time = time.time() - TIME
print('TIME to run init function %s' % time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
sys.stdout.flush()

# ----------------------------------------------------------------------------------
TIME = time.time()
process()
elapsed_time = time.time() - TIME
print('TIME to run process function %s' % time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
sys.stdout.flush()

# ----------------------------------------------------------------------------------
TIME = time.time()
execBCP()
elapsed_time = time.time() - TIME
print('TIME to run execBCP function %s' % time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
sys.stdout.flush()

total_time = time.time() - START_TIME

print('Total run time: %s' %  time.strftime("%H:%M:%S", time.gmtime(total_time)))

print('End time: %s'  % mgi_utils.date())
