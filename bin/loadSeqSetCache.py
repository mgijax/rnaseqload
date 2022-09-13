##########################################################################
#
# Purpose:
# Load the GXD_HTSample_RNASeqSet_Cache
#
# A temp table is loaded in the wrapper script because the data set
# is so large that python can't handle it -  process is killed as follows:
#    loadSeqSetCache.sh: line 62: 19085 Killed                  ${PYTHON} ${RNASEQLOAD}/bin/loadSeqSetCache.py
#
# In this script we query the temp table and write to the bcp file
# in two parts. Then we execute bcp
#
# Usage: 
# Env Vars:
#      1. OUTPUTDIR
#      2. PG_DBUTILS
#
# Inputs:
#	1. mgd database
#	2. Configuration (see rnaseqload.config)
#
# Outputs:
#	 1. bcp files
#        2. records in the database
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
import db
import sys 	 # to flush stdout
import time	 # used for its time.time() function (for timestamps)

# paths to input and two output files
outputDir = os.getenv('OUTPUTDIR')

# bcp stuff
bcpCommand = os.getenv('PG_DBUTILS') + '/bin/bcpin.csh'

cacheTable = 'GXD_HTSample_RNASeqSet_Cache'

cacheBcp = '%s.bcp' % cacheTable

fpCache = None	# created in init()

# Constants
TAB = '\t'
CRT = '\n'

# load date
loaddate = loadlib.loaddate

# rnaseqload MGI_User
createdByKey = 1613

assocKey = None

def init():
    global assocKey, fpCache

    fpCache = open('%s/%s' % (outputDir, cacheBcp), 'w')

    db.useOneConnection(1)

    assocKey = 1

    return 0

# end init() -------------------------------------------------


def process():
    global assocKey

    TIME = time.time()

    results = db.sql('''select max(_rnaseqset_key) as ct
            from gxd_htsample_rnaseqset''', 'auto')
    
    ct = int(results[0]['ct'])
    print ('ct: %s' % ct)

    ct1 = ct//2

    # if odd number ct2 will be ct1 + 1
    ct2 = ct - ct1

    print ('count: %s ct1: %s ct2: %s' % (ct, ct1, ct2))

    results = db.sql('''select * 
            from rnaseqsetcombined
            where _rnaseqset_key between 1 and %s''' % ct1, 'auto')

    elapsed_time = time.time() - TIME
    print('TIME to run query %s' % time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    sys.stdout.flush()

    #
    # Write to the bcp file part 1
    #
    for r in results:
        combinedKey = r['_rnaseqcombined_key']
        seqSetKey = r['_rnaseqset_key']
        fpCache.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (assocKey, TAB, combinedKey, TAB, seqSetKey, TAB, createdByKey, TAB, createdByKey, TAB, loaddate, TAB, loaddate, CRT))
        assocKey +=1

    results = db.sql('''select *
            from rnaseqsetcombined
            where _rnaseqset_key between %s and %s''' % (ct1 + 1, ct), 'auto')

    elapsed_time = time.time() - TIME
    print('TIME to run query %s' % time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    sys.stdout.flush()

    #
    # Write to the bcp file part 2
    # 
    for r in results:
        combinedKey = r['_rnaseqcombined_key']
        seqSetKey = r['_rnaseqset_key']
        fpCache.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (assocKey, TAB, combinedKey, TAB, seqSetKey, TAB, createdByKey, TAB, createdByKey, TAB, loaddate, TAB, loaddate, CRT))
        assocKey +=1

    fpCache.close()

    return 0

# end process() -------------------------------------------------

def execBCP ():

    bcpCmd = '%s %s %s %s %s %s "\\t" "\\n" mgd' % \
    (bcpCommand, db.get_sqlServer(), db.get_sqlDatabase(), cacheTable, outputDir, cacheBcp)
    print('%s\n' % bcpCmd)
    os.system(bcpCmd)

    # reset the rnaseq primary key sequence
    db.sql(''' select setval('gxd_htsample_rnaseqset_cache_seq', (select max(_assoc_key) from gxd_htsample_rnaseqset_cache)) ''', None)

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
