#!/opt/python/bin/python
# python2.7 required for pandas
##########################################################################
#
# Purpose:
# Compare the 'RNA Seq Load Experiment' set to the RNA SEQ experiments
# Loaded in the database. See returns below
#
# Usage: checkSet.py
# Env Vars:
# Inputs:
#       1. INPUTFILE - Connie's file of experiment IDs
#       2. ArrayExpress files by experiment
#       3. Expression Atlas files by experiment
#       4. Configuration (see rnaseqload.config)
#
# Outputs:
#        1. aes and eae files from source
#        2. preprocessed file for each experiment aes and eae
#        3. joined aes and eae file for each experiment
#        4. bcp files - one per experiment
#        5. curator and diagnostic log
#
# Exit Codes:
#
#      0:  No changes between RNA Seq Load Experiment set and RNA Seq data loaded
#      1:  RNA Seq Load Experiment set is different than the RNA Seq data loaded
#      2:  RNA Seq Load Experiment set is empty/does not exist
#
#  Assumes:  Nothing
#
#  Notes:  None
#
###########################################################################

import os               # for system to execute bcp, getenv
import string
import db
import sys              # to flush stdout
import time             # used for its time.time() function (for timestamps)
import Set

rnaSeqSetSet = set()
rnaSeqExpSet = set()

# select the IDs from 'RNA-Seq Load' MGI_Set
results = db.sql('''select a.accid
    from ACC_Accession a, MGI_Set s, MGI_SetMember sm
    where s.name = 'RNASeq Load Experiments'
    and s._Set_key = sm._Set_key
    and sm._Object_key = a._Object_key
    and a._MGIType_key = 42 --GXD_HTExperiment
    and a._LogicalDB_key = 189
    and a.preferred = 1
    order by a.accid''', 'auto')

for r in results:
    rnaSeqSetSet.add(r['accid'])

results = db.sql('''select distinct a.accid
    from GXD_HTSample_RNASeq r, GXD_HTSample s, ACC_Accession a
    where r._Sample_key = s._Sample_key
    and s._Experiment_key = a._Object_key
    and a._MGIType_key = 42 --GXD_HTExperiment
    and a._LogicalDB_key = 189
    and a.preferred = 1
    order by a.accid''', 'auto')
for r in results:
    rnaSeqExpSet.add(r['accid'])

rc = 0
if len(rnaSeqSetSet) < 1:
    rc = 2
#print rnaSeqSetSet.difference(rnaSeqExpSet)
#print  rnaSeqExpSet.difference(rnaSeqSetSet)
if rnaSeqSetSet != rnaSeqExpSet:
    rc = 1

sys.exit(rc)
# select loaded IDs
