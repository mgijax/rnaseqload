#!/usr/local/bin/python
#
#  rnaseqQC.py
###########################################################################
#
#  Purpose:
#
#	This script will generate a QC report for a feature relationship
#	    input file
#
#  Usage:
#
#      rnaseqQC.py  filename
#
#      where:
#          filename = path to the input file
#
#  Inputs:
#      - input file as parameter - see USAGE
#
#  Outputs:
#
#      - QC report (${QC_RPT})
#
#  Exit Codes:
#
#      0:  Successful completion
#      1:  An exception occurred
#      2:  Fatal QC errors detected and written to report
#
#  Assumes:
#
#  Implementation:
#
#      This script will perform following steps:
#
#      1) Validate the arguments to the script.
#      2) Perform initialization steps.
#      3) Open the input/output files.
#      4) Generate the QC reports.
#      5) Close the input/output files.
#
#  Notes:  None
#
###########################################################################
#
#  Modification History:
#
#  Date        SE   Change Description
#  ----------  ---  -------------------------------------------------------
#
#  07/17/2019  sc  Initial development
#
###########################################################################

import sys
import os
import string
import mgi_utils
import db
import time

#
#  CONSTANTS
#
TAB = '\t'
CRT = '\n'

USAGE = 'Usage: rnaseqQC.py inputFile'

#
#  GLOBALS
#

# Report file names
qcRptFile = os.getenv('QC_RPT')

# 1 if any QC errors in the input file
hasFatalErrors = 0

#  valid expIDs in the database
expIdList = []

# list of bad IDs in the input file
badIdList = []

#
# Purpose: Validate the arguments to the script.
# Returns: Nothing
# Assumes: Nothing
# Effects: sets global variable, exits if incorrect # of args
# Throws: Nothing
#
def checkArgs ():
    global inputFile

    if len(sys.argv) != 2:
        print USAGE
        sys.exit(1)

    inputFile = sys.argv[1]
    #print 'inputFile: %s' % inputFile
    return

# end checkArgs() -------------------------------

# Purpose: create lookups, open files
# Returns: Nothing
# Assumes: Nothing
# Effects: Sets global variables, exits if a file can't be opened,
#  creates files in the file system, creates connection to a database

def init ():
    global expIdList

    # open input/output files
    openFiles()

    db.useOneConnection(1)

    #
    # create lookups
    #

    # experiment ID lookup
    results = db.sql('''select a.accID as expID
	from ACC_Accession a
	where a._MGIType_key = 42 -- experiment
        and a._LogicalDB_key = 189 --ArrayExpress
	and a.preferred = 1''', 'auto')

    for r in results:
        expIdList.append(r['expID'])

    return

# end init() -------------------------------

#
# Purpose: Open input and output files.
# Returns: Nothing
# Assumes: Nothing
# Effects: Sets global variables.
# Throws: Nothing
#
def openFiles ():
    global fpInput, fpQcRpt

    #
    # Open the input file
    #
    try:
        fpInput = open(inputFile, 'r')
    except:
        print 'Cannot open input file: %s' % inputFile
        sys.exit(1)

    #
    # Open QC report file
    #
    try:
        fpQcRpt = open(qcRptFile, 'w')
    except:
        print 'Cannot open report file: %s' % qcRptFile
        sys.exit(1)

    return

# end openFiles() -------------------------------

#
# Purpose: run all QC checks
# Returns: Nothing
# Assumes: Nothing
# Effects: writes reports to the file system
# Throws: Nothing
#
def runQcChecks ():
    global badIdList

    # current line number we are parsing
    lineCt = 0
    
    #
    # Iterate through the input file to do the remaining QC checks
    #
    for line in fpInput.readlines():
	#print 'line: %s' % line
	lineCt += 1
	expID = string.strip(line)
        if expID not in expIdList:
	    badIdList.append('%-12s   %-68s' % (lineCt, \
            expID))

    # Now write any errors to the report
    #
    writeReport()

    return

# end runQcChecks() -------------------------------

#
# Purpose: writes out errors to the qc report
# Returns: Nothing
# Assumes: Nothing
# Effects: writes report to the file system
# Throws: Nothing
#

def writeReport():
    global hasFatalErrors
    #
    # Now write any errors to the report
    #
    if len(badIdList):
	hasFatalErrors = 1
	fpQcRpt.write(CRT + CRT + string.center('Experiment IDs not in database',60)+\
	    CRT)
	fpQcRpt.write('%-12s  %-68s %s' % ('Line#','Line', CRT))
	fpQcRpt.write(12*'-' + '  ' + 68*'-' + CRT)
	fpQcRpt.write(string.join( badIdList, CRT))
    
    return

# end writeReport() -------------------------------

#
# Purpose: Close the files.
# Returns: Nothing
# Assumes: Nothing
# Effects: Modifies global variables
# Throws: Nothing
#
def closeFiles ():
    fpInput.close()
    fpQcRpt.close()

    return

# end closeFiles) -------------------------------

#
# Main
#
print 'checkArgs(): %s' % time.strftime("%H.%M.%S.%m.%d.%y", \
    time.localtime(time.time()))
sys.stdout.flush()
checkArgs()

print 'init(): %s' % time.strftime("%H.%M.%S.%m.%d.%y", \
    time.localtime(time.time()))
sys.stdout.flush()
init()

print 'runQcChecks(): %s' % time.strftime("%H.%M.%S.%m.%d.%y", \
    time.localtime(time.time()))
sys.stdout.flush()
runQcChecks()

print 'closeFiles(): %s' % time.strftime("%H.%M.%S.%m.%d.%y", 
    time.localtime(time.time()))
sys.stdout.flush()
closeFiles()

db.useOneConnection(0)
print 'done: %s' % time.strftime("%H.%M.%S.%m.%d.%y", 
    time.localtime(time.time()))

if hasFatalErrors == 1 : 
    sys.exit(2)
else:
    sys.exit(0)
