##########################################################################
# 
# Purpose:
#       Download ArrayExpress and Expression Atlas files by experiment
#
# Usage: downloadBaselineFiles.py
# Env Vars:
#	 1. LOGDIR 
#	 2. BASELINERAW_INPUTDIR - files are downloaded to this directory
#	 3a. EAE_EAT_URL_TEMPLATE - url template for Expression Atlas
#	 3b. EAE_EAG_URL_TEMPLATE - url template for Expression Atlas
#	 3c. EAE_AES_URL_TEMPLATE - url template for Expression Atlas
#	 4. DOWNLOAD_OK - if exists then error-free download
#
# Inputs:
#	1. Database - the experiments in the 'Baseline RNA Seq Load Experiment' Set
#		and the experiments loaded 
#	2. Configuration (see rnaseqload.config)
#
# Outputs:
#	 1.  Expression Atlas file for each experiment
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
import os
import mgi_utils
import subprocess
import sys
import db

print('%s' % mgi_utils.date())

# paths to input and two output files
logDir =  os.getenv('LOGDIR')
rawInputDir =  os.getenv('BASELINERAW_INPUTDIR')

# download log
fpLog = None

# constants
CRT = '\n'

# Expression Atlas Experiment file URL Templage
eatTemplate  = os.getenv('EAE_EAT_URL_TEMPLATE')
eatLocalFileTemplate = os.getenv('EAE_EAT_LOCAL_FILE_TEMPLATE')
eagTemplate  = os.getenv('EAE_EAG_URL_TEMPLATE')
eagLocalFileTemplate = os.getenv('EAE_EAG_LOCAL_FILE_TEMPLATE')
aesTemplate  = os.getenv('EAE_AES_URL_TEMPLATE')
aesLocalFileTemplate = os.getenv('EAE_AES_LOCAL_FILE_TEMPLATE')

# number of files unable to be downloaded
errorCt = 0

# list of files that did not download by type
failedList = []

def init():
    global fpLog, rnaSeqSetResults

    # curation log
    fpLog = open (os.getenv('BASELINELOG_DOWNLOAD'), 'w')

    cmd = 'rm %s/*.eae.*' % rawInputDir
    rc = os.system(cmd)
    if rc != 0:
        msg = 'rm cmd did not succeed: %s%s' % (cmd, CRT)
        fpLog.write(msg)

    # create the result set of ids to load
    rnaSeqSetResults = db.sql('''
    select a.accid
    from ACC_Accession a, MGI_Set s, MGI_SetMember sm
    where s.name = 'Baseline RNASeq Load Experiment'
    and s._Set_key = sm._Set_key
    and sm._Object_key = a._Object_key
    and a._MGIType_key = 42 --GXD_HTExperiment
    and a._LogicalDB_key = 189
    and a.preferred = 1
    ''', 'auto')

    return 0

# end init() -------------------------------------------------------------

def downloadEAT(expID):

    eaeURL = eatTemplate % (expID, expID)
    outputFile = rawInputDir + '/' + expID + '-tpms.tsv'
    cmd = ['wget', eaeURL]
    cmd.extend(['-O', outputFile])
    result = subprocess.run(cmd, check=True)
    stdout = result.stdout
    stderr = result.stderr
    statusCode = result.returncode

    if statusCode != 0:
        fpLog.write('%s statusCode: %s stderr: %s%s' % (eaeURL, statusCode, stderr, CRT))

    return statusCode

# end downloadEAT -------------------------------------------------------------

def downloadEAG(expID):

    eaeURL = eagTemplate % (expID, expID)
    outputFile = rawInputDir + '/' + expID + '-configuration.xml'
    cmd = ['wget', eaeURL]
    cmd.extend(['-O', outputFile])
    result = subprocess.run(cmd, check=True)
    stdout = result.stdout
    stderr = result.stderr
    statusCode = result.returncode

    if statusCode != 0:
        fpLog.write('%s statusCode: %s stderr: %s%s' % (eaeURL, statusCode, stderr, CRT))

    return statusCode

# end downloadEAG -------------------------------------------------------------

def downloadAES(expID):

    # E-GEOD-22131 = E, GEOD, 22131 -> 131 last 3 characters
    tokens = expID.split('-')
    print(tokens)
    a = tokens[0] + '-' + tokens[1] + '-'
    b = tokens[2][-3:]
    eaeURL = aesTemplate % (a, b, expID, expID)
    outputFile = rawInputDir + '/' + expID + '.sdrf.txt'
    cmd = ['wget', eaeURL]
    cmd.extend(['-O', outputFile])
    print(cmd)
    result = subprocess.run(cmd, check=True)
    stdout = result.stdout
    stderr = result.stderr
    statusCode = result.returncode

    if statusCode != 0:
        fpLog.write('%s statusCode: %s stderr: %s%s' % (eaeURL, statusCode, stderr, CRT))

    return statusCode

# end downloadAES -------------------------------------------------------------

def downloadFiles():
    global totalCt, successCt, errorCt, failed

    totalCt = 0
    successCt = 0
    errorCt = 0

    for r in rnaSeqSetResults:
        totalCt += 1
        expID = str.strip(r['accid'])

        fpLog.write('%sDownload files for experiment ID: %s%s' % (CRT, expID, CRT))

        e_rc = downloadEAT(expID)
        if e_rc != 0:
            errorCt += 1
            fpLog.write('%s skipping EAE file for %s with wget return code %s %s' % (CRT, expID, e_rc, CRT))
            failedList.append(expID)
            continue

        e_rc = downloadEAG(expID)
        if e_rc != 0:
            errorCt += 1
            fpLog.write('%s skipping EAE file for %s with wget return code %s %s' % (CRT, expID, e_rc, CRT))
            failedList.append(expID)
            continue

        e_rc = downloadAES(expID)
        if e_rc != 0:
            errorCt += 1
            fpLog.write('%s skipping EAE file for %s with wget return code %s %s' % (CRT, expID, e_rc, CRT))
            failedList.append(expID)
            continue

        successCt += 1

    if errorCt > 0:
        return 1

    return 0

# end downloadFiles -----------------------------------------------------------

#
# Main
#
init()

rc = downloadFiles()

fpLog.write('%sTotal experiments in input file: %s%s' % (CRT, totalCt, CRT))
fpLog.write('%sTotal files unable to be downloaded: %s%s' % (CRT, errorCt, CRT))
fpLog.write('%sTotal files successfully downloaded:  %s%s' % (CRT, successCt, CRT))
if rc > 0:
    fpLog.write("%sDownload did not succeed on one or more experiments%s" % (CRT, CRT))
    if failedList:
        fpLog.write('EAE files not downloaded:%s' % CRT)
        for e in failedList:
            fpLog.write('%s%s' % (e, CRT))
        fpLog.write('Total: %s%s' % (len(failedList), CRT))

fpLog.close()
sys.exit(rc)

print('%s' % mgi_utils.date())

