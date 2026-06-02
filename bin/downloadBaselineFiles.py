##########################################################################
# 
# Purpose:
#       Download ArrayExpress and Expression Atlas files by experiment
#
# Usage: downloadBaselineFiles.py
# Env Vars:
#	 1. LOGDIR 
#	 2. BASELINERAW_INPUTDIR - files are downloaded to this directory
#	 3a. EAE_TPMS_URL_TEMPLATE - url template for Expression Atlas
#	 3b. EAE_GROUP_URL_TEMPLATE - url template for Expression Atlas
#	 3c. AES_SDRF_URL_TEMPLATE - url template for Expression Atlas
#	 4. DOWNLOAD_OK - if exists then error-free download
#
# Inputs:
#	1. Database - the experiments in the 'Baseline RNASeq Load Experiment' Set
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
import sys
import subprocess
import mgi_utils
import db

# paths to logs, input
rawInputDir = os.getenv('BASELINERAW_INPUTDIR')

# Expression Atlas Experiment file URL Templage
eatTemplate  = os.getenv('EAE_TPMS_URL_TEMPLATE')
eagTemplate  = os.getenv('EAE_GROUP_URL_TEMPLATE')
aesTemplate  = os.getenv('AES_SDRF_URL_TEMPLATE')

# number of files unable to be downloaded
errorCt = 0

# list of files that did not download by type
failedList = []

def init():
    global rnaSeqSetResults

    cmd = 'rm %s/*.eae.*' % rawInputDir
    rc = os.system(cmd)
    if rc != 0:
        msg = 'rm cmd did not succeed: %s\n' % (cmd)
        print(msg)

    # create the result set of ids to load
    rnaSeqSetResults = db.sql('''
        select a.accid
        from MGI_Set s, MGI_SetMember m , ACC_Accession a
        where s.name = 'Baseline RNASeq Load Experiments'
        and s._set_key = m._set_key
        and s._mgitype_key = a._mgitype_key
        and m._object_key = a._object_key
        and a._logicaldb_key = 189
        and a.preferred = 1
        ''', 'auto')

    return 0

# end init() -------------------------------------------------------------

def downloadTPMS(expID):

    eaeURL = eatTemplate % (expID, expID)
    outputFile = rawInputDir + '/' + expID + '-tpms.tsv'
    cmd = ['wget', eaeURL]
    cmd.extend(['-O', outputFile])
    result = subprocess.run(cmd, check=True)
    stdout = result.stdout
    stderr = result.stderr
    statusCode = result.returncode

    if statusCode != 0:
        print('%s statusCode: %s stderr: %s\n' % (eaeURL, statusCode, stderr))

    return statusCode

# end downloadTPMS -------------------------------------------------------------

def downloadGROUP(expID):

    eaeURL = eagTemplate % (expID, expID)
    outputFile = rawInputDir + '/' + expID + '-configuration.xml'
    cmd = ['wget', eaeURL]
    cmd.extend(['-O', outputFile])
    result = subprocess.run(cmd, check=True)
    stdout = result.stdout
    stderr = result.stderr
    statusCode = result.returncode

    if statusCode != 0:
        print('%s statusCode: %s stderr: %s\n' % (eaeURL, statusCode, stderr))

    return statusCode

# end downloadGROUP -------------------------------------------------------------

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
        print('%s statusCode: %s stderr: %s\n' % (eaeURL, statusCode, stderr))

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

        print('Download files for experiment ID: %s\n' % (expID))

        e_rc = downloadTPMS(expID)
        if e_rc != 0:
            errorCt += 1
            print('skipping EAE file for %s with wget return code %s\n' % (expID, e_rc))
            failedList.append(expID)
            continue

        e_rc = downloadGROUP(expID)
        if e_rc != 0:
            errorCt += 1
            print('skipping EAE file for %s with wget return code %s\n' % (expID, e_rc))
            failedList.append(expID)
            continue

        e_rc = downloadAES(expID)
        if e_rc != 0:
            errorCt += 1
            print('skipping EAE file for %s with wget return code %s' % (expID, e_rc))
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

print(mgi_utils.date())
init()

rc = downloadFiles()

print('Total experiments in input file: %s\n' % (totalCt))
print('Total files unable to be downloaded: %s\n' % (errorCt))
print('Total files successfully downloaded:  %s\n' % (successCt))
if rc > 0:
    print("Download did not succeed on one or more experiments\n")
    if failedList:
        print('EAE files not downloaded:\n')
        for e in failedList:
            print('%s\n' % (e))
        print('Total: %s\n' % (len(failedList)))

print(mgi_utils.date())

