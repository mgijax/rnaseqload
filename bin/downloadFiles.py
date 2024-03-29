##########################################################################
# 
# Purpose:
#       Download ArrayExpress and Expression Atlas files by experiment
#
# Usage: downloadFiles.py
# Env Vars:
#	 1. LOGDIR 
#	 2. RAW_INPUTDIR - files are downloaded to this directory
#	 3a. AES_URL_TEMPLATE - url template for Array Express
#	 3b. AES_LOCAL_FILE_TEMPLATE - full path to the local file
#	 4a. EAE_URL_TEMPLATE - url template for Expression Atlas
#	 4b. EAE_LOCAL_FILE_TEMPLATE - full path to the local file
#	 5. DOWNLOAD_OK - if exists then error-free download
#
# Inputs:
#	1. Database - the experiments in the 'RNA Seq Load Experiment' Set
#		and the experiments loaded 
#	2. Configuration (see rnaseqload.config)
#
# Outputs:
#	 1.  ArrayExpress file for each experiment
#	 2.  Expression Atlas file for each experiment
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
import string
import subprocess
import sys
import db

print('%s' % mgi_utils.date())

# paths to input and two output files
logDir =  os.getenv('LOGDIR')
rawInputDir =  os.getenv('RAW_INPUTDIR')

# download log
fpLog = None

# constants
CRT = '\n'

# ArrayExpress Sample File URL  Template
aesTemplate = os.getenv('AES_URL_TEMPLATE')
aesLocalFileTemplate = os.getenv('AES_LOCAL_FILE_TEMPLATE')

# Expression Atlas Experiment file URL Templage
eaeTemplate  = os.getenv('EAE_URL_TEMPLATE')
eaeLocalFileTemplate = os.getenv('EAE_LOCAL_FILE_TEMPLATE')

# flag set if the file download is successful
downloadOK = os.getenv('DOWNLOAD_OK')

# number of files unable to be downloaded
errorCt = 0

# list of files that did not download by type
failedAESList = []
failedEAEList = []

def init():
    global fpLog, rnaSeqSetResults

    # curation log
    fpLog = open (os.getenv('LOG_DOWNLOAD'), 'w')


    # remove all aes and eae input files a
    cmd = 'rm %s/*.aes.*' % rawInputDir
    rc = os.system(cmd)
    if rc != 0:
        msg = 'rm cmd did not succeed: %s%s' % (cmd, CRT)
        fpLog.write(msg)

    cmd = 'rm %s/*.eae.*' % rawInputDir
    rc = os.system(cmd)
    if rc != 0:
        msg = 'rm cmd did not succeed: %s%s' % (cmd, CRT)
        fpLog.write(msg)

    # create the result set of ids to load
    rnaSeqSetResults = db.sql('''select a.accid
    from ACC_Accession a, MGI_Set s, MGI_SetMember sm
    where s.name = 'RNASeq Load Experiments'
    and s._Set_key = sm._Set_key
    and sm._Object_key = a._Object_key
    and a._MGIType_key = 42 --GXD_HTExperiment
    and a._LogicalDB_key = 189
    and a.preferred = 1''', 'auto')

    return 0

# end init() -------------------------------------------------------------

def downloadAES(expID):

    aesFile = aesLocalFileTemplate % expID 
    fpAes = open(aesFile, 'w')
    aesURL = aesTemplate % (expID, expID)
    msg = ''
    # --retry-max-time 0, don't time out retries
    #--max-time - max time in seconds to allow the whole operation to take.
    #  Use "-C -" to tell curl to automatically find out  where/how  to
    #      resume  the  transfer. It then uses the given output/input files
    #      to figure that out. 
    cmd = "curl -C - --retry 5 --retry-delay 5 --retry-max-time 0 --max-time 40  %s" % aesURL
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    stdout = result.stdout
    stderr = result.stderr
    statusCode = result.returncode

    if statusCode != 0:
        msg = '%s statusCode: %s stderr: %s%s' % (aesURL, statusCode, stderr, CRT)
        fpLog.write(msg)
        return statusCode
    else:
        fpAes.write(stdout)
        fpAes.close()

    return 0

# end downloadAES -------------------------------------------------------------

def downloadEAE(expID):

    eaeFile = eaeLocalFileTemplate % expID 
    fpEae = open(eaeFile, 'w')
    eaeURL = eaeTemplate % expID
    msg = ''
    cmd = "curl --compressed --max-time 40 --retry 5 --retry-delay 5 --retry-max-time 0 %s" % eaeURL
    #print('cmd: %s' % cmd)
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    stdout = result.stdout
    stderr = result.stderr
    statusCode = result.returncode

    if statusCode != 0:
        msg = '%s statusCode: %s stderr: %s%s' % (eaeURL, statusCode, stderr, CRT)
        fpLog.write(msg)
        return statusCode
    else:
        fpEae.write(stdout)
        fpEae.close()
        checkEAEFile(eaeFile)
        
    return statusCode

# end downloadEAE -------------------------------------------------------------

def checkEAEFile(file):

    # 
    # First find the length of the file
    #
    fpLog.write('FINDING length of file\n')
    cmd = "cat %s | wc -l" % file

    fpLog.write(cmd + '\n')

    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    stdout = result.stdout
    stderr = result.stderr
    statusCode = result.returncode

    lastLineNum = str.strip(stdout)
    msg = '%sstatusCode: %s stderr: %s lastLineNum: %s%s' % (CRT, statusCode, stderr, lastLineNum, CRT)
    fpLog.write(msg + '\n')

    #
    # Now find the last index of 'Gene ID' in the file
    # if there is > 1 this will be our starting line number
    #
    fpLog.write('FINDING line numbers of Gene ID\n')
    cmd = "grep -n 'Gene ID' %s | cut -d: -f1" % file
    fpLog.write(cmd + '\n')

    startLineNum = 0 
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    stdout = result.stdout
    stderr = result.stderr
    statusCode = result.returncode

    if stdout:
        stdout = str.strip(stdout)
        startLineNum = stdout.split('\n')[-1]

    msg = '%sstatusCode: %s stderr: %s startLineNum: %s%s' % (CRT, statusCode, stderr, startLineNum, CRT)
    fpLog.write(msg  + '\n')

    #
    # Check to see if we have multi 'Gene ID" in file, and if so extract
    # the record we want
    #
    if int(startLineNum) > 1:
        intermediateFile = '%s.int' % file

        fpLog.write('FINDING proper part of the file\n')
        cmd = "sed -n '%sq;%s,%sp' %s > %s" % (int(lastLineNum)+1, int(startLineNum)-1, int(lastLineNum)+1, file, intermediateFile)
        fpLog.write(cmd + '\n')
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        stdout = result.stdout
        stderr = result.stderr
        statusCode = result.returncode
        msg = '%sstatusCode: %s stderr: %s %s' % (CRT, statusCode, stderr, CRT)
        fpLog.write(msg + '\n')

        fpLog.write('SAVING copy of the file\n')
        cmd = 'mv %s %s.bak' % (file, file)
        fpLog.write(cmd + '\n')
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        stdout = result.stdout
        stderr = result.stderr
        statusCode = result.returncode
        msg = '%sstatusCode: %s stderr: %s %s' % (CRT, statusCode, stderr, CRT)
        fpLog.write(msg + '\n')

        fpIn = open(intermediateFile, 'r')
        fpOut = open(file, 'w') # overwrite the bad file with the good file

        firstLine = 1
        for line in fpIn.readlines():
            if firstLine:
                index = str.find('Gene', line)
                line = line[str.find(line, 'Gene'):]
                firstLine = 0
            fpOut.write(line)
        fpIn.close()
        fpOut.close()
    else:
        fpLog.write('Did not find dupe records\n')

# end checkEAEFile ------------------------------------------------------------

def downloadFiles():
    global totalCt, successCt, errorCt, failedAESList, failedEAE

    totalCt = 0
    successCt = 0
    errorCt = 0

    for r in rnaSeqSetResults:
        totalCt += 1
        expID = str.strip(r['accid'])

        fpLog.write('%sDownload files for experiment ID: %s%s' % (CRT, expID, CRT))

        a_rc = downloadAES(expID)
        if a_rc == 9:
            fpLog.write('%s skipping AES file for %s with curl return code 9 Server denied you to change to the given directory%s' % (CRT, expID, CRT))
            failedAESList.append(expID)
            continue
        elif a_rc == 78:
            fpLog.write('%s skipping AES file for %s with curl return code 78 RETR response: 550%s' % (CRT, expID, CRT))
            failedAESList.append(expID)
            continue
        while a_rc != 0 and a_rc != 9 and a_rc != 78:
            fpLog.write('\ntry AES again\n')
            a_rc = downloadAES(expID)

        e_rc = downloadEAE(expID)
        if e_rc == 9:
            errorCt += 1
            fpLog.write('%s skipping EAE file for %s with curl return code 9 Server denied you to change to the given directory%s' % (CRT, expID, CRT))
            failedEAEList.append(expID)
            continue
        elif e_rc == 78:
            errorCt += 1
            fpLog.write('%s skipping EAE file for %s with curl return code 78 RETR response: 550%s' % (CRT, expID, CRT))
            failedEAEList.append(expID)
            continue
        elif e_rc != 0:
            errorCt += 1
            fpLog.write('%s skipping EAE file for %s with curl return code %s %s' % (CRT, expID, e_rc, CRT))
            failedEAEList.append(expID)
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
    if failedAESList:
        fpLog.write('AES files not downloaded:%s' % CRT)
        for e in failedAESList:
            fpLog.write('%s%s' % (e, CRT))
        fpLog.write('Total: %s%s' % (len(failedAESList), CRT))
    if failedEAEList:
        fpLog.write('EAE files not downloaded:%s' % CRT)
        for e in failedEAEList:
            fpLog.write('%s%s' % (e, CRT))
        fpLog.write('Total: %s%s' % (len(failedEAEList), CRT))

fpLog.close()
sys.exit(rc)

print('%s' % mgi_utils.date())
