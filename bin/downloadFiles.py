#!/usr/bin/python
# /usr/bin/python is 2.6.* or 2.7.* depending on server. It is required for urllib2
# /opt/python2.7/bin/python
##########################################################################
# 
# Purpose:
#       Download ArrayExpress and Expression Atlas files by experiment
#
# Usage: downloadFiles.py
# Env Vars:
#	 1. INPUT_FILE - Connie's file of experiment IDs
#	 2. LOGDIR 
#	 3. INPUTDIR - files are downloaded to this directory
#	 4. AES_URL_TEMPLATE - url template for Array Express
#	 5. EAE_URL_TEMPLATE - url template for Expression Atlas
#
# Inputs:
#	1. INPUTFILE - Connie's file of experiment IDs
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
import urllib2
import runCommand

print '%s' % mgi_utils.date()

# paths to input and two output files
inFilePath= os.getenv('INPUT_FILE')
fpInfile = open(inFilePath, 'r')
logDir =  os.getenv('LOGDIR')
inputDir =  os.getenv('INPUTDIR')

# curation log
fpCur = open (os.environ['LOG_CUR'], 'a')
fpDiag = open (os.environ['LOG_DIAG'], 'a')
# constants
CRT = '\n'

# ArrayExpress Sample File URL  Template
aesTemplate = os.getenv('AES_URL_TEMPLATE')
aesLocalFileTemplate = os.getenv('AES_LOCAL_FILE_TEMPLATE')

# Expression Atlas Experiment file URL Templage
eaeTemplate  = os.getenv('EAE_URL_TEMPLATE')
eaeLocalFileTemplate = os.getenv('EAE_LOCAL_FILE_TEMPLATE')

# remove all aes and eae input files and joined files; 
# we don't want to remove Connie's published file, so we do separately
cmd = 'rm %s/*.aes.*' % inputDir
print cmd
rc = os.system(cmd)
if rc != 0:
    msg = 'rm cmd failed: %s%s' % (cmd, CRT)
    fpDiag.write(msg)

cmd = 'rm %s/*.eae.*' % inputDir
print cmd
rc = os.system(cmd)
if rc != 0:
    msg = 'rm cmd failed: %s%s' % (cmd, CRT)
    fpDiag.write(msg)

cmd = 'rm %s/*.joined.*' % inputDir
print cmd
rc = os.system(cmd)
if rc != 0:
    msg = 'rm cmd failed: %s%s' % (cmd, CRT)
    fpDiag.write(msg)

# number of files unable to be downloaded
errorCt = 0

fpCur.write('%sFiles unable to be downloaded%s' % (CRT, CRT))
fpCur.write('------------------------------%s' % CRT)

# parse the input file and download files for each experiment
for line in fpInfile.readlines():
    expID = string.strip(string.split(line)[0])
    fpDiag.write('Download files for experiment ID: %s%s' % (expID, CRT))

    # AES
    aesFile = aesLocalFileTemplate % expID 
    fpAes = open(aesFile, 'w')
    aesURL = aesTemplate % (expID, expID)
    msg = ''
    # --retry-max-time 0, don't time out retries
    #--max-time - max time in seconds to allow the whole operation to take.
    stdout, stderr, statusCode = runCommand.runCommand("curl --max-time 10 --retry 5 --retry-delay 5 --retry-max-time 0 '%s'" % aesURL)
    #print 'statusCode: %s' % statusCode
    #print 'stderr: %s' % stderr
    if statusCode != 0:
        msg = '%s stderr: %s%s' % (aesURL, stderr, CRT)
	errorCt += 1
        fpDiag.write(msg)
        fpCur.write(msg)
    else:
	fpAes.write(stdout)

    # EAE
    eaeFile = eaeLocalFileTemplate % expID 
    fpEae = open(eaeFile, 'w')
    eaeURL = eaeTemplate % expID
    msg = ''
    #stdout, stderr, statusCode = runCommand.runCommand("curl --compressed --max-time 20 --retry 5 --retry-delay 5 --retry-max-time 0 '%s'" % eaeURL)
    stdout, stderr, statusCode = runCommand.runCommand("curl -v --retry 5 --retry-delay 5 '%s'" % eaeURL)
    print '%s%s statusCode: %s%s' % (CRT, expID, statusCode, CRT)
    print 'stderr: %s%s' % (stderr, CRT)
    if statusCode != 0:
        msg = '%s stderr: %s%s' % (eaeURL, stderr, CRT)
	errorCt += 1
	fpDiag.write(msg)
        fpCur.write(msg)
    else:
        fpEae.write(stdout)

fpCur.write('%sTotal: %s%s' % (CRT, errorCt, CRT))
fpDiag.write('%sTotal files unable to be downloaded: %s%s' % (CRT, errorCt, CRT))
fpCur.close()
fpDiag.close()
fpInfile.close()

print '%s' % mgi_utils.date()
