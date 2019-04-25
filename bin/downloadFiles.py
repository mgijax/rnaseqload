#!/usr/bin/python
# /usr/bin/python is python2.6.* or 2.7.* depending on server. It is required for urllib2
##########################################################################
#
# Purpose:
#       Download ArrayExpress and Expression Atlas files by experiment
#
# Usage: downloadFiles.py
# Env Vars:
#	 1. 
#
# Inputs:
#	1. ${INPUTFILE} - Connie's file of experiment IDs
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

import sys
import os
import mgi_utils
import string
import db
import urllib2

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
TAB= '\t'
CRT = '\n'
SPACE = ' '

# ArrayExpress Sample File URL  Template
aesTemplate = os.getenv('AES_URL_TEMPLATE')

# Expression Atlas Experiment file URL Templage
eaeTemplate  = os.getenv('EAE_URL_TEMPLATE')

# wget commands

# ArrayExpress
# 'wget --no-check-certificate --tries=10 -nd -o E-MTAB-7277.tpm.htseq2.tsv.log -O E-MTAB-7277.tpm.htseq2.tsv ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/studies/arrayexpress/E-MTAB-7277/mus_musculus/genes.tpm.htseq2.tsv'

# Expression Atlas
# 'wget --no-check-certificate --tries=10 -nd -o log -O E-MTAB-7277.tpm.htseq2.tsv ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/rnaseq/studies/arrayexpress/E-MTAB-7277/mus_musculus/genes.tpm.htseq2.tsvdrf.txt'

wgetCmd = 'wget --no-check-certificate --tries=10 -nd -o %s -O %s %s'

# remove all aes and eae input files; we don't want to remove Connie's
#  published file
cmd = 'rm %s/*.aes.*' % inputDir
rc = os.system(cmd)
if rc != 0:
    print 'rm cmd failed: %s%s' % (cmd, CRT)
    fpDiag.write(cmd)
    fpCur.write(cmd)

cmd = 'rm %s/*.eae.*' % inputDir
rc = os.system(cmd)
if rc != 0:
    print 'rm cmd failed: %s%s' % (cmd, CRT)
    fpDiag.write(cmd)
    fpCur.write(cmd)

for line in fpInfile.readlines():
    expID = string.strip(string.split(line)[0])
    print 'ID: "%s"' % expID

    # AES
    aesLog = '%s/%s.aes.log' % (logDir, expID) 
    aesFile = '%s/%s.aes.txt' % (inputDir, expID)
    fpAes = open(aesFile, 'w')
    aesURL = aesTemplate % (expID, expID)
    #aesCmd = wgetCmd % (aesLog, aesFile, aesURL )
    #print aesCmd
    #rc = os.system(aesCmd)
    #if rc != 0:
    #    print 'aesCmd failed: %s' % aesCmd
    try:
	request = urllib2.Request(aesURL)
	result = urllib2.urlopen(request)
	fpAes.write(result.read())
	fpAes.close()
    except:
	msg = 'Unable to download %s%s' % (aesURL, CRT)
	fpDiag.write(msg)
	fpCur.write(msg)
    
    # EAE
    eaeLog = '%s/%s.eae.log' % (logDir, expID)
    eaeFile = '%s/%s.eae.txt' % (inputDir, expID)
    fpEae = open(eaeFile, 'w')
    eaeURL = eaeTemplate % expID
    #eaeCmd = wgetCmd % (eaeLog, eaeFile, eaeURL)
    #print eaeCmd
    #rc = os.system(eaeCmd)
    #if rc != 0:
    #    print 'eaeCmd failed: %s' % eaeCmd
    try:
	request = urllib2.Request(eaeURL)
	result = urllib2.urlopen(request)
	fpEae.write(result.read())
	fpEae.close()
    except:
	msg = 'Unable to download %s%s' % (eaeURL, CRT)
  	print msg
	fpDiag.write(msg)
	fpCur.write(msg)

fpCur.close()
fpDiag.close()
fpInfile.close()

print '%s' % mgi_utils.date()
