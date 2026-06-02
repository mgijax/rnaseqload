##########################################################################
#
# Purpose: Create pre-prcoessing files for each Experiment in Baseline set
#
#	 BASELINERAW_INPUTDIR
#	 BASELINEINPUTDIR
#
#   input files read from BASELINERAW_INPUTDIR
#	 BASELINE_TPMS_LOCAL_FILE_TEMPLATE
#    BASELINE_GROUP_LOCAL_FILE_TEMPLATE
#    BASELINE_SDRF_LOCAL_FILE_TEMPLATE 
#
#   generated pre-processing files created in BASELINEINPUTDIR
#	 BASELINE_TPMS_PP_FILE_TEMPLATE
#    BASELINE_GROUP_PP_FILE_TEMPLATE
#
# For each Experiment (xxx) from Baseline RNASeq MGI_Set
#   for Experiment file in BASELINERAW_INPUTDIR
#       process the tpms file (ppEAETpmsFile())
#           -> BASELINEINPUTDIR/xxx.tpms.txt
#       prcoess the sdrf file (ppAESSdrfFile())
#           -> runToSampleDict
#       process the configuration (ppEAEGroupFile())
#           -> BASELINEINPUTDIR/xxx.group.txt
#
###########################################################################

import os
import sys
import xml.etree.ElementTree as ET
import db
import mgi_utils

# Expression Atlas Experiment file Template - name of file stored locally
tpmsTemplate = '%s' % os.getenv('BASELINE_TPMS_LOCAL_FILE_TEMPLATE')
tpmsPPTemplate = '%s' % os.getenv('BASELINE_TPMS_PP_FILE_TEMPLATE')
groupTemplate = '%s' % os.getenv('BASELINE_GROUP_LOCAL_FILE_TEMPLATE')
groupPPTemplate = '%s' % os.getenv('BASELINE_GROUP_PP_FILE_TEMPLATE')
aesTemplate = '%s' % os.getenv('BASELINE_SDRF_LOCAL_FILE_TEMPLATE')

# unique set of raw samples
rawRunList = []

# which samples to this run belong to
# run -> sample
runToSampleDict = {}

#
# loads a lookup of samples in the db for the given experiment
#
def loadSamples(expID):
    global sampleInMGI
    sampleInMGI = []

    results = db.sql('''
        select hts.name, hts._Sample_key
        from GXD_HTSample hts, ACC_Accession a
        where hts._Experiment_key = a._Object_key
        and a._MGIType_key = 42 -- experiment
        and a._LogicalDB_key = 189 --ArrayExpress
        and a.preferred = 1
        and a.accID = '%s' 
        ''' % expID, 'auto')
    for r in results:
        key = str.strip(r['name'])
        sampleInMGI.append(key)

    return 0

# end loadSamples()

#
# input  : BASELINE_TPMS_LOCAL_FILE_TEMPLATE
# output : BASELINE_TPMS_PP_FILE_TEMPLATE
#
# format:
#   ensembm ID
#   marker key
#   marker symbol
#   each group (g1, g2, etc. value = 3rd value avg QN TPM)
#
def ppEAETpmsFile(expID):

    print('in ppEAETpmsFile(expID): %s' % expID)

    #  read the input file
    eaeFile = tpmsTemplate % expID
    print('eaeFile: %s' % eaeFile)
    try:
        fpEae = open(eaeFile, 'r')
    except:
        print('skipping: missing -tpms.tsv file: %s' % (expID))
        return 1 # file does not exist

    #  create the output file
    ppFile = tpmsPPTemplate % expID
    try:
        fpPP = open(ppFile, 'w')
    except:
        return 1 # file does not exist

    ensemblDict = {}
    results = db.sql('''
        select a.accid, a._object_key, m.symbol 
        from acc_accession a, mrk_marker m 
        where a._logicaldb_key = 60 
        and a._mgitype_key = 2 
        and a.preferred = 1
        and a._object_key = m._marker_key
        ''', 'auto')
    for r in results:
        key = r['accid']
        value = r
        if key not in ensemblDict:
            ensemblDict[key] = []
        ensemblDict[key].append(value)

    # read the header from fpEae and create header for fpPP
    headerList = str.split(fpEae.readline(), '\t')
    groupSet = []
    for h in headerList[2:]:
        groupSet.append(str.strip(h))
    fpPP.write('ensembl_id\t_marker_key\tsymbol\t' + '\t'.join(groupSet) + '\n')

    # iterate thru the fpEae input file
    for line in fpEae.readlines():

        tokens = str.split(line, '\t')
        ensemblID = str.strip(tokens[0])

        # if ensemblID is not in MGI, then set markerKey = 0
        # will handle this later during TPMS processing
        if ensemblID in ensemblDict:
            markerKey = ensemblDict[ensemblID][0]['_object_key']
            markerSymbol = ensemblDict[ensemblID][0]['symbol']
        else:
            markerKey = 0
            markerSymbol = ''

        fpPP.write('%s\t%s\t%s' % (ensemblID, markerKey, markerSymbol))

        for g in range(len(groupSet)):
            values = tokens[g+2].split(',')
            fpPP.write('\t' + values[2])

        fpPP.write('\n')

    fpEae.close();
    fpPP.close();

    return 0

# end ppEAETpmsFile()

#
# input  : BASELINE_SDRF_LOCAL_FILE_TEMPLATE
# output : runToSampleDict
#   store all Source Name's per ENA_RUN
#   E-ERAD-169.sdrf.txt : runToSampleDict['ERS223116'] = ['ERR323395', 'ERR323401']
#
# format:
#   Source Name
#   ENA_SAMPLE
#   ENA_RUN
#
def ppAESSdrfFile(expID, objectKey):

    global rawRunList, runToSampleDict

    print('in ppAESSdrfFile(expID): %s' % expID)

    #  read the input file
    aesFile = aesTemplate % expID
    try:
        fpAes = open(aesFile, 'r')
    except:
        print('skiping: missing .sdrf.txt file: %s' % (expID))
        return 1 # file does not exist

    # process the header line
    #
    headerList = str.split(fpAes.readline(), '\t')
    if headerList == ['']: # means file is empty
        print ('skipping: missing header: %s' % (expID))
        return 1

    # load sampleInMGI() for expID
    loadSamples(expID)

    # find the idx of the columns we want - they are not ordered
    sourceSampleIDX = None
    enaSampleIDX = None
    enaRunIDX = None
    for idx, colName in enumerate(headerList):
        colName = str.strip(colName)
        if str.find(colName, 'Source Name') != -1:
            sourceSampleIDX = idx
        elif str.find(colName, 'ENA_SAMPLE') != -1:
            enaSampleIDX = idx
        elif str.find(colName, 'ENA_RUN') != -1:
            enaRunIDX = idx
    if sourceSampleIDX == None:
        print('skipping: missing Source Name column: %s' % (expID))
        return 1
    if enaRunIDX == None:
        print('skipping: missing ENA_RUN column: %s' % (expID))
        return 1

    # iterate thru the fpEae input file
    for line in fpAes.readlines():

        tokens = str.split(line, '\t')

        sourceSample = None
        enaSample = None
        enaRun = None

        if sourceSampleIDX != None:
            sourceSample = str.strip(tokens[sourceSampleIDX])

        # not every sample file contais enaSampleIDX column
        if enaSampleIDX != None:
            enaSample = str.strip(tokens[enaSampleIDX])

        if enaRunIDX != None:
            enaRun = str.strip(tokens[enaRunIDX])

        # skip if this source is a duplicate; but don't report
        if enaRun in rawRunList:
            #print('skipping: enaRun already processed: %s,%s' % (expID, enaRun))
            continue

        if sourceSample not in sampleInMGI:
            if enaSample != None and enaSample in sampleInMGI:
                sourceSample = enaSample
            else:
                print('skipping: sample is not in MGI: %s, sourceSample = %s, enaSample = %s' % (expID, sourceSample, str(enaSample)))
                continue

        # if sourceSample exists in MGI, is genotype = J:DO (_genotype_key = 90560), 
        #   or Relevance != Yes (_relevance_key != 20475450), 
        # then skip
        ignoreResults = db.sql('''
            select * from GXD_HTSample where (_genotype_key = 90560 or _relevance_key != 20475450)
                and _experiment_key = %s and name = '%s' 
            ''' % (objectKey, sourceSample), 'auto')
        if len(ignoreResults) > 0:
            #print('skipping: sample is J:DO or Relevance != Yes')
            continue

        rawRunList.append(enaRun)

        key = enaRun
        value = sourceSample
        if key not in runToSampleDict:
            runToSampleDict[key] = []
        runToSampleDict[key].append(value)

    fpAes.close();

    return 0

# end ppAESSdrfFile()

#
# input  : BASELINE_GROUP_LOCAL_FILE_TEMPLATE
# output : BASELINE_GROUP_PP_FILE_TEMPLATE
#
# format:
# 	group ID 
# 	label 
# 	Run IDs
#   Sample IDs
#
def ppEAEGroupFile(expID):

    print('in ppEAEGroupFile(expID): %s' % expID)

    #  read the input file
    try:
        eaeFile = groupTemplate % expID
    except:
        print('skipping: missing -configuration.xml: %s' % (expID))
        return 1 # file does not exist

    #  create the output file
    ppFile = groupPPTemplate % expID
    try:
        fpPP = open(ppFile, 'w')
    except:
        return 1 # file does not exist

    #
    # eaeFile is in XML format
    #
    # iterate thru the eaeFile xml file
    #        <assay_group id="g1" label="brown adipose tissue">
    #            <assay>ERR4193656</assay>
    #            <assay>ERR4193654</assay>
    #            <assay>ERR4193655</assay>
    #        </assay_group>

    tree = ET.parse(eaeFile)
    root = tree.getroot()
    assay_groups = root.findall('.//assay_group')
    for ag in assay_groups:
        id = ag.get('id')
        label = ag.get('label')
        runids = []
        for child in ag:
            #print(child.tag, child.text)
            runID = child.text
            sampleID = 'missing'
            if runID in runToSampleDict:
                sampleID = runToSampleDict[runID][0]
            fpPP.write('%s\t%s\t%s\t%s\n' % (id, label, runID, sampleID))

    fpPP.close();

    return 0

# end ppEAEGroupFile()
 
#
# pre processing
#   read EAE tpms, group and AES sdrf files
#   generate tpms output file, group output file
#
# inputs:
#	 BASELINE_TPMS_LOCAL_FILE_TEMPLATE
#    BASELINE_GROUP_LOCAL_FILE_TEMPLATE
#    BASELINE_SDRF_LOCAL_FILE_TEMPLATE 
#
# outputs:
#	 BASELINE_TPMS_PP_FILE_TEMPLATE
#    BASELINE_GROUP_PP_FILE_TEMPLATE
#
def process():
    global rawRunList

    results = db.sql('''
        select a.accid, a._object_key
        from MGI_Set s, MGI_SetMember m , ACC_Accession a
        where s.name = 'Baseline RNASeq Load Experiments'
        and s._set_key = m._set_key
        and s._mgitype_key = a._mgitype_key
        and m._object_key = a._object_key
        and a._logicaldb_key = 189
        and a.preferred = 1
        ''', 'auto')

    #
    # for each expID in the MGI_Set:
    #
    for r in results:

        expID = str.strip(r['accid'])
        objectKey = r['_object_key']

        # process the eae/tpms file for this expID
        rc = ppEAETpmsFile(expID)
        if rc != 0:
            print('processing EAE tpms file returned rc %s, skipping file for %s' % (rc, expID))
            continue

        # process the aes/sdrf file for this expID to create the runToSampleDict
        rc = ppAESSdrfFile(expID, objectKey)
        if rc != 0:
            print('processing AES sdrf file returned rc %s, skipping file for %s' % (rc, expID))
            continue

        # process the eae/group file for this expID
        rc = ppEAEGroupFile(expID)
        if rc != 0:
            print('processing EAE group file returned rc %s, skipping file for %s' % (rc, expID))
            continue

    return 0

# end process()

#
# Main
#

print('start time: %s' %  mgi_utils.date())
if process() != 0:
     exit(1, 'Error in process()\n')
print('end time: %s' %  mgi_utils.date())
