##########################################################################
#
# Purpose: Create pre-prcoessing files for each Experiment in Differential set
#
#	 DIFFRAW_INPUTDIR
#	 DIFFINPUTDIR
#
#   input files read from DIFFRAW_INPUTDIR
#	 DIFF_RAWCOUNTS_LOCAL_FILE_TEMPLATE
#        DIFF_SDRF_LOCAL_FILE_TEMPLATE 
#
#   generated pre-processing files created in DIFFINPUTDIR
#	 DIFF_RAWCOUNTS_PP_FILE_TEMPLATE
#        DIFF_SDRF_PP_FILE_TEMPLATE
#
# For each Experiment (xxx) from RNASeq MGI_Set
#   for Experiment file in DIFFRAW_INPUTDIR
#       prcoess the sdrf file (ppAESSdrfFile())
#           -> DIFFINPUTDIR/xxx.sdrf.txt
#       process the raw counts file (ppEAERawCountsFile())
#           -> DIFFINPUTDIR/xxx.raw_counts.txt
#           -> only process genes that appear in *all* raw_counts.txt files
#
###########################################################################

import os
import sys
import xml.etree.ElementTree as ET
import db
import mgi_utils

db.setTrace(True)

# Expression Atlas Experiment file Template - name of file stored locally
rawcountsTemplate = '%s' % os.getenv('DIFF_RAWCOUNTS_LOCAL_FILE_TEMPLATE')
rawcountsPPTemplate = '%s' % os.getenv('DIFF_RAWCOUNTS_PP_FILE_TEMPLATE')
sdrfTemplate = '%s' % os.getenv('DIFF_SDRF_LOCAL_FILE_TEMPLATE')
sdrfPPTemplate = '%s' % os.getenv('DIFF_SDRF_PP_FILE_TEMPLATE')
groupTemplate = '%s' % os.getenv('DIFF_GROUP_LOCAL_FILE_TEMPLATE')

# excluded genes
tsvGenesExcluded = os.getenv('DIFFRAW_INPUTDIR') + '/tsvGenesExcluded'
excludedGenes = []

# set of runs that exist in configuration for given experiment
rawRunConfigList = []

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
# input  : DIFF_RAWCOUNTS_LOCAL_FILE_TEMPLATE
#
# genes that do not appear in all raw-counts.txt files
# all genes must appear in all raw-counts.txt files
#
# input:
#   ensembl ID
#
def loadExcludedGenes():
    global excludedGenes

    print('in loadExcludedGenes()')

    #  read the tsvGenesExcluded
    print('tsvGenesExcluded: %s' % tsvGenesExcluded)

    # iterate thru the fpTsvGenesExcluded input file
    with open(tsvGenesExcluded, "r") as file:
        excludedGenes = file.read().strip().split()

    #print(excludedGenes)
    print('excludedGenes: ', str(len(excludedGenes)))

    return 0

# end loadExcludedGenes()

#
# input  : DIFF_GROUP_LOCAL_FILE_TEMPLATE
# output : rawRunConfigList
#
# list of run ids that exist in DIFF_GROUP_LOCAL_FILE_TEMPLATE
#
def ppEAEConfigurationFile(expID):
    global rawRunConfigList

    print('in ppEAEConfigurationFile(expID): %s' % expID)

    rawRunConfigList = []

    #  read the input file
    try:
        eaeFile = groupTemplate % expID
    except:
        print('skipping: missing -configuration.xml: %s' % (expID))
        return 1 # file does not exist
    print(eaeFile)

    #
    # eaeFile is in XML format
    # exclude <assay_group id="xxx_yyy" label...
    #
    # iterate thru the eaeFile xml file
    #        <assay_group id="g1" label="brown adipose tissue">
    #            <assay>ERR4193656</assay>
    #            <assay>ERR4193654</assay>
    #            <assay>ERR4193655</assay>
    #        </assay_group>

    print(eaeFile)
    tree = ET.parse(eaeFile)
    root = tree.getroot()
    assay_groups = root.findall('.//assay_group')
    for ag in assay_groups:
        id = ag.get('id')
        if str.find(id, '_') > -1:
            continue
        label = ag.get('label')
        runids = []
        for child in ag:
            #print(child.tag, child.text)
            runID = child.text
            rawRunConfigList.append(runID)

    print(rawRunConfigList)

    return 0

# end ppEAEConfigurationFile()

#
# input  : DIFF_SDRF_LOCAL_FILE_TEMPLATE
# output : DIFF_SDRF_PP_FILE_TEMPLATE
#
# format:
#   Source Name
#   > 1 ENA_RUN 
#
def ppAESSdrfFile(expID, objectKey):

    print('in ppAESSdrfFile(expID, object_key): %s,%s' % (expID, objectKey))

    enaRuns = []
    rawSampleRunList = {}

    #  read the input file
    aesFile = sdrfTemplate % expID
    try:
        fpAes = open(aesFile, 'r')
    except:
        print('skiping: missing .sdrf.txt file: %s' % (expID))
        return 1 # file does not exist

    #  create the output file
    ppFile = sdrfPPTemplate % expID
    try:
        fpPP = open(ppFile, 'w')
    except:
        return 1 # file does not exist

    # process the header line
    #
    headerList = str.split(fpAes.readline(), '\t')
    if headerList == ['']: # means file is empty
        print('skipping: missing header: %s' % (expID))
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
    if sourceSampleIDX == None and enaSampleIDX == None:
        print('skipping: missing Source Name/ENA_SAMPLE column: %s' % (expID))
        return 1
    if enaRunIDX == None:
        print('skipping: missing ENA_RUN column: %s' % (expID))
        return 1

    # iterate thru the fpAes input file
    for line in fpAes.readlines():

        tokens = str.split(line, '\t')

	# match MGI Sample to either sourceSampleIDX or enaSampleIDX
        sourceSample = str.strip(tokens[sourceSampleIDX])
        if sourceSample not in sampleInMGI:
            sourceSample = str.strip(tokens[enaSampleIDX])
        if sourceSample not in sampleInMGI:
            print('skipping: sample is not in MGI: %s, sourceSample = %s, enaSample = %s' % (expID, sourceSample, enaSample))
            continue

        enaRun = str.strip(tokens[enaRunIDX])

        # enaRun must exist in raw config file
        if enaRun not in rawRunConfigList:
            print('skipping: enaRun not found in rawRunConfigList: %s, %s' % (expID, enaRun))
            continue

        # skip duplicate enaRun
        if enaRun in enaRuns:
            continue
        enaRuns.append(enaRun)

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

        if sourceSample not in rawSampleRunList:
            rawSampleRunList[sourceSample] = []
        rawSampleRunList[sourceSample].append(enaRun)
    
    print(rawSampleRunList)
    for s in rawSampleRunList:
        fpPP.write('%s\t%s\n' % (s, ','.join(rawSampleRunList[s])))

    fpPP.close();
    fpAes.close();

    return 0

# end ppAESSdrfFile()

#
# input  : DIFF_RAWCOUNTS_LOCAL_FILE_TEMPLATE
# output : DIFF_RAWCOUNTS_PP_FILE_TEMPLATE
#
# input:
#   ensembl ID
#   marker symbol
#   samples
#
# format:
#   ensembl ID
#   marker key
#   marker symbol
#   sample tpm value
#
def ppEAERawCountsFile(expID):

    print('in ppEAERawCountsFile(expID): %s' % expID)

    #  read the input file
    eaeFile = rawcountsTemplate % expID
    print('eaeFile: %s' % eaeFile)
    try:
        fpEae = open(eaeFile, 'r')
    except:
        print('skipping: missing -rawcounts.tsv file: %s' % (expID))
        return 1 # file does not exist

    #  create the output file
    ppFile = rawcountsPPTemplate % expID
    try:
        fpPP = open(ppFile, 'w')
    except:
        return 1 # file does not exist

    #
    # read the header from fpEae and create header for fpPP
    # each "run" is in its own column
    # some raw-counts have duplicate headers & columns
    # example: E-MTAB-8161.raw-counts.txt
    # logic added:  uniqueHeaderList, dupSkip
    #

    headerList = str.split(fpEae.readline(), '\t')
    uniqueHeaderList = list(dict.fromkeys(headerList))
    columnList = []
    col = 1
    for h in uniqueHeaderList[2:]:

        h = str.strip(h)

        if len(headerList) > len(uniqueHeaderList):
            dupSkip = 1
        else:
            dupSkip = 0

        if dupSkip:
            col += 2
        else:
            col += 1

        if h not in rawRunConfigList:
            print('skipping: raw run not found in rawRunConfigList: %s, %s' % (expID, h))
            continue

        columnList.append(col)

    fpPP.write('ensembl_id\t_marker_key\tsymbol\t')
    fpPP.write('\t'.join(uniqueHeaderList))
    print(headerList)
    print(columnList)

    # iterate thru the fpEae input file
    for line in fpEae.readlines():

        tokens = str.split(line[:-1], '\t')
        ensemblID = str.strip(tokens[0])

        # if ensemblID is in excludedGenes, then skip
        if ensemblID in excludedGenes:
            #print('skipping: ensemblid is in the exclude set: %s, %s' % (expID, ensemblID))
            continue

        # if ensemblID is not in MGI, then set markerKey = 0
        # will handle this later during RAWCOUNTS processing
        if ensemblID in ensemblDict:
            markerKey = ensemblDict[ensemblID][0]['_object_key']
            markerSymbol = ensemblDict[ensemblID][0]['symbol']
        else:
            markerKey = 0
            markerSymbol = ''

        fpPP.write('%s\t%s\t%s' % (ensemblID, markerKey, markerSymbol))

        # for each column in this row
        for col in columnList:
            fpPP.write('\t' + str.strip(tokens[col]))
        fpPP.write('\n')

    fpEae.close();
    fpPP.close();

    return 0

# end ppEAERawCountsFile()

#
# pre processing
#   read EAE rawcounts, AES sdrf files
#   generate rawcounts output file, sdrf output file
#
# inputs:
#    DIFF_RAWCOUNTS_LOCAL_FILE_TEMPLATE
#    DIFF_SDRF_LOCAL_FILE_TEMPLATE 
#
# outputs:
#    DIFF_RAWCOUNTS_PP_FILE_TEMPLATE
#    DIFF_SDRF_PP_FILE_TEMPLATE
#
def process():
    global ensemblDict

    # set the excludedGenes list
    loadExcludedGenes()

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

    results = db.sql('''
        select a.accid, a._object_key
        from MGI_Set s, MGI_SetMember m , ACC_Accession a
        where s.name = 'RNASeq Load Experiments'
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

        #
        # order is important!
        #

        # process the eae/configuration file for this expID
        rc = ppEAEConfigurationFile(expID)
        if rc != 0:
            print('processing EAE configuration file returned rc %s, skipping file for %s' % (rc, expID))
            continue

        # process the aes/sdrf file for this expID
        rc = ppAESSdrfFile(expID, objectKey)
        if rc != 0:
            print('processing AES sdrf file returned rc %s, skipping file for %s' % (rc, expID))
            continue

        # process the eae/rawcounts file for this expID
        rc = ppEAERawCountsFile(expID)
        if rc != 0:
            print('processing EAE rawcounts file returned rc %s, skipping file for %s' % (rc, expID))
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
