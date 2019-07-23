#!/bin/sh
#
#  rnaseqQC.sh
###########################################################################
#
#  Purpose:
#
#      This script does sanity checks and is a wrapper around the process 
#	that does QC checks for the RNA-Seq load
#
#  Usage:
#
#      rnaseqQC.sh  filename  
#
#      where
#          filename = full path to the input file
#
#  Env Vars:
#
#      See the configuration file
#
#  Inputs:
#	Experiment ID file
#
#  Outputs:
#
#      - sanity report for the input file.
#      - QC report for the input file 	
#      - Log file (${QC_LOGFILE})
#
#  Exit Codes:
#
#      0:  Successful completion
#      1:  Sanity error (prompt to view sanity report)
#      2:  Unexpected error occured running rnaseqQC.py (prompt to view log)
#      3:  QC errors (prompt to view qc report)
#
#  Assumes:  Nothing
#
#  Implementation:
#
#      This script will perform following steps:
#
#      1) Validate the arguments to the script
#      2) Validate & source the configuration files to establish the environment
#      3) Verify that the input file exists
#      4) Update path to sanity/QC reports if this is not a 'live' run 
#	     i.e. curators running the scripts 
#      5) Initialize the log file
#      6) Call rnaseqQC.py to generate the QC report
#
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
CURRENTDIR=`pwd`
BINDIR=`dirname $0`

CONFIG=`cd ${BINDIR}/..; pwd`/rnaseqload.config
USAGE='Usage: rnaseqQC.sh  filename'

# set LIVE_RUN  to sanity/QC check only as the default
LIVE_RUN=0; export LIVE_RUN

#
# Make sure an input file was passed to the script. If the optional "live"
# argument is given, that means that the output files are located in the
# /data/loads/... directory, not in the current directory.
#
if [ $# -eq 1 ]
then
    INPUT_FILE=$1
elif [ $# -eq 2 -a "$2" = "live" ]
then
    INPUT_FILE=$1
    LIVE_RUN=1
else
    echo ${USAGE}; exit 1
fi

#
# Make sure the configuration file exists and source it.
#
if [ -f ${CONFIG} ]
then
    . ${CONFIG}
else
    echo "Missing configuration file: ${CONFIG}"
    exit 1
fi

#
# If the QC check is being run by a curator, the mgd_dbo password needs to
# be in a password file in their HOME directory because they won't have
# permission to read the password file in the pgdbutilities product.
#
if [ "${USER}" != "mgiadmin" ]
then
    PGPASSFILE=$HOME/.pgpass
fi

#
# If this is not a "live" run, the output, log and report files should reside
# in the current directory, so override the default settings.
#
if [ ${LIVE_RUN} -eq 0 ]
then
	SANITY_RPT=${CURRENTDIR}/`basename ${SANITY_RPT}`
	QC_RPT=${CURRENTDIR}/`basename ${QC_RPT}`
	QC_LOGFILE=${CURRENTDIR}/`basename ${QC_LOGFILE}`

fi

#
# Initialize the log file.
#
LOG=${QC_LOGFILE}
rm -rf ${LOG}
touch ${LOG}

#
# Convert the input file into a QC-ready version that can be used to run
# the sanity/QC reports against.
#
dos2unix ${INPUT_FILE} ${INPUT_FILE} 2>/dev/null

#
# FUNCTION: Check for lines with missing columns in input file and
#           write the line numbers to the sanity report.
#
checkColumns ()
{
    FILE=$1         # The input file to check
    REPORT=$2       # The sanity report to write to
    NUM_COLUMNS=$3  # The number of columns expected in each input record
    echo "" >> ${REPORT}
    echo "Lines With Missing Columns or Data" >> ${REPORT}
    echo "-----------------------------------" >> ${REPORT}
    echo " ${RNASEQLOAD}/bin/checkColumns.py ${FILE} ${NUM_COLUMNS} > ${TMP_FILE}"
    ${RNASEQLOAD}/bin/checkColumns.py ${FILE} ${NUM_COLUMNS} > ${TMP_FILE}
    cat ${TMP_FILE} >> ${REPORT}
    if [ `cat ${TMP_FILE} | wc -l` -eq 0 ]
    then
        return 0
    else
        return 1
    fi
}

#
# Initialize the report file(s) to make sure the current user can write to them.
#
RPT_LIST="${SANITY_RPT}"

for i in ${RPT_LIST}
do
    rm -f $i; >$i
done

#
# Create a temporary file and make sure it is removed when this script
# terminates.
#
TMP_FILE=/tmp/`basename $0`.$$
trap "rm -f ${TMP_FILE}" 0 1 2 15

#
# FUNCTION: Check an input file to make sure it has a minimum number of lines.
#
checkLineCount ()
{
    FILE=$1        # The input file to check
    REPORT=$2      # The sanity report to write to
    NUM_LINES=$3   # The minimum number of lines expected in the input file

    COUNT=`cat ${FILE} | wc -l | sed 's/ //g'`
    if [ ${COUNT} -lt ${NUM_LINES} ]
    then
        return 1
    else
        return 0
    fi
}

#
# FUNCTION: Check for duplicate lines in an input file and write the lines
#           to the sanity report.
#
checkDupLines ()
{
    FILE=$1    # The input file to check
    REPORT=$2  # The sanity report to write to

    echo "" >> ${REPORT}
    echo "Duplicate Lines" >> ${REPORT}
    echo "---------------" >> ${REPORT}
    sort ${FILE} | uniq -d > ${TMP_FILE}
    cat ${TMP_FILE} >> ${REPORT}
    if [ `cat ${TMP_FILE} | wc -l` -eq 0 ]
    then
        return 0
    else
        return 1
    fi
}

echo "" >> ${LOG}
date >> ${LOG}
echo "Run sanity checks on the input file" >> ${LOG}
SANITY_ERROR=0

#
# Make sure the input files exist (regular file or symbolic link).
#
echo "input file exists INPUT_FILE: ${INPUT_FILE}"
if [ "`ls -L ${INPUT_FILE} 2>/dev/null`" = "" ]
then
    echo "" | tee -a ${LOG}
    echo "Input file does not exist: ${INPUT_FILE}" | tee -a ${LOG}
    echo "" | tee -a ${LOG}
    exit 1
fi

checkLineCount ${INPUT_FILE} ${SANITY_RPT} ${MIN_LINES}

if [ $? -ne 0 ]
then
    echo "" | tee -a ${LOG}
    echo "Input file has less than min number lines:  ${MIN_LINES} ${INPUT_FILE}" | tee -a ${LOG}
    echo "" | tee -a ${LOG}
    exit 1
fi

checkColumns ${INPUT_FILE} ${SANITY_RPT} ${NUM_COLUMNS}
if [ $? -ne 0 ]
then
    SANITY_ERROR=1
fi

checkDupLines ${INPUT_FILE} ${SANITY_RPT}
if [ $? -ne 0 ]
then
    SANITY_ERROR=1
fi

if [ ${SANITY_ERROR} -ne 0 ]
then
    if [ ${LIVE_RUN} -eq 0 ]
    then
	echo "Sanity errors detected. See ${SANITY_RPT}" | tee -a ${LOG}
    fi
    exit 1
fi


#
# Generate the QC reports.
#
echo "" >> ${LOG}
date >> ${LOG}
echo "Generate the QC reports" >> ${LOG}
echo "rnaseqQC.py  INPUT_FILE: ${INPUT_FILE}"
{ ${RNASEQLOAD}/bin/rnaseqQC.py ${INPUT_FILE} 2>&1; echo $? > ${TMP_FILE}; } >> ${LOG}

if [ `cat ${TMP_FILE}` -eq 1 ]
then
    echo "An error occurred while generating the QC reports"
    echo "See log file (${LOG})"
    RC=2
elif [ `cat ${TMP_FILE}` -eq 2 ]
then
    if [ ${LIVE_RUN} -eq 0 ]
    then
	echo ""
	echo "QC errors detected. See ${QC_RPT} " | tee -a ${LOG}
    fi
    RC=3
else
    if [ ${LIVE_RUN} -eq 0 ]
    then
	echo "No QC errors detected"
    fi
    RC=0
fi

echo "" >> ${LOG}
date >> ${LOG}
echo "Finished running QC checks on the input file" >> ${LOG}

exit ${RC}
