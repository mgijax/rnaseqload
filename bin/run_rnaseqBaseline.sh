#!/bin/sh
#
# Purpose:
#	Run Baselinne Processing
#

cd `dirname $0`
LOG=`pwd`/rnaseqBaseline.log
rm -rf ${LOG}

CONFIG_LOAD=../rnaseqload.config

#
# verify & source the configuration file
#

if [ ! -r ${CONFIG_LOAD} ]
then
    echo "Cannot read configuration file: ${CONFIG_LOAD}"
    exit 1
fi

. ${CONFIG_LOAD}

rm -rf ${BASELINELOG}

cat - <<EOSQL | ${PG_DBUTILS}/bin/doisql.csh $0 | tee -a ${BASELINELOG}

select a.accid, a._object_key, gs.name, gs._sample_key
into temp toDelete
from ACC_Accession a, MGI_Set s, MGI_SetMember sm, GXD_HTSample gs
where s.name = 'Baseline RNASeq Load Experiment'
and s._Set_key = sm._Set_key
and sm._Object_key = a._Object_key
and a._MGIType_key = 42
and a._LogicalDB_key = 189
and a.preferred = 1
and a._Object_key = gs._Experiment_key
;

select * from toDelete;
--delete from gxd_htsample_rnaseq using toDelete where toDelete._sample_key = gxd_htsample_rnaseq._sample_key;


EOSQL

echo "Baseline Processing" 
${PYTHON} ${RNASEQLOAD}/bin/rnaseqBaseline.py >> ${>> ${BASELINELOG}

date | tee -a ${BASELINELOG}
