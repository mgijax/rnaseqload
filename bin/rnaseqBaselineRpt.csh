#!/bin/csh -f

if ( ${?MGICONFIG} == 0 ) then
        setenv MGICONFIG /usr/local/mgi/live/mgiconfig
endif

source ${MGICONFIG}/master.config.csh

cd `dirname $0`

setenv LOG $0.log
rnamem -rf $LOG
touch $LOG
 
date | tee -a $LOG
 
cat - <<EOSQL | ${PG_DBUTILS}/bin/doisql.csh $0 | tee -a $LOG

select a.accid, a._object_key as _experiment_key, m.*
into temp table experiments
from MGI_Set s, MGI_SetMember m , ACC_Accession a
where s.name = 'Baseline RNASeq Load Experiment'
and s._set_key = m._set_key
and s._mgitype_key = a._mgitype_key
and m._object_key = a._object_key
and a._logicaldb_key = 189
and a.preferred = 1
;

select e.*, hts._sample_key, hts.name
into temp table samples
from experiments e, gxd_htsample hts
where hts._genotype_key != 90560
and hts._relevance_key = 20475450
and hts._experiment_key = e._experiment_key
;

--select a.accid, count(rnaset._rnaseqset_key)
--from experiments a, GXD_HTSample_RNASeqSet rnaset, GXD_HTSample_RNASeqSetMember rnamem, VOC_Term t
--where a._experiment_key = rnaset._experiment_key
--and rnaset._rnaseqset_key = rnamem._rnaseqset_key
--and rnaset._sex_key = t._Term_key
--and a.accid = 'E-ERAD-352'
--group by a.accid
--order by a.accid
--;

--select a.accid
--from experiments a
--where not exists (select 1 from GXD_HTSample_RNASeqSet rnaset where a._object_key = rnaset._experiment_key)
--order by a.accid
--;

--select distinct a.accid
--from experiments a, gxd_htsample hts
--where hts._genotype_key != 90560 --J:DO
--and hts._relevance_key = 20475450 --Yes
--and not exists (select 1 from GXD_HTSample_RNASeqSet rnaset where a._object_key = rnaset._experiment_key)
--order by a.accid
--;

--select distinct a.accid,  hts._genotype_key, hts._relevance_key
--from experiments a, gxd_htsample hts
--where hts._genotype_key != 90560
--and hts._relevance_key = 20475450 --Yes
--and a.accid = 'E-MTAB-7637'
--order by a.accid
--;

--select distinct a.accid
--from experiments a, GXD_HTSample_RNASeqSet rnaset, GXD_HTSample_RNASeqSetMember rnamem, VOC_Term t
--where a._experiment_key = rnaset._experiment_key
--and rnaset._rnaseqset_key = rnamem._rnaseqset_key
--and rnaset._sex_key = t._Term_key
--and exists (select 1 from GXD_HTSample_RNASeqCombined ss where rnamem._Sample_key = ss._Sample_key)
--order by a.accid, rnamem._Sample_key
--;

--select s.accid, rm._sample_key, rm._rnaseqsetmember_key, rs._rnaseqset_key
--from samples s, gxd_htsample_rnaseqset rs, gxd_htsample_rnaseqsetmember rm
--where s._sample_key = rm._sample_key
--and rm._rnaseqset_key = rs._rnaseqset_key
--and s.accid = 'E-GEOD-55966'
--order by s.accid
--;

-- number of bioreplicates
select distinct s.accid, rm._rnaseqset_key, rs.groupset, count(rm._rnaseqsetmember_key)
from samples s, gxd_htsample_rnaseqsetmember rm, gxd_htsample_rnaseqset rs
where s._sample_key = rm._sample_key
and rm._rnaseqset_key = rs._rnaseqset_key
--and s.accid in ('E-ERAD-401')
group by s.accid, rm._rnaseqset_key, rs.groupset
order by s.accid
;

select a.accid, count(a._object_key)
from acc_accession a
where a._mgitype_key = 2
and a._logicaldb_key = 60
and a.preferred = 1
group by a.accid having count(a._object_key) > 1
;

WITH ensembls AS (
select a._object_key, m.symbol, count(a.accid)
from acc_accession a, mrk_marker m
where a._mgitype_key = 2
and a._logicaldb_key = 60
and a.preferred = 1
and a._object_key = m._marker_key
group by a._object_key, m.symbol having count(a.accid) > 1
)
select a.accid, e._object_key, e.symbol
from ensembls e, acc_accession a
where e._object_key = a._object_key
and a._mgitype_key = 2
and a._logicaldb_key = 60
and a.preferred = 1
order by e.symbol
;

select count(*) from gxd_htsample_rnaseqset where _createdby_key = 1673;
select count(*) from gxd_htsample_rnaseqsetmember where _createdby_key = 1673;
select count(*) from gxd_htsample_rnaseqcombined where _createdby_key = 1673;
select _createdby_key, count(*) from gxd_htsample_rnaseqset group by _createdby_key;
select _createdby_key, count(*) from gxd_htsample_rnaseqcombined group by _createdby_key;

select r._experiment_key, count(*) as numResults
into temp table a
from gxd_htsample_rnaseqset_cache rc, gxd_htsample_rnaseqset r, gxd_htsample_rnaseqcombined m, mrk_marker mm
where rc._rnaseqset_key = r._rnaseqset_key
and rc._rnaseqcombined_key = m._rnaseqcombined_key
and m._marker_key = mm._marker_key
and mm._marker_status_key = 1
group by r._experiment_key
;

select r._experiment_key, count(*) as numResults
into temp table b
from gxd_htsample_rnaseqset r, gxd_htsample_rnaseqcombined m, mrk_marker mm
where r._rnaseqset_key = m._rnaseqset_key
and m._marker_key = mm._marker_key
and mm._marker_status_key = 1
group by r._experiment_key
;

EOSQL

date |tee -a $LOG

