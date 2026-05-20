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
where s.name = 'Baseline RNASeq Load Experiments'
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

select distinct s.accid from samples s;

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
from gxd_htsample_rnaseqset r, gxd_htsample_rnaseqcombined m, mrk_marker mm
where r._rnaseqset_key = m._rnaseqset_key
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

-- are there any experiments in the 1
select old._rnaseqcombined_key, oldset._experiment_key, newset._rnaseqset_key
from GXD_HTSample_RNASeqCombined old, GXD_HTSample_RNASeqSet oldset,
        GXD_HTSample_RNASeqSet newset
where old._CreatedBy_key = 1613
and old._rnaseqset_key = oldset._rnaseqset_key
and oldset._experiment_key = newset._experiment_key
and newset._createdby_key = 1673
;

select distinct _marker_key
from gxd_htsample_rnaseqcombined
where _marker_key in (79559, 51389, 193172, 54345, 63868, 457016, 59010, 446540, 426272, 633827, 633841, 462799, 50578, 13004, 58264, 193898, 460129, 324682, 55282, 462572, 50504, 57533, 633847, 196238, 31865, 53493, 384828, 324182, 80937, 11813, 14883, 322170, 306160, 787424, 197065, 461805, 331138, 426758, 426360, 198858, 43740, 326110, 198450, 322243, 101137, 59342, 305633, 322024, 185382, 639640, 645479, 47349, 970754, 52982, 193804, 13010, 62139, 457957, 680005, 308615, 103731, 306046, 463209, 323939, 426036, 633475, 119570, 85736, 322862, 463131, 93082, 460538, 200254, 36952, 191299, 329064, 77087, 83365, 193084, 58457, 196999, 458901, 321036, 62777, 971427, 31282, 10819, 324103, 198196, 57165, 101920, 195887, 27129, 1006839, 192952, 45063, 196951, 459565, 55108, 646219, 61680, 313703, 32465, 461578, 104546, 31159, 58216, 461955, 415020, 191973, 459463, 307455, 105444, 104749, 57247, 307292, 425817, 309155, 331850, 58412, 463165, 62442, 198581, 463017, 41615, 72997, 9724, 306809, 194936, 101389, 46669, 415085, 158981, 198529, 51834, 429091, 57246, 600676, 444261, 25551, 382381, 426427, 307935, 972761, 56025, 461241, 667407, 192905, 12935, 324762, 195920, 101035, 459649, 193872, 460122)
;

-- fe check
--select count(e.*)
--from expression_ht_consolidated_sample e,
--expression_ht_consolidated_sample_measurement c
--where e.experiment_key in (17734,17831)
--and e.consolidated_sample_key = c.consolidated_sample_key
--;

EOSQL

date |tee -a $LOG

