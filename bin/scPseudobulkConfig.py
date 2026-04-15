import logging as log

log.basicConfig(level=log.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

class SCPseudobulkConfig:

    PSEUDOBULK_ID_SEX_2_group = 1
    PSEUDOBULK_ID_SEX_4_group = 2

    sex_2_group_query ='''
        WITH data AS (
            SELECT DISTINCT ON (individual)
                sex,
                _sample_key,
                individual,
                sex AS replicateKey
            FROM gxd_htsample
            WHERE _Experiment_key = 99680 AND individual NOT in('3_57_F', '3_11_M')
        )
        SELECT replicateKey,
            string_agg(DISTINCT _sample_key::text, ',') _sample_keys, count(*) num_of_samples,
            string_agg(DISTINCT (_sample_key::text || '|' || individual), ',') sample_key_names
        FROM data
        GROUP BY replicateKey
    '''
    sex_2_group_samples = {
        "name": "Option_A",
        "joined_file_name": "sex_7_samples.csv",
        "replicate_query": sex_2_group_query
    }

    sex_4_group_query ='''
        WITH data AS (
            SELECT DISTINCT ON (bioreplicate_name)
                sex,
                _sample_key,
                CASE
                    WHEN individual IN ('3_38_F', '3_56_F') THEN '3_38_F_and_3_56_F'
                    WHEN individual IN ('3_8_M', '3_9_M') THEN '3_8_M_and_3_9_M'
                    ELSE individual
                END AS bioreplicate_name,
				sex AS replicateKey
            FROM gxd_htsample
            WHERE _Experiment_key = 99680 AND individual NOT in('3_57_F', '3_11_M')
        )
        SELECT replicateKey,
            string_agg(DISTINCT _sample_key::text, ',') _sample_keys, count(*) num_of_samples,
            string_agg(DISTINCT (_sample_key::text || '|' || bioreplicate_name), ',') sample_key_names
        FROM data
        GROUP BY replicateKey
    '''
    sex_4_group_samples = {
        "name": "Option_B",
        "joined_file_name": "sex_4_samples.csv",
        "replicate_query": sex_4_group_query
    }

    PSEUDOBULK_LIST = {
        PSEUDOBULK_ID_SEX_2_group: sex_2_group_samples,
        PSEUDOBULK_ID_SEX_4_group: sex_4_group_samples
    }

    @staticmethod
    def get_pseudobulk_id(pseudobulk_id):
        return SCPseudobulkConfig.PSEUDOBULK_LIST.get(pseudobulk_id) 