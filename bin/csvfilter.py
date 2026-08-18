import os
import pandas as pd
import numpy as np
import pg_db
import os, glob
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations
import numpy as np
import pandas as pd

class CSVFileMerger:

    def __init__(self, input_dir, output_dir, merged_file_name):
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.merged_file_name = merged_file_name

    def run(self, tissue):
        files = glob.glob(os.path.join(self.input_dir, tissue, "*.csv"))
        dfs = []

        left_null_files = {"remove_null_TPM.csv"}
        left_zero_files = {"saver_all_TPM.csv", "saver_replace_null_only_TPM.csv"}

        for f in files:
            tag = os.path.splitext(os.path.basename(f))[0].replace("_TPM", "")

            df_tmp = pd.read_csv(f).rename(
                columns=lambda c:
                    "gene" if c == "gene"
                    else f"{tag}_A" if "_A_" in c
                    else f"{tag}_B" if "_B_" in c
                    else tag
            )

            dfs.append((os.path.basename(f), df_tmp))

        # Inner join all files except the three left-join files
        inner_dfs = [
            df for name, df in dfs
            if name not in left_null_files and name not in left_zero_files
        ]
        left_null_dfs = [df for name, df in dfs if name in left_null_files]
        left_zero_dfs = [df for name, df in dfs if name in left_zero_files]

        df = inner_dfs[0]

        for df_tmp in inner_dfs[1:]:
            df = df.merge(df_tmp, on="gene", how="inner")

        for df_tmp in left_null_dfs:
            df = df.merge(df_tmp, on="gene", how="left")

        for df_tmp in left_zero_dfs:
            cols = df_tmp.columns.drop("gene")
            df = df.merge(df_tmp, on="gene", how="left")
            df[cols] = df[cols].fillna(0) 

        sort_order = [
            "gene","E-MTAB-4035","E-MTAB-8573","E-GEOD-65775","E-GEOD-74747","E-MTAB-2801","E-MTAB-4035",
            "average_A","null_to_zero_A","metabimpute_rf_A",
            "remove_null_A","saver_all_A","saver_replace_null_only_A",
            "average_B","null_to_zero_B","metabimpute_rf_B","remove_null_B",
            "saver_all_B","saver_replace_null_only_B"
        ]

        sort_position = {col: i for i, col in enumerate(sort_order)}

        df = df[sorted(df.columns, key=lambda c: sort_position.get(c, len(sort_order)))]

        output_tissue_dir = os.path.join(self.output_dir, tissue)
        os.makedirs(output_tissue_dir, exist_ok=True)

        out = os.path.join(output_tissue_dir, self.merged_file_name)
        df.to_csv(out, index=False)

        return df

if __name__ == "__main__":
    d = CSVFileMerger("/data/loads/liangh/rnaseqload/tpm/", "/data/loads/liangh/rnaseqload/bland_altman", "merged_tpm.csv")
    tissues = ["heart", "intestine", "lung", "spleen"]
    for tissue in tissues:
        result = d.run(tissue)
        print(f"\n{tissue}")
        print(result.head()) 