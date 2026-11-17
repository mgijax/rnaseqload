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
import shutil

class MatchGenes:

    def __init__(self, tpm_dir, mismatch_data_dir, merged_dir, merged_file_name):
        self.tpm_dir = tpm_dir
        self.mismatch_data_dir = mismatch_data_dir
        self.merged_dir = merged_dir
        self.merged_file_name = merged_file_name

    def prepareData(self, gxd, mismatchRef, tm, trueRef):
        fileNames = ["average_TPM.csv", "metabimpute_rf_TPM.csv", "null_to_zero_TPM.csv", "saver_all_TPM.csv", "saver_imput_TPM.csv"]

        fromDir = os.path.join(self.tpm_dir, tm)
        toDir = os.path.join(self.mismatch_data_dir, self.getDirName(gxd, mismatchRef, tm, trueRef))
        # print (f'  {toDir}')
        os.makedirs(toDir, exist_ok=True)

        for filename in os.listdir(fromDir):
            if not filename in fileNames:
                continue
            src = os.path.join(fromDir, filename)
            dst = os.path.join(toDir, filename)
            if os.path.isfile(src):
                self.copy(src, dst)

        fileName = f"{mismatchRef}_TPM.csv"
        self.copy(os.path.join(self.tpm_dir, gxd, fileName), os.path.join(toDir, fileName))
        fileName = f"{trueRef}_TPM.csv"
        self.copy(os.path.join(self.tpm_dir, tm, fileName), os.path.join(toDir, fileName))

    def getDirName(self, gxd, mismatchRef, tm, trueRef):
        return f"GXD_{gxd}_{mismatchRef}__TM_{tm}_{trueRef}"

    def copy(self, src, dst):
        # print(f"    {src} => {dst}")
        shutil.copy2(src, dst)

    def mergeFile(self, gxd, mismatchRef, tm, trueRef):
        tissueDir = self.getDirName(gxd, mismatchRef, tm, trueRef)
        files = glob.glob(os.path.join(self.mismatch_data_dir, tissueDir, "*_TPM.csv"))
        dfs = []

        left_null_files = {"remove_null_TPM.csv"}
        left_zero_files = {"saver_all_TPM.csv", "saver_imput_TPM.csv"}

        for f in files:
            tag = os.path.splitext(os.path.basename(f))[0].replace("_TPM", "")

            df_tmp = pd.read_csv(f).rename(
                columns=lambda c:
                    "gene" if c == "gene"
                    else f"{tag}_A" if "_A_" in c
                    else f"{tag}_B" if "_B_" in c
                    else tag
            )

            print(f'    {os.path.basename(f)}: {len(df_tmp)} genes')

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
            "gene", "E-MTAB-2801", "E-MTAB-4035", "E-MTAB-5020", "E-MTAB-7279", "E-MTAB-8573", "E-MTAB-8838", 
            "E-GEOD-65775","E-GEOD-74747",
            "average_A","null_to_zero_A","metabimpute_rf_A",
            "remove_null_A","saver_all_A","saver_imput_A",
            "average_B","null_to_zero_B","metabimpute_rf_B","remove_null_B",
            "saver_all_B","saver_imput_B"
        ]

        sort_position = {col: i for i, col in enumerate(sort_order)}

        df = df[sorted(df.columns, key=lambda c: sort_position.get(c, len(sort_order)))]

        output_tissue_dir = os.path.join(self.merged_dir, tissueDir)
        os.makedirs(output_tissue_dir, exist_ok=True)

        out = os.path.join(self.merged_dir, tissueDir, self.merged_file_name)
        df.to_csv(out, index=False)
        print(f'    merged file: {len(df)} genes')
        return df

if __name__ == "__main__":  
    d = MatchGenes("/data/loads/liangh/rnaseqload/tpm_11_tissues", "/data/loads/liangh/rnaseqload/match_genes/",
                   "/data/loads/liangh/rnaseqload/match_genes_out/",  "matched_genes.csv")

    gxdToTM = {}
    gxdToTM["b_striatum"] = "b_cerebellum"
    gxdToTM["b_cerebellum"] = "b_cerebral"
    gxdToTM["b_cerebral"] = "b_hippocampus"
    gxdToTM["b_hippocampus"] = "b_striatum"
    gxdToTM["trachea"] = "fat"
    gxdToTM["fat"] = "heart"
    gxdToTM["heart"] = "intestine"
    gxdToTM["intestine"] = "lung"
    gxdToTM["lung"] = "spleen"
    gxdToTM["spleen"] = "thymus"
    gxdToTM["thymus"] = "trachea"

    refs = {}
    refs["b_cerebellum"] = ["E-MTAB-8838"]
    refs["b_cerebral"] = ["E-MTAB-7279", "E-MTAB-8838"]
    refs["b_hippocampus"] = ["E-MTAB-7279"]
    refs["b_striatum"] = ["E-MTAB-5020", "E-MTAB-7279"]
    refs["fat"] = ["E-MTAB-4035"]
    refs["heart"] = ["E-GEOD-65775", "E-GEOD-74747"]
    refs["intestine"] = ["E-MTAB-2801"]
    refs["lung"] = ["E-MTAB-4035","E-MTAB-8573"]
    refs["spleen"] = ["E-MTAB-4035"]
    refs["thymus"] = ["E-GEOD-74747"]
    refs["trachea"] = ["E-MTAB-8573"]

    for gxd, tm in gxdToTM.items():
        trueRefs = refs[tm]
        mismatchRefs = refs[gxd]
        print(f"{gxd} {mismatchRefs} (GXD): {tm} {trueRefs} (TM)")

        for mismatchRef in mismatchRefs:
            for trueRef in trueRefs:
                print(f"  {gxd} {mismatchRef} (GXD): {tm} {trueRef} (TM)")

                d.prepareData(gxd, mismatchRef, tm, trueRef)
                d.mergeFile(gxd, mismatchRef, tm, trueRef)

               
        # ba.process(tissue, "merged_tpm.csv", refs)
        # ba.splitFile(tissue, "merged_tpm.csv", refs)



    # tissues = ["fat_E-MTAB-4035_heart_E-GEOD-65775", "fat_E-MTAB-4035_heart_E-GEOD-74747"]
    # for tissue in tissues:
    #     result = d.run(tissue)
    #     print(f"\n{tissue}")
    #     print(result.head()) 

   # mismatchTissueRefs = {}
    # mismatchTissueRefs["b_cerebellum"] =  ["E-MTAB-5020", "E-MTAB-7279"]
    # mismatchTissueRefs["b_cerebral"] =["E-MTAB-8838"]
    # mismatchTissueRefs["b_hippocampus"] = ["E-MTAB-7279", "E-MTAB-8838"]
    # mismatchTissueRefs["b_striatum"] = ["E-MTAB-7279"]
    # mismatchTissueRefs["fat"] =  ["E-MTAB-8573"]
    # mismatchTissueRefs["heart"] = ["E-MTAB-4035"]
    # mismatchTissueRefs["intestine"] = ["E-GEOD-65775", "E-GEOD-74747"]
    # mismatchTissueRefs["lung"] = ["E-MTAB-2801"]
    # mismatchTissueRefs["spleen"] = ["E-MTAB-4035", "E-MTAB-8573"]
    # mismatchTissueRefs["thymus"] = ["E-MTAB-4035"]
    # mismatchTissueRefs["trachea"] = ["E-GEOD-74747"]
    #       