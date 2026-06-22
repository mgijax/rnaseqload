import os
import pandas as pd
import pg_db
import os, glob

class CsvFilter:
    def __init__(self, data_dir):
        self.data_dir = data_dir
        self.sample_columns = [
            "3_8_M",
            "3_9_M",
            "3_10_M",
            "3_11_M"
        ]

    def remove_all_null_rows(self, input_file, output_file):
        input_path = os.path.join(self.data_dir, input_file)
        output_path = os.path.join(self.data_dir, output_file)

        df = pd.read_csv(input_path)

        # remove rows where all sample columns are NULL
        #df = df.dropna(subset=self.sample_columns, how="all")

        # replace remaining NULL values with zero
        #df[self.sample_columns] = df[self.sample_columns].fillna(0)

        # remove rows where ANY sample column is NULL
        df = df.dropna(subset=self.sample_columns, how="any")

        # replace NULL with row average
        # df[self.sample_columns] = df[self.sample_columns].apply(
        #     lambda row: row.fillna(row.mean()),
        #     axis=1
        # )

        df.to_csv(output_path, index=False)

        print(f"Output written to {output_path}")
        print(f"Rows remaining: {len(df)}")

    def merge_columns(self, input_file, output_file):
        input_path = os.path.join(self.data_dir, input_file)
        output_path = os.path.join(self.data_dir, output_file)

        df = pd.read_csv(input_path)
    # fill NULL with row mean
        # df[self.sample_columns] = df[self.sample_columns].apply(
        #     lambda row: row.fillna(row.mean()),
        #     axis=1
        # )

        # merge columns by sum
        df["3_8_M_n_3_9_M"] = df["3_8_M"] + df["3_9_M"]
        df["3_10_M_n_3_11_M"] = df["3_10_M"] + df["3_11_M"]

        # optional: drop original columns
        df = df.drop(columns=self.sample_columns)

        # write output
        df.to_csv(output_path, index=False)

        print(f"Output written to {output_path}")
        print(f"Rows remaining: {len(df)}")        




class CsvToPgDB:
    def __init__(self, data_dir, table_name):
        self.data_dir = data_dir
        self.table_name = table_name

    def transform(self, input_file, temp_file):
        input_path = os.path.join(self.data_dir, input_file)
        temp_path = os.path.join(self.data_dir, temp_file)

        df = pd.read_csv(input_path)

        # wide -> long
        value_cols = [c for c in df.columns if c != "gene"]

        df_long = df.melt(
            id_vars=["gene"],
            value_vars=value_cols,
            var_name="individual",
            value_name="count_sum"
        )

        df_long["bioreplicate_name"] = df_long["individual"]

        # write TSV for COPY
        df_long.to_csv(temp_path, sep="\t", index=False, header=False)

        return temp_path

    def load(self, temp_file):
        temp_path = os.path.join(self.data_dir, temp_file)

        with open(temp_path, "r") as f:
            pg_db.executeCopyFrom(
                f,
                self.table_name,
                sep="\t",
                null=r"\N"
            )

        pg_db.commit()

        print(f"Loaded into table: {self.table_name}")

    def create_table(self):
        sql = f"""
        DROP TABLE IF EXISTS {self.table_name};
        CREATE TABLE IF NOT EXISTS {self.table_name}
        (
            gene text COLLATE pg_catalog."default",
            individual character varying COLLATE pg_catalog."default",
            count_sum real,
            bioreplicate_name character varying COLLATE pg_catalog."default"
        )
        """
        pg_db.sql(sql)
        pg_db.commit()
        print(f"Table ensured: {self.table_name}")        

    def run(self, input_file):
        self.create_table()
        temp_file = "tmp_long.tsv"
        self.transform(input_file, temp_file)
        self.load(temp_file)


class DataDiff:

    def __init__(self, data_dir):
        self.files = glob.glob(os.path.join(data_dir, "*.csv"))
        self.df = None


    def run(self, out="data_diff.csv"):

        dfs = []
        for f in self.files:
            tag = os.path.splitext(os.path.basename(f))[0]

            df = pd.read_csv(f)
            df = df.rename(columns=lambda c: "gene" if c=="gene"
                           else f"{tag}_{c.replace('_replicate_group_avg','')}")

            dfs.append(df)

        df = dfs[0]
        for x in dfs[1:]:
            df = df.merge(x, on="gene", how="outer")

        colsA = [c for c in df if "_A_" in c]
        colsB = [c for c in df if "_B_" in c]

        same = lambda r, c: int(r[c].dropna().nunique() == 1)

        df["diffA"] = df.apply(lambda r: same(r, colsA), axis=1)
        df["diffB"] = df.apply(lambda r: same(r, colsB), axis=1)

        df = df[["gene"] + colsA + colsB +
                [c for c in df if c not in (["gene"]+colsA+colsB)]]

        df.to_csv(out, index=False)
        self.df = df
        return df


if __name__ == "__main__":
    # csv_filter = CsvFilter("/data/loads/liangh/rnaseqload/input")
    # csv_filter.remove_all_null_rows(
    #     "Lung_A_with_null.csv",
    #     "Lung_A_remove_null.csv"
    # )

    # csv_filter = CsvFilter("/data/loads/liangh/rnaseqload/input")
    # csv_filter.merge_columns(
    #     "Lung_A_metabimpute_rf.csv",
    #     "Lung_B_metabimpute_rf.csv"
    # )
    


    # pg_db.set_sqlLogin(
    #     user="xxx",
    #     password="xxx",
    #     server="bhmgidb06ld.jax.org",
    #     database="liangh_pub2"
    # )
    # pg_db.useOneConnection(True)

    # pg_db.sql("select 1")
    # print("DB connection OK")


    # files = ["Lung_A_average.csv","Lung_A_eucknn.csv","Lung_A_metabimpute_rf.csv",
    #          "Lung_A_original.csv","Lung_A_remove_null.csv","Lung_B_average.csv",
    #          "Lung_B_eucknn.csv","Lung_B_metabimpute_rf.csv","Lung_B_original.csv","Lung_B_remove_null.csv"]
    # tables = ["tm_im_lung_average_a","tm_im_lung_eucknn_a","tm_im_lung_metabimpute_rf_a",
    #         "tm_im_lung_original_a","tm_im_lung_remove_null_a","tm_im_lung_average_b",
    #         "tm_im_lung_eucknn_b","tm_im_lung_metabimpute_rf_b","tm_im_lung_original_b","tm_im_lung_remove_null_b"]
    # for i, file in enumerate(files):
    #     tableName = tables[i]
    #     loader = CsvToPgDB(
    #         data_dir="/data/loads/liangh/rnaseqload/input",
    #         table_name=tableName
    #     )
    #     loader.run(file)

    d = DataDiff("/data/loads/liangh/rnaseqload/imput")
    result = d.run("/data/loads/liangh/rnaseqload/imput/merged_diff.csv")
    print(result.head())    