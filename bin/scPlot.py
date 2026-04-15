import logging as log
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

log.basicConfig(level=log.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

class SCPlot:
    def __init__(self, exp_id):    
        self.exp_id = exp_id
        self.input_dir = os.getenv('INPUTDIR', '.')
        self.output_dir = "/usr/local/mgi/live/rnaseqload/data"

    def boxplot(self, title, csv_file_name):
        file_path = os.path.join(self.input_dir, csv_file_name)
        log.info(f"Reading file: {file_path}")

        df = pd.read_csv(file_path)
        if "gene" in df.columns:
            df = df.set_index("gene")

        df.boxplot(rot=45)
        plt.title(f"{title})")
        plt.ylabel("Expression")

        output_file = os.path.join(self.output_dir, f"{title.replace(' ', '_')}.png")
        plt.savefig(output_file, dpi=150)
        log.info(f"Plot saved to: {output_file}")

        # optional: still show if environment supports it
        # plt.show()

        plt.close()

    def plot_before_after(self, raw_file, qn_file, title, output_file):

        df_raw = pd.read_csv(os.path.join(self.input_dir, raw_file)).set_index("gene")
        df_qn = pd.read_csv(os.path.join(self.input_dir, qn_file)).set_index("gene")

        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        df_raw.boxplot(ax=axes[0])
        axes[0].set_title("Before QN")

        df_qn.boxplot(ax=axes[1])
        axes[1].set_title("After QN")

        # ✅ Main title
        fig.suptitle(f"{title} Before vs After Quantile Normalization", fontsize=14)

        output = os.path.join(self.output_dir, output_file)

        # leave space for suptitle
        plt.tight_layout(rect=[0, 0, 1, 0.95])

        plt.savefig(output, dpi=150)
        plt.close()

        log.info(f"Saved plot to: {output}")

    def plot_heatmap(self, qn_files, heatmap_file_name, top_n=50, title=None):
        import numpy as np
        import seaborn as sns
        import matplotlib.pyplot as plt
        import pandas as pd
        import os
        import logging as log

        # allow single file input
        if isinstance(qn_files, str):
            qn_files = [qn_files]

        dfs = []

        # ---- Load and prepare each file ----
        for qn_file in qn_files:
            qn_file_path = os.path.join(self.input_dir, qn_file)

            df = pd.read_csv(qn_file_path).set_index("gene")

            # add prefix (file name) to columns
            base = os.path.basename(qn_file).replace(".csv", "")
            df.columns = [f"{base}_{col}" for col in df.columns]

            dfs.append(df)

        # ---- Combine all files ----
        combined = pd.concat(dfs, axis=1)
        combined = combined.fillna(0)

        # ---- Clean column labels ----
        cleaned_cols = [c.split(")")[-1].strip("_") for c in combined.columns]

        # ⚠️ ensure uniqueness (important!)
        seen = {}
        final_cols = []
        for c in cleaned_cols:
            if c in seen:
                seen[c] += 1
                final_cols.append(f"{c}_{seen[c]}")
            else:
                seen[c] = 0
                final_cols.append(c)

        combined.columns = final_cols

        # ---- Select top variable genes ----
        top_var = combined.var(axis=1).sort_values(ascending=False).head(top_n).index
        df_sub = combined.loc[top_var]

        # ---- Log transform ----
        df_sub = np.log2(df_sub + 1)

        # ---- Plot ----
        g = sns.clustermap(
            df_sub,
            cmap="viridis",
            figsize=(12, 10),
            standard_scale=0
        )

        # rotate labels (important)
        plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha="right")

        # ---- Add title ----
        if title:
            g.fig.suptitle(title, fontsize=14)
        else:
            g.fig.suptitle("Combined Heatmap (Quantile Normalized Data)", fontsize=14)

        g.fig.subplots_adjust(top=0.92)

        # ---- Output ----
        output_file = os.path.join(self.output_dir, heatmap_file_name)
        g.savefig(output_file, dpi=150)
        plt.close()

        log.info(f"Saved: {output_file}")    

    def plot_distribution_compare(self, raw_file, qn_file, title, output_file, log_transform=True):
        # ---- Load data ----
        raw_path = os.path.join(self.input_dir, raw_file)
        qn_path = os.path.join(self.input_dir, qn_file)

        df_raw = pd.read_csv(raw_path).set_index("gene")
        df_qn = pd.read_csv(qn_path).set_index("gene")

        plt.figure(figsize=(10, 6))

        # ---- Plot RAW (per sample) ----
        for col in df_raw.columns:
            values = df_raw[col].values
            if log_transform:
                values = np.log2(values + 1)

            pd.Series(values).plot(
                kind="density",
                alpha=0.3,
                linestyle="--",
                label=f"RAW_{col}"
            )

        # ---- Plot QN (per sample) ----
        for col in df_qn.columns:
            values = df_qn[col].values
            if log_transform:
                values = np.log2(values + 1)

            pd.Series(values).plot(
                kind="density",
                alpha=0.8,
                linestyle="-",
                label=f"QN_{col}"
            )

        # ---- Labels ----
        plt.title(f"{title} Distribution: Raw vs Quantile Normalized")
        plt.xlabel("Expression (log2 scale)" if log_transform else "Expression")
        plt.ylabel("Density")

        plt.legend(fontsize=8, ncol=2)
        plt.tight_layout()

        # ---- Save ----
        output = os.path.join(self.output_dir, output_file)
        plt.savefig(output, dpi=150)
        plt.close()

        print(f"Saved: {output}")                    

def main():
    scPlot = SCPlot('') 
    scPlot.plot_before_after('Option_A_QN_input_female_(3_samples).csv', 'Option_A_QN_output_female_(3_samples).csv', 'Option A Female', 'Option_A_Female.png')
    scPlot.plot_before_after('Option_A_QN_input_male_(3_samples).csv', 'Option_A_QN_output_male_(3_samples).csv', 'Option A Male', 'Option_A_Male.png')
    scPlot.plot_before_after('Option_B_QN_input_female_(2_samples).csv', 'Option_B_QN_output_female_(2_samples).csv', 'Option B Female', 'Option_B_Female.png')
    scPlot.plot_before_after('Option_B_QN_input_male_(2_samples).csv', 'Option_B_QN_output_male_(2_samples).csv', 'Option B Male', 'Option_B_Male.png')

    scPlot.plot_heatmap(["Option_A_QN_output_female_(3_samples).csv", "Option_A_QN_output_male_(3_samples).csv"], 'Option_A_HeatMap.png',
                        top_n=1000,
                        title="Option A Heatmap (Quantile Normalized Data)")
    scPlot.plot_heatmap(["Option_B_QN_output_female_(2_samples).csv", "Option_B_QN_output_male_(2_samples).csv"], 'Option_B_HeatMap.png',
                        top_n=1000,
                        title="Option B Heatmap (Quantile Normalized Data)")

    scPlot.plot_distribution_compare('Option_A_QN_input_female_(3_samples).csv', 'Option_A_QN_output_female_(3_samples).csv', 'Option A Female', 'Option_A_Female_Distribution.png')
    scPlot.plot_distribution_compare('Option_A_QN_input_male_(3_samples).csv', 'Option_A_QN_output_male_(3_samples).csv', 'Option A Male', 'Option_A_Male_Distribution.png')
    scPlot.plot_distribution_compare('Option_B_QN_input_female_(2_samples).csv', 'Option_B_QN_output_female_(2_samples).csv', 'Option B Female', 'Option_B_Female_Distribution.png')
    scPlot.plot_distribution_compare('Option_B_QN_input_male_(2_samples).csv', 'Option_B_QN_output_male_(2_samples).csv', 'Option B Male', 'Option_B_Male_Distribution.png')
 

    










    # scPlot.boxplot('Raw for 2 Group Female (3 Samples)', 'sex_2_group_QN_input_female_3.csv') 
    # scPlot.boxplot('Boxplot for 2 Group Female (3 Samples)', 'sex_2_group_QN_output_female_3.csv') 
    # scPlot.boxplot('Boxplot for 2 Group Male (3 Samples)', 'sex_2_group_QN_output_male_3.csv') 

    # scPlot.boxplot('Boxplot for 4 Group Female 1 (2 Samples)', 'sex_4_group_QN_output_female_1_2.csv') 
    # scPlot.boxplot('Boxplot for 4 Group Female 2 (1 Samples)', 'sex_4_group_QN_output_female_2_1.csv') 
    # scPlot.boxplot('Boxplot for 4 Group Male 1 (2 Samples)', 'sex_4_group_QN_output_male_1_2.csv') 

if __name__ == "__main__":
    main()