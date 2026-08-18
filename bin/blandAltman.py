import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


class BlandAltman:

    def __init__(self,dataDir):
        self.dataDir=dataDir

    def process(self, tissue, inputFile, refsColumns):

        df=pd.read_csv(os.path.join(self.dataDir, tissue, inputFile))
        methods=[c for c in df.columns if c != "gene"]

        outdir=os.path.join(self.dataDir, tissue, "bland_altman")
        os.makedirs(outdir,exist_ok=True)

        summary=[]
        for ref in refsColumns:
            refdir=os.path.join(outdir,ref)
            os.makedirs(refdir,exist_ok=True)

            methods_copy = methods.copy()
            methods_copy.remove(ref)
            for method in methods_copy:

                d=df[[ref,method]].dropna()
                d=d[(d[ref].notna())&(d[method].notna())]
                if d.empty: continue

                d["ref_log2"]=np.log2(d[ref]+1)
                d["test_log2"]=np.log2(d[method]+1)

                d["mean"]=(d["ref_log2"]+d["test_log2"])/2
                d["diff"]=d["test_log2"]-d["ref_log2"]

                # debug_file = os.path.join(refdir, f"debug_{ref}_{method}.csv")
                # d.to_csv(debug_file, index=False)                

                bias=d["diff"].mean()
                sd=d["diff"].std(ddof=1)
                loa=1.96*sd
                upper=bias+loa
                lower=bias-loa

                loa20=loa*.2
                upper20=bias+loa20
                lower20=bias-loa20

                total=len(d)
                genes20=d["diff"].between(lower20,upper20).sum()
                pct20=genes20/total*100

                summary.append({
                    "Reference":ref,
                    "Method":method,
                    "Bias":bias,
                    "SD":sd,
                    "Upper_LoA":upper,
                    "Lower_LoA":lower,
                    "20%_Upper":upper20,
                    "20%_Lower":lower20,
                    "Total_Genes":total,
                    "Genes_20%_LoA":genes20,
                    "%Genes_20%_LoA":pct20
                })

                plt.figure(figsize=(7,6))
                # plt.scatter(d["mean"],d["diff"],s=8,alpha=.5)
                plt.scatter(
                    d["mean"],
                    d["diff"],
                    s=7,          # smaller dots
                    alpha=1.0,    # slightly more transparent
                    linewidths=0, # no marker edge
                    marker="."
                )

                plt.axhline(bias,color="red",lw=2,label=f"Bias={bias:.3f}")
                plt.axhline(upper,color="blue",ls="--")
                plt.axhline(lower,color="blue",ls="--")
                plt.axhline(upper20,color="green",ls=":",label="20% LoA")
                plt.axhline(lower20,color="green",ls=":")
                plt.xlabel("Mean log₂(TPM)")
                plt.ylabel(f"{method}-{ref} log₂(TPM)")
                plt.title(f"{ref} vs {method}")

                plt.text(.02,.98,
                    f"Bias={bias:.3f}\nSD={sd:.3f}\nGenes={total}\n20% LoA={genes20} ({pct20:.1f}%)",
                    transform=plt.gca().transAxes,
                    va="top",
                    bbox=dict(facecolor="white",alpha=.85))

                plt.legend()
                plt.tight_layout()
                plt.savefig(os.path.join(refdir,f"{ref}_{method}.png"),dpi=300)
                plt.close()

        report=pd.DataFrame(summary)
        report.to_csv(os.path.join(outdir,"bland_altman_summary.csv"),index=False,float_format="%.9f")
        print(report.round(9))


if __name__=="__main__":
    ba=BlandAltman("/data/loads/liangh/rnaseqload/bland_altman")

    tissueRefs = {}
    # tissueRefs["heart"] = ["E-GEOD-65775", "E-GEOD-74747"]
    # tissueRefs["intestine"] = ["E-MTAB-2801"]
    # tissueRefs["lung"] = ["E-MTAB-4035","E-MTAB-8573"]
    # tissueRefs["spleen"] = ["E-MTAB-4035"]

    # miss match
    ba=BlandAltman("/data/loads/liangh/rnaseqload/bland_altman_mismatch")
    tissueRefs["heart"] = ["E-MTAB-4035"] 
    tissueRefs["intestine"] = ["E-GEOD-65775", "E-GEOD-74747"] 
    tissueRefs["lung"] = ["E-MTAB-2801"]
    tissueRefs["spleen"] = ["E-MTAB-4035","E-MTAB-8573"]


    for tissue, refs in tissueRefs.items():
        print(f"\n{tissue}: {refs}")
        ba.process(tissue, "merged_tpm.csv", refs)