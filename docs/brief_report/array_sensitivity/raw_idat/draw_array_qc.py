"""Plot IDAT-derived control and intensity metrics by Sentrix row."""
import argparse
from pathlib import Path
from textwrap import fill

import matplotlib
matplotlib.use("Agg")

from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs",type=Path)
    parser.add_argument("pdf",type=Path)
    a=parser.parse_args()
    qc=pd.read_csv(a.inputs/"array_qc.tsv",sep="\t").sort_values("sample")
    controls=pd.read_csv(a.inputs/"control_metrics_long.tsv",sep="\t")
    labels=[f"{x}\n{y[:3]}" for x,y in zip(qc.label,qc.position,strict=True)]
    cols=list(qc.columns[4:])
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"pdf.fonttype":42})
    color=["#21688B"]*3+["#A24B1D"]*3
    group_labels=[Line2D([],[],marker="o",linestyle="",color=color[0],label="Control"),
                  Line2D([],[],marker="o",linestyle="",color=color[3],label="DNMT3L overexpression"),
                  Line2D([],[],linestyle="--",color="#627383",label="ewastools cutoff")]

    settings=[]
    settings.append(("Raw fluorescence and detection",[(n,qc[n].to_numpy(),np.nan) for n in cols]))
    groups=list(controls.groupby("metric",sort=False))
    for i in range(0,len(groups),12):
        heading="Illumina control metrics" if i==0 else "Illumina control metrics continued"
        settings.append((heading,[(n,part.sort_values("sample").value.to_numpy(),part.threshold.iloc[0])
                                 for n,part in groups[i:i+12]]))
    with PdfPages(a.pdf, metadata={"Title":"GSE157078 IDAT quality control by chip row","CreationDate":None}) as pdf:
        for i,(title,panels) in enumerate(settings,1):
            fig, grid=plt.subplots(4,3,figsize=(11.1,8.4))
            fig.subplots_adjust(left=.07,right=.982,bottom=.115,top=.86,hspace=.78,wspace=.34)
            fig.text(.07,.95,title,fontsize=20,weight="bold",color="#173046")
            fig.text(.07,.908,"GSE157078 · Sentrix chip 203774260082 · Six arrays",fontsize=10,color="#445460")
            for ax,(name,values,threshold) in zip(grid.flat,panels):
                ax.axvspan(2.5,2.51,color="#C9D0D4",alpha=.8)
                ax.plot(range(6),values,color="#98A4AC",linewidth=1.1,zorder=2)
                ax.scatter(range(6),values,c=color,s=24,zorder=3)
                if np.isfinite(threshold):
                    low, high=ax.get_ylim()
                    if low<=threshold<=high or np.any(values<threshold):
                        ax.axhline(threshold,color="#627383",linestyle="--",linewidth=.9)
                    ax.text(.98,.94,f"Cutoff {threshold:g}",transform=ax.transAxes,
                            ha="right",va="top",fontsize=7,color="#445460",
                            bbox={"facecolor":"white","alpha":.8,"edgecolor":"none","pad":1})
                ax.set_title(fill(name,35),loc="left",fontsize=8.3,weight="bold",pad=5)
                ax.set_xticks(range(6),labels,fontsize=7)
                ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
                ax.tick_params(axis="y",labelsize=7)
                ax.tick_params(axis="both",length=2,color="#8898A3")
                ax.spines[["top","right"]].set_visible(False)
                ax.spines[["bottom","left"]].set_color("#A0AEB8")
                ax.set_xlim(-.25,5.25)
            for ax in list(grid.flat)[len(panels):]:
                ax.set_axis_off()
            fig.legend(handles=group_labels[:2] if i==1 else group_labels,
                       loc="lower center",ncol=3,frameon=False,bbox_to_anchor=(.5,.024))
            fig.text(.978,.043,str(i),ha="right",fontsize=8,color="#445460")
            pdf.savefig(fig)
            plt.close(fig)


if __name__=="__main__":
    main()
