import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FormatStrFormatter
import os
import sys
import re
import glob

def plot_meanfield_E(resfolder):
    # python plot_wf.py ./../res/0122/125200/meanfield_Energy.txt
    filename = os.path.join(resfolder, "meanfield_Energy.txt")

    # filenameをパスのコンポーネントに分割
    path_components = filename.split('/')
    index = -2

    # 新しいパスを作成
    savefolder = '/'.join(path_components[:index+1] + ['plots'])+"/"
    if not os.path.exists(savefolder):
        os.makedirs(savefolder)

    df = pd.read_csv(filename,sep="\s+")
    iter = df["iter"].values[1:]
    Kinetic = df["Kinetic"].values[1:]
    Skyrme = df["Skyrme"].values[1:]
    Coulomb = df["Coulomb"].values[1:]
    CM = df["CM"].values[1:]
    Total = df["Total"].values[1:]

    fig, ax = plt.subplots(figsize=(8,8))
    ax.plot(iter,Kinetic,label="Kinetic",ls = "-", lw=2)
    ax.plot(iter,Skyrme,label="Skyrme"  ,ls = "-", lw=2)
    ax.plot(iter,Coulomb,label="Coulomb",ls = "-", lw=2)
    ax.plot(iter,CM,label="CM"          ,ls = "-", lw=2)
    ax.plot(iter,Total,label="Total"    ,ls = "-", lw=2)
    ax.legend()
    ax.set_xlabel("Iteration Nubmer",fontsize=14)
    ax.set_ylabel("Energy [MeV]",fontsize=14)
    ymin,ymax = plt.ylim()
    #ax.set_ylim(ymin,0)
    ax.set_title("Mean Field Energy",fontsize=20)
    ax.grid(True, which="major", ls="--", color='0.65',lw=0.5)
    ax.grid(False, which="minor")
    plt.tight_layout()
    fig.savefig(savefolder+"meanfield_Energy.png",dpi=300)
    plt.close(fig)


    fig, ax = plt.subplots(figsize=(8,8))
    ax.plot(iter,df["Total"].diff().abs()[1:],label="$|E_{\\text{tot}}-E_{\\text{tot}_{\\text{prev}}}|$"
            ,  linestyle="-", linewidth=2)
    ax.legend(fontsize=12)
    ax.set_xlabel("Iteration Number",fontsize=14)
    ax.set_ylabel(r"$\Delta E_{tot}$ [MeV]",fontsize=14)
    ax.set_yscale("log")
    ax.set_title("Total Energy Change",fontsize=20)
    ax.grid(True, which="major", ls="--", color='0.65',lw=0.5)
    ax.grid(False, which="minor")
    plt.tight_layout()
    fig.savefig(savefolder+"meanfield_Energy_change.png",dpi=300)
    return

def main():
    filename=os.path.dirname(sys.argv[1])
    plot_meanfield_E(filename)
    return

if __name__ == '__main__':
    main()