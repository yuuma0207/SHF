import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import os
import sys
import re
import glob

def plot_d_wf(filename, z_val):
    # python plot_wf.py ./../res/0117/132937/densities/wavefunc0000.txt 10
    pattern = r'(\d{1,})\.txt$'
    num = re.search(pattern,filename).group(1)
    # filenameをパスのコンポーネントに分割
    path_components = filename.split('/')

    # 'orbital_wavefunction' が存在するインデックスを探す
    index = -3

    # 新しいパスを作成
    savefolder = '/'.join(path_components[:index+1] + ['plots'] + [path_components[index+1]])+"/"
    if not os.path.exists(savefolder):
        os.makedirs(savefolder)
    for file in os.listdir(savefolder):
        if(file.startswith(num)):
            return None
    savefolder = savefolder+ num

    df = pd.read_csv(filename,sep="\s+")
    eps = 0.0001
    r_val = df["r"].min()
    z_val = min(df["z"].values, key=lambda x:abs(x-z_val))
    df_p  = df[abs(df["z"] - z_val) < eps]
    df_p2 = df[abs(df["r"] - r_val) < eps]
    r = df_p.iloc[:,0].values
    z = df_p2.iloc[:,1].values
    dfs = [df_p,df_p2]
    xs = [r,z]
    # 逆転させる
    xlabels = ["r","z"]
    names = ["z","r"]
    vals = [z_val,r_val]
    idxs = [0,1]

    #"r z dr_wfn1u dz_wfn1u dr_wfn1d dz_wfn1d &dr_wfn2u dz_wfn2u dr_wfn2d dz_wfn2d &dr_wfp1u dz_wfp1u dr_wfp1d dz_wfp1d &dr_wfp2u dz_wfp2u dr_wfp2d dz_wfp2d "
    fig, axes = plt.subplots(8,2,figsize=(30,30),sharex=True)
    fig.suptitle("one times derivative of wave function\n"+num+" Steps",fontsize=30,fontweight="bold")
    for ax in axes.flat:
        #ax.xaxis.set_major_locator(MaxNLocator(nbins=5))  # X軸の主要目盛りを5つに限定
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))  # Y軸の主要目盛りを5つに限定
    for idx,df,x,name,x_val,xlabel in zip(idxs,dfs,xs,names,vals,xlabels):
        dr_wfn1u = df.iloc[:,2].values
        dz_wfn1u = df.iloc[:,3].values
        dr_wfn1d = df.iloc[:,4].values
        dz_wfn1d = df.iloc[:,5].values
        dr_wfn2u = df.iloc[:,6].values
        dz_wfn2u = df.iloc[:,7].values
        dr_wfn2d = df.iloc[:,8].values
        dz_wfn2d = df.iloc[:,9].values
        dr_wfp1u = df.iloc[:,10].values
        dz_wfp1u = df.iloc[:,11].values
        dr_wfp1d = df.iloc[:,12].values
        dz_wfp1d = df.iloc[:,13].values
        dr_wfp2u = df.iloc[:,14].values
        dz_wfp2u = df.iloc[:,15].values
        dr_wfp2d = df.iloc[:,16].values
        dz_wfp2d = df.iloc[:,17].values

        axes[0][idx].plot(x,dr_wfn1u,label="dr_wfn1u")
        axes[0][idx].plot(x,dz_wfn1u,label="dz_wfn1u")
        axes[0][idx].legend()
        axes[0][idx].set_xlabel(xlabel+" [fm]")
        axes[0][idx].set_ylabel("d_wfn1 [fm^-5/2]")
        axes[0][idx].set_title("d_wfn1u at "+name+" "+str(x_val)+"fm")

        axes[1][idx].plot(x,dr_wfn1d,label="dr_wfn1d")
        axes[1][idx].plot(x,dz_wfn1d,label="dz_wfn1d")
        axes[1][idx].legend()
        axes[1][idx].set_xlabel(xlabel+" [fm]")
        axes[1][idx].set_ylabel("d_wfn1 [fm^-5/2]")
        axes[1][idx].set_title("d_wfn1d at "+name+" "+str(x_val)+"fm")

        axes[2][idx].plot(x,dr_wfn2u,label="dr_wfn2u")
        axes[2][idx].plot(x,dz_wfn2u,label="dz_wfn2u")
        axes[2][idx].legend()
        axes[2][idx].set_xlabel(xlabel+" [fm]")
        axes[2][idx].set_ylabel("d_wfn2 [fm^-5/2]")
        axes[2][idx].set_title("d_wfn2u at "+name+" "+str(x_val)+"fm")

        axes[3][idx].plot(x,dr_wfn2d,label="dr_wfn2d")
        axes[3][idx].plot(x,dz_wfn2d,label="dz_wfn2d")
        axes[3][idx].legend()
        axes[3][idx].set_xlabel(xlabel+" [fm]")
        axes[3][idx].set_ylabel("d_wfn2 [fm^-5/2]")
        axes[3][idx].set_title("d_wfn2d at "+name+" "+str(x_val)+"fm")

        axes[4][idx].plot(x,dr_wfp1u,label="dr_wfp1u")
        axes[4][idx].plot(x,dz_wfp1u,label="dz_wfp1u")
        axes[4][idx].legend()
        axes[4][idx].set_xlabel(xlabel+" [fm]")
        axes[4][idx].set_ylabel("d_wfp1 [fm^-5/2]")
        axes[4][idx].set_title("d_wfp1u at "+name+" "+str(x_val)+"fm")

        axes[5][idx].plot(x,dr_wfp1d,label="dr_wfp1d")
        axes[5][idx].plot(x,dz_wfp1d,label="dz_wfp1d")
        axes[5][idx].legend()
        axes[5][idx].set_xlabel(xlabel+" [fm]")
        axes[5][idx].set_ylabel("d_wfp1 [fm^-5/2]")
        axes[5][idx].set_title("d_wfp1d at "+name+" "+str(x_val)+"fm")

        axes[6][idx].plot(x,dr_wfp2u,label="dr_wfp2u")
        axes[6][idx].plot(x,dz_wfp2u,label="dz_wfp2u")
        axes[6][idx].legend()
        axes[6][idx].set_xlabel(xlabel+" [fm]")
        axes[6][idx].set_ylabel("d_wfp2 [fm^-5/2]")
        axes[6][idx].set_title("d_wfp2u at "+name+" "+str(x_val)+"fm")

        axes[7][idx].plot(x,dr_wfp2d,label="dr_wfp2d")
        axes[7][idx].plot(x,dz_wfp2d,label="dz_wfp2d")
        axes[7][idx].legend()
        axes[7][idx].set_xlabel(xlabel+" [fm]")
        axes[7][idx].set_ylabel("d_wfp2 [fm^-5/2]")
        axes[7][idx].set_title("d_wfp2d at "+name+" "+str(x_val)+"fm")


    fig.savefig(savefolder+"wf_d.png")
    plt.close(fig)

def plot_d_wf_all(resfolder, z_val):
    folder_path = os.path.join(resfolder, 'orbital_wavefunction/')
    txt_files = glob.glob(os.path.join(folder_path, 'wavefunc_d[0-9]*.txt'))
    for file in txt_files:
        plot_d_wf(file, z_val)
    return

def main():
    resfolder = os.path.dirname(sys.argv[1])
    z_val = float(sys.argv[2])
    plot_d_wf_all(resfolder, z_val)

if __name__ == "__main__":
    main()

