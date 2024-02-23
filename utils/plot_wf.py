import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FormatStrFormatter
import os
import sys
import re
import glob

def plot_wf(filename, z_val):
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
    
    eps = 0.0001
    df = pd.read_csv(filename,sep="\s+")
    r_val = df["r"].min()
    z_val = min(df["z"].values, key=lambda x:abs(x-z_val))

    # conter図
    r = df.iloc[:,0].values
    z = df.iloc[:,1].values
    # 各列のユニークな値の数を計算
    r = np.unique(df['r'])
    z = np.unique(df['z'])

    Nr = len(r)
    Nz = len(z)

    # データを2次元配列にリシェイプ
    wfn1u = df.iloc[:,2].values.reshape(Nr, Nz).T
    wfn1d = df.iloc[:,3].values.reshape(Nr, Nz).T
    wfn2u = df.iloc[:,4].values.reshape(Nr, Nz).T
    wfn2d = df.iloc[:,5].values.reshape(Nr, Nz).T
    wfp1u = df.iloc[:,6].values.reshape(Nr, Nz).T
    wfp1d = df.iloc[:,7].values.reshape(Nr, Nz).T
    wfp2u = df.iloc[:,8].values.reshape(Nr, Nz).T
    wfp2d = df.iloc[:,9].values.reshape(Nr, Nz).T

    Zlist = [wfn1u,wfn1d,wfn2u,wfn2d,wfp1u,wfp1d,wfp2u,wfp2d]
    names = ["wfn1u","wfn1d","wfn2u","wfn2d","wfp1u","wfp1d","wfp2u","wfp2d"]

    fig, axes = plt.subplots(4,2,figsize=(12,20))
    fig.suptitle("wave function contour map\n"+num+" Steps",fontsize=30,fontweight="bold")
    for ax in axes.flat:
        ax.xaxis.set_major_locator(MaxNLocator(nbins=6))  # X軸の主要目盛りを5つに限定
        ax.yaxis.set_major_locator(MaxNLocator(nbins=6))  # Y軸の主要目盛りを5つに限定
        ax.set_aspect("equal", adjustable="box")  # アスペクト比を1:1にする

    for wf,ax,name in zip(Zlist,axes.flat,names):
        levels = np.linspace(np.min(wf), np.max(wf), 200)
        cntr = ax.contourf(r,z,wf,levels=levels,cmap="plasma")
        # カラーバーを追加
        cbar = fig.colorbar(cntr, ax=ax,format=FormatStrFormatter('%.3f'))
        levels = np.linspace(np.min(wf), np.max(wf), 5)[1:-1]
        ax.contour(r,z,wf,levels,colors="black",linestyles="dotted")
        ax.scatter(r_val,z_val,marker="x",color="green",s=200)
        ax.set_xlabel("r"+" [fm]")
        ax.set_ylabel("z"+" [fm]")
        ax.set_title(name)
        ax.set_xlim(0,4)
        ax.set_ylim(8,12)

    fig.savefig(savefolder+"wf_contour.png")
    plt.close(fig)

    df_p  = df[abs(df["z"] - z_val) < eps]
    df_p2 = df[abs(df["r"] - r_val) < eps]
    r = df_p.iloc[:,0].values
    z = df_p2.iloc[:,1].values
    dfs = [df_p,df_p2]
    xs = [r,z]
    # 逆転させる
    names = ["z","r"]
    xlabels = ["r","z"]
    vals = [z_val,r_val]
    idxs = [0,1]

    # write(fi,*) "r z wfn1u wfn1d wfn2u wfn2d wfp1u wfp1d wfp2u wfp2d"
    fig, axes = plt.subplots(8,2,figsize=(30,30),sharex=True)
    fig.suptitle("wave function\n"+num+" Steps",fontsize=30,fontweight="bold")
    for ax in axes.flat:
        #ax.xaxis.set_major_locator(MaxNLocator(nbins=5))  # X軸の主要目盛りを5つに限定
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))  # Y軸の主要目盛りを5つに限定
    for idx,df,x,name,x_val,xlabel in zip(idxs,dfs,xs,names,vals,xlabels):
        wfn1u = df.iloc[:,2].values
        wfn1d = df.iloc[:,3].values
        wfn2u = df.iloc[:,4].values
        wfn2d = df.iloc[:,5].values
        wfp1u = df.iloc[:,6].values
        wfp1d = df.iloc[:,7].values
        wfp2u = df.iloc[:,8].values
        wfp2d = df.iloc[:,9].values

        axes[0][idx].plot(x,wfn1u,label="wfn1u")
        axes[0][idx].legend()
        axes[0][idx].set_xlabel(xlabel+" [fm]")
        axes[0][idx].set_ylabel("wfn1 [fm^-3/2]")
        axes[0][idx].set_title("wfn1u at "+name+" "+str(x_val)+"fm")

        axes[1][idx].plot(x,wfn1d,label="wfn1d")
        axes[1][idx].legend()
        axes[1][idx].set_xlabel(xlabel+" [fm]")
        axes[1][idx].set_ylabel("wfn1 [fm^-3/2]")
        axes[1][idx].set_title("wfn1d at "+name+" "+str(x_val)+"fm")

        axes[2][idx].plot(x,wfn2u,label="wfn2u")
        axes[2][idx].legend()
        axes[2][idx].set_xlabel(xlabel+" [fm]")
        axes[2][idx].set_ylabel("wfn2 [fm^-3/2]")
        axes[2][idx].set_title("wfn2u at "+name+" "+str(x_val)+"fm")

        axes[3][idx].plot(x,wfn2d,label="wfn2d")
        axes[3][idx].legend()
        axes[3][idx].set_xlabel(xlabel+" [fm]")
        axes[3][idx].set_ylabel("wfn2 [fm^-3/2]")
        axes[3][idx].set_title("wfn2d at "+name+" "+str(x_val)+"fm")

        axes[4][idx].plot(x,wfp1u,label="wfp1u")
        axes[4][idx].legend()
        axes[4][idx].set_xlabel(xlabel+" [fm]")
        axes[4][idx].set_ylabel("wfp1 [fm^-3/2]")
        axes[4][idx].set_title("wfp1u at "+name+" "+str(x_val)+"fm")

        axes[5][idx].plot(x,wfp1d,label="wfp1d")
        axes[5][idx].legend()
        axes[5][idx].set_xlabel(xlabel+" [fm]")
        axes[5][idx].set_ylabel("wfp1 [fm^-3/2]")
        axes[5][idx].set_title("wfp1d at "+name+" "+str(x_val)+"fm")

        axes[6][idx].plot(x,wfp2u,label="wfp2u")
        axes[6][idx].legend()
        axes[6][idx].set_xlabel(xlabel+" [fm]")
        axes[6][idx].set_ylabel("wfp2 [fm^-3/2]")
        axes[6][idx].set_title("wfp2u at "+name+" "+str(x_val)+"fm")

        axes[7][idx].plot(x,wfp2d,label="wfp2d")
        axes[7][idx].legend()
        axes[7][idx].set_xlabel(xlabel+" [fm]")
        axes[7][idx].set_ylabel("wfp2 [fm^-3/2]")
        axes[7][idx].set_title("wfp2d at "+name+" "+str(x_val)+"fm")
    fig.savefig(savefolder+"wf.png")
    plt.close(fig)
    return

def plot_wf_all(resfolder, z_val):
    folder_path = os.path.join(resfolder, 'orbital_wavefunction/')
    txt_files = glob.glob(os.path.join(folder_path + 'wavefunc[0-9]*.txt'))
    for file in txt_files:
        plot_wf(file, z_val)
    return

def main():
    resfolder = os.path.dirname(sys.argv[1])
    z_val = float(sys.argv[2])
    plot_wf_all(resfolder, z_val)

if __name__ == "__main__":
    main()