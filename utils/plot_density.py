import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FormatStrFormatter
import os
import sys
import re
import glob

def plot_density(filename, z_val):
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

    # contour図
    r = df.iloc[:,0].values
    z = df.iloc[:,1].values
    # 各列のユニークな値の数を計算
    r = np.unique(df['r'])
    z = np.unique(df['z'])
    r_val = df["r"].min()
    z_val = min(z, key=lambda x:abs(x-z_val))

    Nr = len(r)
    Nz = len(z)
    # データを2次元配列にリシェイプ
    rho_n = df.iloc[:,2].values.reshape(Nr, Nz).T
    #rho_n = np.ma.masked_where(rho_n < 0.01, rho_n)
    rho_p = df.iloc[:,8].values.reshape(Nr, Nz).T
    tau_n = df.iloc[:,14].values.reshape(Nr, Nz).T
    tau_p = df.iloc[:,15].values.reshape(Nr, Nz).T
    jr_n = df.iloc[:,16].values.reshape(Nr, Nz).T
    jz_n = df.iloc[:,17].values.reshape(Nr, Nz).T
    jr_p = df.iloc[:,19].values.reshape(Nr, Nz).T
    jz_p = df.iloc[:,20].values.reshape(Nr, Nz).T

    Zlist = [rho_n,rho_p,tau_n,tau_p,jr_n,jz_n,jr_p,jz_p]
    names = ["rho_n","rho_p","tau_n","tau_p","jr_n","jz_n","jr_p","jz_p"]

    fig, axes = plt.subplots(4,2,figsize=(12,20))
    fig.suptitle("density contour map\n"+num+" Steps",fontsize=30,fontweight="bold")
    plot_range = 9
    for ax in axes.flat:
        ax.xaxis.set_major_locator(MaxNLocator(nbins=6))  # X軸の主要目盛りを5つに限定
        ax.yaxis.set_major_locator(MaxNLocator(nbins=6))  # Y軸の主要目盛りを5つに限定
        ax.set_aspect("equal", adjustable="box")  # アスペクト比を1:1にする
    for dens,ax,name in zip(Zlist,axes.flat,names):
        levels = np.linspace(np.min(dens), np.max(dens), 100)
        if(np.min(dens) >= np.max(dens)):
            min_value = -1e-10  # ここに適切な最小値を設定
            max_value = np.max(1e-10)  # densの最大値を使用
            levels = np.linspace(min_value, max_value, 200)

        cntr = ax.contourf(r,z,dens,levels=levels,cmap="plasma")
        ax.scatter(r_val,z_val,marker="x",color="green",s=200)
        # カラーバーを追加
        fig.colorbar(cntr, ax=ax,format=FormatStrFormatter('%.3f'))
        levels = np.linspace(np.min(dens), np.max(dens), 5)[1:-1]
        if(np.min(dens) >= np.max(dens)):
            min_value = -1e-10  # ここに適切な最小値を設定
            max_value = np.max(1e-10)  # densの最大値を使用
            levels = np.linspace(min_value, max_value, 5)[1:-1]
        ax.contour(r,z,dens,levels,colors="black",linestyles="dotted")
        ax.set_xlabel("r"+" [fm]")
        ax.set_ylabel("z"+" [fm]")
        ax.set_title(name)
        ax.set_xlim(0,plot_range)
        ax.set_ylim(z_val - plot_range/2,z_val + plot_range/2)

    fig.savefig(savefolder+"dens_contour.png")
    plt.close(fig)



    r_val = df["r"].min()
    df_p  = df[abs(df["z"] - z_val) < eps]
    df_p2 = df[abs(df["r"] - r_val) < eps]
    r = df_p.iloc[:,0].values
    z = df_p2.iloc[:,1].values
    dfs = [df_p,df_p2]
    xs = [r,z]
    # 逆転させる
    names = ["z","r"]
    vals = [z_val,r_val]
    xlabels = ["r","z"]
    idxs = [0,1]

    fig_r, axes_r = plt.subplots(4,2,figsize=(16,16),sharex=True)
    fig_t, axes_t = plt.subplots(1,2,figsize=(16,10),sharex=True)
    fig_j, axes_j = plt.subplots(2,2,figsize=(16,12),sharex=True)
    for idx,df,x,name,x_val,xlabel in zip(idxs,dfs,xs,names,vals,xlabels):
        rho_n     = df.iloc[:,2].values
        dr_rho_n  = df.iloc[:,3].values
        dz_rho_n  = df.iloc[:,4].values
        ddr_rho_n = df.iloc[:,5].values
        ddz_rho_n = df.iloc[:,6].values
        lap_rho_n = df.iloc[:,7].values
        rho_p     = df.iloc[:,8].values
        dr_rho_p  = df.iloc[:,9].values
        dz_rho_p  = df.iloc[:,10].values
        ddr_rho_p = df.iloc[:,11].values
        ddz_rho_p = df.iloc[:,12].values
        lap_rho_p = df.iloc[:,13].values
        tau_n     = df.iloc[:,14].values
        tau_p     = df.iloc[:,15].values
        jr_n      = df.iloc[:,16].values
        jz_n      = df.iloc[:,17].values
        divj_n    = df.iloc[:,18].values
        jr_p      = df.iloc[:,19].values
        jz_p      = df.iloc[:,20].values
        divj_p    = df.iloc[:,21].values
        tho_n     = df.iloc[:,22].values
        tho_p     = df.iloc[:,23].values

        axes_r[0][idx].plot(x,rho_n,label="rho_n")
        axes_r[0][idx].plot(x,rho_p,label="rho_p")
        axes_r[0][idx].legend()
        axes_r[0][idx].set_xlabel(xlabel+" [fm]")
        axes_r[0][idx].set_ylabel("rho [fm^-3]")
        axes_r[0][idx].set_title("rho at "+name+" "+str(x_val)+"fm")
        axes_r[1][idx].plot(x,dr_rho_n,label="dr_rho_n")
        axes_r[1][idx].plot(x,dr_rho_p,label="dr_rho_p")
        axes_r[1][idx].plot(x,dz_rho_n,label="dz_rho_n")
        axes_r[1][idx].plot(x,dz_rho_p,label="dz_rho_p")
        axes_r[1][idx].legend()
        axes_r[1][idx].set_xlabel(xlabel+" [fm]")
        axes_r[1][idx].set_ylabel("d_rho [fm^-4]")
        axes_r[1][idx].set_title("d_rho at "+name+" "+str(x_val)+"fm")
        axes_r[2][idx].plot(x,ddr_rho_n,label="ddr_rho_n")
        axes_r[2][idx].plot(x,ddr_rho_p,label="ddr_rho_p")
        axes_r[2][idx].plot(x,ddz_rho_n,label="ddz_rho_n")
        axes_r[2][idx].plot(x,ddz_rho_p,label="ddz_rho_p")
        axes_r[2][idx].legend()
        axes_r[2][idx].set_xlabel(xlabel+" [fm]")
        axes_r[2][idx].set_ylabel("dd_rho [fm^-5]")
        axes_r[2][idx].set_title("dd_rho at "+name+" "+str(x_val)+"fm")
        axes_r[3][idx].plot(x,lap_rho_n,label="lap_rho_n")
        axes_r[3][idx].plot(x,lap_rho_p,label="lap_rho_p")
        axes_r[3][idx].set_xlabel(xlabel+" [fm]")
        axes_r[3][idx].set_ylabel("lap_rho [fm^-5]")
        axes_r[3][idx].set_title("lap_rho at "+name+" "+str(x_val)+"fm")
        axes_r[3][idx].legend()

        axes_t[idx].plot(x,tau_n,label="tau_n")
        axes_t[idx].plot(x,tau_p,label="tau_p")
        axes_t[idx].plot(x,tho_n,label="tho_n")
        axes_t[idx].set_xlabel(xlabel+" [fm]")
        axes_t[idx].set_ylabel("tau [fm^-4]")
        axes_t[idx].legend()
        axes_t[idx].set_title("tau at "+name+" "+str(x_val)+"fm")

        axes_j[0][idx].plot(x,jr_n,label="jr_n")
        axes_j[0][idx].plot(x,jz_n,label="jz_n")
        axes_j[0][idx].plot(x,jr_p,label="jr_p")
        axes_j[0][idx].plot(x,jz_p,label="jz_p")
        axes_j[0][idx].set_xlabel(xlabel+" [fm]")
        axes_j[0][idx].set_ylabel("j [fm^-2]")
        axes_j[0][idx].legend()
        axes_j[0][idx].set_title("Spin-Orbit dens at "+name+" "+str(x_val)+"fm")

        axes_j[1][idx].plot(x,divj_n,label="divj_n")
        axes_j[1][idx].plot(x,divj_p,label="divj_p")
        axes_j[1][idx].set_xlabel(xlabel+" [fm]")
        axes_j[1][idx].set_ylabel("j [fm^-2]")
        axes_j[1][idx].legend()
        axes_j[1][idx].set_title("divJ at "+name+" "+str(x_val)+"fm")

    fig_r.savefig(savefolder+"rho.png")
    fig_j.savefig(savefolder+"Spin-Orbit.png")
    fig_t.savefig(savefolder+"tau.png")
    fig_r.clf()
    fig_j.clf()
    fig_t.clf()
    plt.close(fig_r)
    plt.close(fig_j)
    plt.close(fig_t)
    return


def plot_density_all(resfolder, z_val):
    folder_path = os.path.join(resfolder, 'densities/')
    txt_files = glob.glob(os.path.join(folder_path, 'density_and_derivative[0-9]*.txt'))
    for file in txt_files:
        plot_density(file, z_val)
    return

def main():
    resfolder = os.path.dirname(sys.argv[1])
    z_val = float(sys.argv[2])
    plot_density_all(resfolder, z_val)

if __name__ == "__main__":
    main()







