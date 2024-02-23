import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import re
import glob

def plot_potential(filename, z_val):
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
    z_val = min(df["z"].values, key=lambda x:abs(x-z_val))
    df_p = df[abs(df["z"] - z_val) < eps]
    #"r z Coulomb U_n U_p B_n B_p Wr_n Wz_n, Wr_p Wz_p"
    x = df_p.iloc[:,0].values
    U_c = df_p.iloc[:,2].values
    U_n = df_p.iloc[:,3].values
    U_p = df_p.iloc[:,4].values
    B_n = df_p.iloc[:,5].values
    B_p = df_p.iloc[:,6].values
    Wr_n = df_p.iloc[:,7].values
    Wz_n = df_p.iloc[:,8].values
    Wr_p = df_p.iloc[:,9].values
    Wz_p = df_p.iloc[:,10].values

    df_p2 = df[abs(df["r"] - df["r"].min()) < eps]
    z = df_p2.iloc[:,1].values
    zU_c = df_p2.iloc[:,2].values
    zU_n = df_p2.iloc[:,3].values
    zU_p = df_p2.iloc[:,4].values
    zB_n = df_p2.iloc[:,5].values
    zB_p = df_p2.iloc[:,6].values
    zWr_n = df_p2.iloc[:,7].values
    zWz_n = df_p2.iloc[:,8].values
    zWr_p = df_p2.iloc[:,9].values
    zWz_p = df_p2.iloc[:,10].values

    fig,axes = plt.subplots(2,1,figsize=(11,11),sharex=True)
    axes[0].plot(x,U_n + B_n + Wr_n,label="Total potential_n")
    axes[0].plot(x,U_p + B_p + Wr_p,label="Total potential_p")
    axes[0].legend()
    axes[0].set_xlabel("r [fm]")
    axes[0].set_ylabel("Total Pot[MeV]")
    axes[0].set_title("Total potential at z= "+str(z_val)+"fm "+num)

    axes[1].plot(x,zU_n + zB_n + zWr_n,label="Total potential_n")
    axes[1].plot(x,zU_p + zB_p + zWr_p,label="Total potential_p")
    axes[1].legend()
    axes[1].set_xlabel("z [fm]")
    axes[1].set_title("Total potential at r= "+str(df["r"].min())+"fm "+num)

    fig.savefig(savefolder+"Total_potential.png")
    plt.clf()
    plt.close()




    fig, axes = plt.subplots(2,2,figsize=(11,11),sharex=True)
    axes[0][0].plot(x,U_p,label="U_p")
    axes[0][0].plot(x,U_n,label="U_n")
    axes[0][0].legend()
    #axes[0][0].set_xlabel("r [fm]")
    axes[0][0].set_ylabel("U [MeV]")
    axes[0][0].set_title("U at z= "+str(z_val)+"fm"+num)


    axes[0][1].plot(x,B_n,label="B_n")
    axes[0][1].plot(x,B_p,label="B_p")
    #axes[0][1].set_xlabel("r [fm]")
    axes[0][1].set_ylabel("B [MeV.fm^2]")
    axes[0][1].set_title("B at z= "+str(z_val)+"fm "+num)
    axes[0][1].legend()
    
    axes[1][0].plot(x,Wr_n,label="Wr_n")
    axes[1][0].plot(x,Wz_n,label="Wz_n")
    axes[1][0].plot(x,Wr_p,label="Wr_p")
    axes[1][0].plot(x,Wz_p,label="Wz_p")
    axes[1][0].set_xlabel("r [fm]")
    axes[1][0].set_ylabel("W [MeV.fm]")
    axes[1][0].set_title("W at z= "+str(z_val)+"fm "+num)
    axes[1][0].legend() 
    
    df_p2 = df[abs(df["r"] - df["r"].min()) < eps]
    z = df_p2.iloc[:,1].values
    Wr_n2 = df_p2.iloc[:,7].values
    Wz_n2 = df_p2.iloc[:,8].values
    Wr_p2 = df_p2.iloc[:,9].values
    Wz_p2 = df_p2.iloc[:,10].values
    
    axes[1][1].plot(x,Wr_n2,label="Wr_n")
    axes[1][1].plot(x,Wz_n2,label="Wz_n")
    axes[1][1].plot(x,Wr_p2,label="Wr_p")
    axes[1][1].plot(x,Wz_p2,label="Wz_p")
    axes[1][1].set_xlabel("z [fm]")
    axes[1][1].set_ylabel("W [MeV.fm]")
    axes[1][1].set_title("W at r= "+str(df["r"].min())+"fm "+num)
    axes[1][1].legend()
    fig.savefig(savefolder+"potential.png")
    plt.clf()
    plt.close()
    
    plt.figure()
    plt.plot(x,U_c,label="Direct_term")
    plt.legend()
    plt.title("Direct term")
    plt.xlabel("r [fm]")
    plt.ylabel("U_d [MeV]")
    plt.savefig(savefolder+"Direct_term.png")
    plt.clf()
    plt.close()
    return 

def plot_potential_all(resfolder, z_val):
    folder_path = os.path.join(resfolder,'potentials/')
    txt_files = glob.glob(os.path.join(folder_path, 'pot[0-9]*.txt'))
    for file in txt_files:
        plot_potential(file, z_val)
    return

def main():
    resfolder = os.path.dirname(sys.argv[1])
    z_val = float(sys.argv[2])
    plot_potential_all(resfolder, z_val)
    return

if __name__ == "__main__":
    main()
