import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import re

#python plot_pot_term.py ./../res/0117/132937/densities/density_and_derivative0000.txt 10
filename=sys.argv[1]
z_val = float(sys.argv[2])
pattern = r'(\w\d{4})\.txt$'
num = re.search(pattern,filename).group(1)
# filenameをパスのコンポーネントに分割
path_components = filename.split('/')

# 'orbital_wavefunction' が存在するインデックスを探す
index = -3

# 新しいパスを作成
savefolder = '/'.join(path_components[:index+1] + ['plots'] + [path_components[index+1]])+"/"
if not os.path.exists(savefolder):
    os.makedirs(savefolder)
savefolder = savefolder+ num

eps = 0.0001
df = pd.read_csv(filename,sep="\s+")
df_p = df[abs(df["z"] - z_val) < eps]
#"r z b0 b0p b1 b1p b2 b2p b3 b3p b3p2 b4 b4p dir ex"
x = df_p.iloc[:,0].values
b0 = df_p.iloc[:,2].values
b0p = df_p.iloc[:,3].values
b1 = df_p.iloc[:,4].values
b1p = df_p.iloc[:,5].values
b2 = df_p.iloc[:,6].values
b2p = df_p.iloc[:,7].values
b3 = df_p.iloc[:,8].values
b3p = df_p.iloc[:,9].values
b3p2 = df_p.iloc[:,10].values
b4 = df_p.iloc[:,11].values
b4p = df_p.iloc[:,12].values
dir = df_p.iloc[:,13].values
ex = df_p.iloc[:,14].values

fig, axes = plt.subplots(3,2,figsize=(11,11),sharex=True)

axes[0][0].plot(x,b0,label="b0")
axes[0][0].plot(x,b0p,label="b0p")
#axes[0][0].set_xlabel("r [fm]")
axes[0][0].set_ylabel("b0_term [MeV]")
axes[0][0].set_title("b0_term at z= "+str(z_val)+"fm")
axes[0][0].legend()

axes[0][1].plot(x,b1,label="b1")
axes[0][1].plot(x,b1p,label="b1p")
#axes[0][1].set_xlabel("r [fm]")
axes[0][1].set_ylabel("b1_term [MeV]")
axes[0][1].set_title("b1_term at z= "+str(z_val)+"fm")
axes[0][1].legend()

axes[1][0].plot(x,b2,label="b2")
axes[1][0].plot(x,b2p,label="b2p")
#axes[1][0].set_xlabel("r [fm]")
axes[1][0].set_ylabel("b2_term [MeV]")
axes[1][0].set_title("b2_term at z= "+str(z_val)+"fm")
axes[1][0].legend()

axes[1][1].plot(x,b3,label="b3")
axes[1][1].plot(x,b3p,label="b3p")
axes[1][1].plot(x,b3p2,label="b3p2")
#axes[1][1].set_xlabel("r [fm]")
axes[1][1].set_ylabel("b3_term [MeV]")
axes[1][1].set_title("b3_term at z= "+str(z_val)+"fm")
axes[1][1].legend()

axes[2][0].plot(x,b4,label="b4")
axes[2][0].plot(x,b4p,label="b4p")
axes[2][0].set_xlabel("r [fm]")
axes[2][0].set_ylabel("b4_term [MeV]")
axes[2][0].set_title("b4_term at z= "+str(z_val)+"fm")
axes[2][0].legend()

axes[2][1].plot(x,dir,label="Direct")
axes[2][1].plot(x,ex,label="Ex")
axes[2][1].set_xlabel("r [fm]")
axes[2][1].set_ylabel("Coulomb_term [MeV]")
axes[2][1].set_title("Coulomb_term at z= "+str(z_val)+"fm")
axes[2][1].legend()

plt.savefig(savefolder+"Uterm.png")

plt.clf()
plt.close()


fig, axes = plt.subplots(3,2,figsize=(11,11),sharex=True)

axes[0][0].plot(x,b0-b0p,label="b0-b0p")
#axes[0][0].set_xlabel("r [fm]")
axes[0][0].set_ylabel("b0_term [MeV]")
axes[0][0].set_title("b0_term at z= "+str(z_val)+"fm")
axes[0][0].legend()

axes[0][1].plot(x,b1-b1p,label="b1-b1p")
#axes[0][1].set_xlabel("r [fm]")
axes[0][1].set_ylabel("b1_term [MeV]")
axes[0][1].set_title("b1_term at z= "+str(z_val)+"fm")
axes[0][1].legend()

axes[1][0].plot(x,-b2+b2p,label="-b2+b2p")
#axes[1][0].set_xlabel("r [fm]")
axes[1][0].set_ylabel("b2_term [MeV]")
axes[1][0].set_title("b2_term at z= "+str(z_val)+"fm")
axes[1][0].legend()

axes[1][1].plot(x,b3-b3p-b3p2,label="b3-b3p-b3p2")
#axes[1][1].set_xlabel("r [fm]")
axes[1][1].set_ylabel("b3_term [MeV]")
axes[1][1].set_title("b3_term at z= "+str(z_val)+"fm")
axes[1][1].legend()

axes[2][0].plot(x,-b4-b4p,label="-b4-b4p")
axes[2][0].set_xlabel("r [fm]")
axes[2][0].set_ylabel("b4_term [MeV]")
axes[2][0].set_title("b4_term at z= "+str(z_val)+"fm")
axes[2][0].legend()

axes[2][1].plot(x,dir-ex,label="Direct - Ex")
axes[2][1].set_xlabel("r [fm]")
axes[2][1].set_ylabel("Coulomb_term [MeV]")
axes[2][1].set_title("Coulomb_term at z= "+str(z_val)+"fm")
axes[2][1].legend()

plt.savefig(savefolder+"Uterm_sum.png")

plt.clf()
plt.close()

