import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import re
from PIL import Image

def make_figure_Pot_rdir(z_val):
    df_p = df[abs(df["z"] - z_val) < eps]

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

    figure,axes = plt.subplots(1,2,figsize=(11,11),sharey=True)
    axes[0].plot(x,U_n + B_n + Wr_n,label="Total potential_n")
    axes[0].plot(x,U_p + B_p + Wr_p,label="Total potential_p")
    axes[0].legend()
    axes[0].set_xlabel("r [fm]")
    axes[0].set_ylabel("Total Pot[MeV]")
    axes[0].set_title("Total potential at z= "+str(z_val)+"fm")

    axes[1].plot(x,zU_n + zB_n + zWr_n,label="Total potential_n")
    axes[1].plot(x,zU_p + zB_p + zWr_p,label="Total potential_p")
    axes[1].legend()
    axes[1].set_xlabel("z [fm]")
    axes[1].set_title("Total potential at r= "+str(df["r"].min())+"fm")
    return figure

def plt_to_img(fig):
    # Figure をレンダリングする。
    fig.canvas.draw()
    # 画像をバイト列で取得する。
    data = fig.canvas.tostring_rgb()

    # 画像の大きさを取得する。
    w, h = fig.canvas.get_width_height()
    c = len(data) // (w * h)

    # PIL.Image に変換する
    img = Image.frombytes("RGB", (w, h), data, "raw")
    return img

def make_gifanime(num_gif,duration):
    num_gif = num_gif
    duration = duration
    pictures = []
    for i in range(num_gif):
        skip = len(sort_z)//num_gif
        gif_z = sort_z[i + skip]
        fig = make_figure_Pot_rdir(gif_z)
        im = plt_to_img(fig)
        pictures.append(im)
        plt.close(fig)  # matplotlibの図を閉じる

    # gifアニメを出力する
    pictures[0].save('pot_anime.gif', save_all=True, append_images=pictures[1:], optimize=True, duration=duration, loop=0)
    return()

# python plot_pot.py ./../res/0117/132937/potentials/potential0000.txt 10
filename=sys.argv[1]
pattern = r'(\d{1,})\.txt$'
num = re.search(pattern,filename).group(1)
savefolder = os.path.dirname(filename) + "/plots/" + num

eps = 0.0001
df = pd.read_csv(filename,sep="\s+")
sort_z = np.sort(df["z"].unique())
pic_z = float(sys.argv[2])

make_gifanime(10,100)
figure = make_figure_Pot(pic_z)
figure.savefig(savefolder+"Total_potential.png")
plt.clf()
plt.close()

def make_figure_Pot_zdir(r_val):
    return figure


fig, axes = plt.subplots(2,2,figsize=(11,11),sharex=True)
axes[0][0].plot(x,U_p,label="U_p")
axes[0][0].plot(x,U_n,label="U_n")
axes[0][0].legend()
#axes[0][0].set_xlabel("r [fm]")
axes[0][0].set_ylabel("U [MeV]")
axes[0][0].set_title("U at z= "+str(z_val)+"fm")


axes[0][1].plot(x,B_n,label="B_n")
axes[0][1].plot(x,B_p,label="B_p")
#axes[0][1].set_xlabel("r [fm]")
axes[0][1].set_ylabel("B [MeV.fm^2]")
axes[0][1].set_title("B at z= "+str(z_val)+"fm")
axes[0][1].legend()

axes[1][0].plot(x,Wr_n,label="Wr_n")
axes[1][0].plot(x,Wz_n,label="Wz_n")
axes[1][0].plot(x,Wr_p,label="Wr_p")
axes[1][0].plot(x,Wz_p,label="Wz_p")
axes[1][0].set_xlabel("r [fm]")
axes[1][0].set_ylabel("W [MeV.fm]")
axes[1][0].set_title("W at z= "+str(z_val)+"fm")
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
axes[1][1].set_title("W at r= "+str(df["r"].min())+"fm")
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