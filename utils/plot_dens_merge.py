import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FormatStrFormatter
from scipy.interpolate import griddata
import os
import sys
import re
import glob

def plot_merge_dens_all(resfolder, z_val):
    folder_path = os.path.join(resfolder,'densities/')
    txt_files = glob.glob(os.path.join(folder_path, 'density_and_derivative[0-9]*.txt'))
    for file in txt_files:
        plot_with_slices(file, z_val)
    return

# Function to read the data and return structured DataFrame
def read_data(filename):
    data = pd.read_csv(filename, sep="\s+", header=None)
    data.columns = data.iloc[0] # Using the second row as the header
    data = data[1:] # Dropping the first row which is the header
    data = data.apply(pd.to_numeric) # Converting all data to numeric type
    return data

# 指定されたZ座標で、特定のレベルに対応するR座標を取得
def get_r_coordinates_for_z(cntr, specified_z):
    r_coordinates = []

    for i, seg in enumerate(cntr.allsegs):
        level = cntr.levels[i]
        closest_point = None
        min_distance = float('inf')

        for path in seg:
            for point in path:
                distance = abs(point[1] - specified_z)
                if distance < min_distance:
                    min_distance = distance
                    closest_point = point

        if closest_point is not None:
            r_coordinates.append((level, closest_point[0]))

    return r_coordinates

def get_z_coordinates_for_r(cntr, specified_r, threshold=0.1):
    z_coordinates = []

    for i, seg in enumerate(cntr.allsegs):
        for path in seg:
            # pathは等高線の一部の(x, y)座標のリスト
            # ここではxがrに対応していると仮定
            # 指定されたrに近い全ての点を見つけます
            close_points = [point for point in path if abs(point[0] - specified_r) < threshold]

            # しきい値内にある全てのz座標を追加
            for point in close_points:
                z_coordinates.append(point[1])
    
    return z_coordinates
def plot_merge(U_grid, U_r_slice, U_z_slice, particle_type,
                z_val, r_min, savefolder, range_rmin, range_rmax, range_zmin, range_zmax,
                z_unique, r_unique, name, R, Z, num):
    matches = re.findall(r'\$.*?\$|\S+', name)
    fig = plt.figure(figsize=(8, 8))
    plt.suptitle(f'{particle_type.capitalize()} '+matches[1]+' in '+num+" Steps",fontsize=25)
    gs = gridspec.GridSpec(2, 2, 
                            width_ratios=[2, 5], height_ratios=[5, 2],
                            hspace=0,wspace=0)

    # Create all subplots without sharing axes
    ax_contour = plt.subplot(gs[0, 1])
    ax_r_slice = plt.subplot(gs[1, 1])
    ax_z_slice = plt.subplot(gs[0, 0])
    ax_cbar = plt.subplot(gs[1, 0])
    #ax_r_slice.yaxis.set_major_locator(MaxNLocator(nbins=2))
    #ax_z_slice.xaxis.set_major_locator(MaxNLocator(nbins=3))

    # Contour plot
    levels = np.linspace(np.min(U_grid), np.max(U_grid), 200)
    if(np.min(U_grid) >= np.max(U_grid)):
        min_value = -1e-10  # ここに適切な最小値を設定
        max_value = 1e-10  # densの最大値を使用
        levels = np.linspace(min_value, max_value, 200)
    contour = ax_contour.contourf(R, Z, U_grid, 
                                levels=levels, cmap="plasma")
    levels = np.linspace(np.min(U_grid), np.max(U_grid), 5)[1:-1]
    if(np.min(U_grid) >= np.max(U_grid)):
        min_value = -1e-10  # ここに適切な最小値を設定
        max_value = 1e-10  # densの最大値を使用
        levels = np.linspace(min_value, max_value, 5)[1:-1]
    cntr = ax_contour.contour(R, Z, U_grid, levels=levels, 
                                colors="black", linestyles="dotted")
    ax_contour.scatter(r_min, z_val, marker="x", color="green", s=200)
    r_coordinates = get_r_coordinates_for_z(cntr,z_val)
    z_coordinates = get_z_coordinates_for_r(cntr,r_min)
    r_lim = ax_contour.get_xlim()
    r_min = r_lim[0]
    z_lim = ax_contour.get_ylim()
    z_min = z_lim[0]

    # Plot dashed lines on the contour plot
    pre_r=r_min
    lw = 2
    #dashed_style = "dotted"
    dashed_style = (0, (5, 10))  # 線分の長さと間隔の長さをタプルで定義 (線の長さ, 間隔の長さ)
    #alpha=0.3
    alpha=1
    color="black"
    for level, r in r_coordinates:
        ax_contour.plot([r, r], [z_min, z_val], 
                        linestyle=dashed_style, alpha=alpha, color=color)  # Vertical dashed line at the R coordinate
    #cntr_labels = ax_contour.clabel(cntr, inline=False, fontsize=12, colors="black", fmt='%1.2f')

    #for txt in cntr_labels:
    #    x, y = txt.get_position()
    #    angle_rad = np.arctan2(y - z_val, x)  # z_fixedからの角度を計算
    #    theta = np.pi/4
    #    off = 0.1
    #    dist = np.sqrt(x**2 + (y-z_val)**2)+off
    #    new_x = dist * np.cos(-angle_rad)  # 新しいX位置
    #    new_y = z_val + dist * np.sin(-angle_rad)  # 新しいY位置
    #    new_x = dist * np.cos(theta)  # 新しいX位置
    #    new_y = z_val + dist * np.sin(theta)  # 新しいY位置
    #    txt.set_position((new_x, new_y))
    #    txt.set_rotation(np.degrees(-np.pi/4))

    #ax_contour.title.set_text(f'{particle_type.capitalize()} '+name+' contour')
    #ax_contour.set_aspect("equal")
    ax_contour.set_xlim([range_rmin, range_rmax])
    ax_contour.set_ylim([range_zmin, range_zmax])
    #cbar_position = fig.add_axes([0.9, 0.1, 0.03, 0.8]) # 例: [left, bottom, width, height]
    #cbar = fig.colorbar(contour, ax=cbar_position,format=FormatStrFormatter('%.3f'))
    #cbar.ax.set_position(cbar_position.get_position())
    #ax_contour.set_xlim(0,6)
    #ax_contour.set_ylim(7,13)
    ax_contour.spines["left"].set_visible(False)
    ax_contour.yaxis.set_ticks_position("none")


    # r slice plot
    ax_r_slice.plot(r_unique, U_r_slice, 'b-', label=f"z={round(z_val,2)} fm",lw=lw)
    ax_r_slice.set_xlabel('r [fm]',fontsize=14)
    matches = re.findall(r'\$.*?\$|\S+', name)
    ax_r_slice.set_ylabel(matches[1]+" "+matches[2],fontsize=14)
    ax_r_slice.set_ylim([min(U_r_slice), max(U_r_slice)])
    ax_r_slice.legend()
    ax_r_slice.set_yticks([0,round(max(U_r_slice)/2,2),U_r_slice[-2]])
    ax_r_slice.yaxis.set_label_position("right")
    ax_r_slice.yaxis.tick_right()
    ax_r_slice.legend()
    ax_r_slice.spines['top'].set_visible(False)
    ymax = ax_r_slice.get_ylim()[1]

    # 縦線
    for _, r in r_coordinates:
        idx = np.argmin(np.abs(r_unique-r))
        ymin = min(U_r_slice[idx],ymax)
        r = r_unique[idx]
        ax_r_slice.plot([r, r], [ymin,ymax],
                        color="black", linestyle=dashed_style,alpha=alpha)

    # Remove the gaps by adjusting the axes positions directly
    # Get the bounding box of the contour plot
    bbox = ax_contour.get_position()
    left, bottom, width, height = bbox.bounds

    # Set the position of the r slice plot to directly below the contour plot
    ax_r_slice.set_position([left, bottom - height / 2.5, width, height / 2.5])


    # 横線
    xmax = ax_z_slice.get_xlim()[1]
    for z in z_coordinates:
        idx = np.argmin(np.abs(z_unique-z))
        x_min = min(U_z_slice[idx],xmax)
        z = z_unique[idx]
        ax_z_slice.plot([x_min, xmax], [z,z], 
                        color="black", linestyle=dashed_style,alpha=alpha)
    # z slice plot
    ax_z_slice.plot(U_z_slice, z_unique, 'r-', label=f'r={round(r_min,2)} fm',lw=lw)
    matches = re.findall(r'\$.*?\$|\S+', name)
    ax_z_slice.set_xlabel(matches[1]+" "+matches[2],fontsize=14)
    #ax_z_slice.set_xticks([0,-20,-40])
    ax_z_slice.set_xlim([min(U_z_slice), max(U_z_slice)])
    ax_z_slice.set_ylabel('z [fm]',fontsize=14)
    ax_z_slice.legend()
    #ax_z_slice.title.set_text(f'Potential Slice at r={r_min} (left)')
    #ax_z_slice.invert_xaxis()  # Invert x-axis to align with the contour plot
    # Set the position of the z slice plot to directly to the left of the contour plot
    ax_z_slice.set_position([left - width / 5.0, bottom, width / 5.0, height])

    ax_contour.set_xlim([range_rmin, range_rmax])
    ax_r_slice.set_xlim([range_rmin, range_rmax])

    ax_contour.set_ylim([range_zmin, range_zmax])
    ax_z_slice.set_ylim([range_zmin, range_zmax])
    plt.tight_layout(pad=3.0, h_pad=0.5, w_pad=0.5)
    plt.subplots_adjust(right=0.87)
    cbar = fig.colorbar(contour, cax=ax_cbar,format=FormatStrFormatter('%.3f'))
    

    pos_ax = ax_contour.get_position()
    right = pos_ax.x0 + pos_ax.width
    pos = [right+0.01, 0.4, 0.02, 0.4] # 例: [left, bottom, width, height]
    cbar.ax.set_position(pos)

    # Hide the inner axes numbers of the contour plot
    plt.setp(ax_contour.get_xticklabels(), visible=False)
    plt.setp(ax_contour.get_yticklabels(), visible=False)
    # Adjust layout and spacing
    #plt.tight_layout()
    #plt.subplots_adjust(wspace=0.2, hspace=0.2)
    #plt.show()
    plt.savefig(savefolder+f"{particle_type.capitalize()}_"+name.split()[0]+"_merge.png")
    plt.close(fig)
    return
# Function to create contour and slice plots for both neutron and proton
def plot_with_slices(file_path,z_val):
    pattern = r'(\d{1,})\.txt$'
    num = re.search(pattern,file_path).group(1)
    # filenameをパスのコンポーネントに分割
    path_components = file_path.split('/')

    # 'orbital_wavefunction' が存在するインデックスを探す
    index = -3

    # 新しいパスを作成
    savefolder = '/'.join(path_components[:index+1] + ['plots'] + [path_components[index+1]])+"/"
    if not os.path.exists(savefolder):
        os.makedirs(savefolder)
    for file in os.listdir(savefolder):
        if(file.startswith(num) and file.endswith("merge.png")):
            return None
    savefolder = savefolder+ num
    # Read the data
    data = read_data(file_path)

    # Extracting r, z, and potential values
    r = data['r']
    z = data['z']
    rho_n = data['rho_n']
    rho_p = data['rho_p']
    tau_n = data['tau_n']
    tau_p = data['tau_p']

    # Create grid values for contour plot
    r_unique = np.sort(r.unique())
    z_unique = np.sort(z.unique())
    R, Z = np.meshgrid(r_unique, z_unique)

    # Interpolate unstructured D-dimensional data.
    rho_n_grid = griddata((r, z), rho_n, (R, Z), method='cubic')
    rho_p_grid = griddata((r, z), rho_p, (R, Z), method='cubic')
    tau_n_grid = griddata((r, z), tau_n, (R, Z), method='cubic')
    tau_p_grid = griddata((r, z), tau_p, (R, Z), method='cubic')

    # Select specific z value and r value for the one-dimensional plots
    #z_val = z_val # Selecting the median z value for demonstration
    z_val = min(z_unique, key=lambda x: abs(x - z_val))
    r_min = r_unique[0]  # Selecting the minimum r value for demonstration

    # Extract data at specific z value and r value
    rho_n_z_slice = rho_n_grid[:, r_unique == r_min].flatten()
    rho_p_z_slice = rho_p_grid[:, r_unique == r_min].flatten()
    rho_n_r_slice = rho_n_grid[z_unique == z_val, :].flatten()
    rho_p_r_slice = rho_p_grid[z_unique == z_val, :].flatten()
    tau_n_z_slice = tau_n_grid[:, r_unique == r_min].flatten()
    tau_p_z_slice = tau_p_grid[:, r_unique == r_min].flatten()
    tau_n_r_slice = tau_n_grid[z_unique == z_val, :].flatten()
    tau_p_r_slice = tau_p_grid[z_unique == z_val, :].flatten()

    plot_data=[(rho_n_grid, rho_n_r_slice, rho_n_z_slice, "neutron", z_val, r_min, "rho $\\rho$ [fm$^{-3}$]"),
                (rho_p_grid, rho_p_r_slice, rho_p_z_slice, "proton", z_val, r_min, "rho $\\rho$ [fm$^{-3}$]"),
                (tau_n_grid, tau_n_r_slice, tau_n_z_slice, "neutron", z_val, r_min,"tau $\\tau$ [fm$^{-5}$]"),
                (tau_p_grid, tau_p_r_slice, tau_p_z_slice, "proton", z_val, r_min, "tau $\\tau$ [fm$^{-5}$]")]

    plot_range = 20
    absolute_differences = np.abs(r_unique - 0)
    # Find the index of the smallest difference
    closest_index = np.argmin(absolute_differences)
    range_rmin = r_unique[closest_index]
    # Get the closest value from 'r_unique'
    absolute_differences = np.abs(r_unique - plot_range)
    # Find the index of the smallest difference
    closest_index = np.argmin(absolute_differences)
    # Get the closest value from 'r_unique'
    range_rmax = r_unique[closest_index]
    z_range = [z_val-plot_range/2,z_val+plot_range/2]
    absolute_differences = np.abs(z_unique - z_range[0])
    # Find the index of the smallest difference
    closest_index = np.argmin(absolute_differences)
    # Get the closest value from 'r_unique'
    range_zmin = z_unique[closest_index]

    absolute_differences = np.abs(z_unique - z_range[1])
    # Find the index of the smallest difference
    closest_index = np.argmin(absolute_differences)
    # Get the closest value from 'r_unique'
    range_zmax = z_unique[closest_index]
    
    for U_grid, U_r_slice, U_z_slice, particle_type, z_val, r_min, name in plot_data:
        plot_merge(U_grid, U_r_slice, U_z_slice, particle_type,
                    z_val, r_min, savefolder, range_rmin, range_rmax, range_zmin, range_zmax,
                    z_unique, r_unique, name, R, Z, num)

def main():
    filename=os.path.dirname(sys.argv[1])
    z_val = float(sys.argv[2])
    plot_merge_dens_all(filename,z_val)
    return

# Main execution
if __name__ == "__main__":
    main()
