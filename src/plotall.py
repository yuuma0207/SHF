import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import os
import sys
import argparse
import re
import glob

cur_dir = os.path.dirname(__file__)
utils_path = os.path.join(cur_dir, "./../utils")
sys.path.append(utils_path)

from plot_d_wf import plot_d_wf_all
from plot_wf import plot_wf_all
from plot_density import plot_density_all
from plot_dens_merge import plot_merge_dens_all
from plot_potential import plot_potential_all
from plot_pot_merge import plot_merge_pot_all
from plot_meanfield_E import plot_meanfield_E
from plot_nucleon_data import plot_one_particle_energy

def plot_d_wf(resfolder, z_val, to_write):
    plot_d_wf_all(resfolder, z_val)
    if to_write:
        print("plot_d_wf")
    return

def plot_density(resfolder, z_val, to_write):
    plot_density_all(resfolder, z_val)
    plot_merge_dens_all(resfolder, z_val)
    if to_write:
        print("plot_density")
    return

def plot_wf(resfolder, z_val, to_write):
    plot_wf_all(resfolder, z_val)
    if to_write:
        print("plot_wf")
    return

def plot_potential(resfolder, z_val, to_write):
    plot_potential_all(resfolder, z_val)
    plot_merge_pot_all(resfolder, z_val)
    if to_write:
        print("plot_potential")
    return

def plot_energy(resfolder, to_write):
    plot_meanfield_E(resfolder)
    plot_one_particle_energy(resfolder)
    if to_write:
        print("plot_energy")
    return

# python plotall.py ./../res/0117/132937/ all 10
def main():
    parser = argparse.ArgumentParser(description='Plotting script')
    parser.add_argument('foldername', help='Result folder for the plot data. e.g. ./../res/0117/132937/')
    parser.add_argument("-p",'--plot_type', help='Type of the plot. One of [density,wf,potential,energy,all]',default="all")
    parser.add_argument("-z",'--z_val', type=float, help='Optional z value', default=10.0)
    parser.add_argument("-w",'--write', help='print which type of plot', action='store_true')

    args = parser.parse_args()

    resfolder = args.foldername
    plot_type = args.plot_type
    z_val = args.z_val
    to_write = args.write

    if plot_type == "density":
        plot_density(resfolder, z_val, to_write)
    elif plot_type == "wf":
        plot_wf(resfolder, z_val, to_write)
        plot_d_wf(resfolder, z_val, to_write)
    elif plot_type == "potential":
        plot_potential(resfolder, z_val, to_write)
    elif plot_type == "energy":
        plot_energy(resfolder, to_write)
    elif plot_type == "all":
        plot_wf(resfolder, z_val, to_write)
        plot_density(resfolder, z_val, to_write)
        plot_d_wf(resfolder, z_val, to_write)
        plot_potential(resfolder, z_val, to_write)
        plot_energy(resfolder, to_write)
    else:
        print("Usage: type must be one of [density,wf,potential,energy,all]")
        sys.exit(1)

if __name__ == "__main__":
    main()

