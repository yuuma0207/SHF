import pandas as pd
import matplotlib.pyplot as plt
import os
import sys

# データを解析し、DataFrameに変換するための関数を定義
def parse_data(lines):
    steps = []
    neutron_data = []
    proton_data = []
    step = 0

    for line in lines[2:]:
        if '----' in line:
            step += 1
        else:
            neutron, proton = line.strip().split('|')
            neutron_values = [float(x) for x in neutron.split()]
            proton_values = [float(x) for x in proton.split()]

            steps.append(step)
            neutron_data.append(neutron_values)
            proton_data.append(proton_values)

    columns = ['one_E', 'Kin(u)', 'Kin(d)', 'U(u)', 'U(d)', 'B(u)', 'B(d)', 'W(u)', 'W(d)' , 'm', 'idx']
    df_neutron = pd.DataFrame(neutron_data, columns=columns)
    df_neutron['step'] = steps
    df_neutron['type'] = 'neutron'

    df_proton = pd.DataFrame(proton_data, columns=columns)
    df_proton['step'] = steps
    df_proton['type'] = 'proton'

    return pd.concat([df_neutron, df_proton])

# 各核子ごとに個別のDataFrameを作成し、それぞれのサブプロットにプロットする関数を定義
def plot_data_2x2_subplots(df):
    # パラメータのリスト
    #params = ['one_E', 'Kin(u)', 'Kin(d)', 'U(u)', 'U(d)', 'B(u)', 'B(d)', 'W(u)', 'W(d)']
    params = [('one_E',), ('Kin(u)', 'Kin(d)'), ('U(u)', 'U(d)'), ('B(u)', 'B(d)'), ('W(u)', 'W(d)')]

    # グラフのサイズを調整
    fig, axs = plt.subplots(2, 2, figsize=(14, 14))  # 2x2のサブプロット

    # 核子のタイプとインデックスに基づいてプロット
    for i in range(4):
        # 核子のタイプ（中性子または陽子）とインデックスを決定
        nucleon_type = 'neutron' if i < 2 else 'proton'
        nucleon_idx = 1 if i % 2 == 0 else 2

        # 特定の核子のデータを抽出
        df_nucleon = df[(df['type'] == nucleon_type) & (df['step'] > 0)].groupby("step").first().reset_index()
        # 対応するサブプロットの位置を計算
        subplot_row = i // 2
        subplot_col = i % 2

        for param in params:
            # 各パラメータをプロット
            ax = axs[subplot_row, subplot_col]
            #ax.plot(df_nucleon['step'], df_nucleon[param], ls="-",lw=2, label=param)
            sum_param = df_nucleon[list(param)].sum(axis=1)
            label = ' + '.join(param)
            ax.plot(df_nucleon['step'], sum_param, ls="-", lw=2, label=label)

        ax.set_title(f'{nucleon_type.capitalize()} {nucleon_idx}',fontsize=20)
        ax.legend()
        ax.set_xlabel('Iteration Numbder'           ,fontsize=14)
        ax.set_ylabel('Energy [MeV]',fontsize=14)
        ymin,ymax = plt.ylim()
        #ax.set_ylim(ymin,0)
        ax.grid(True, which="major", ls="--", color='0.65',lw=0.5)
        ax.grid(False, which="minor")
    fig.suptitle('One Particle Energy',fontsize=30)
    plt.tight_layout()
    return fig

def plot_energy_change_log_scale(df):
    fig, axs = plt.subplots(2, 2, figsize=(14, 14))  # 2x2のサブプロット

    # 核子のタイプとインデックスに基づいてプロット
    for i in range(4):
        # 核子のタイプ（中性子または陽子）とインデックスを決定
        nucleon_type = 'neutron' if i < 2 else 'proton'
        nucleon_idx = 1 if i % 2 == 0 else 2

        # 特定の核子のデータを抽出（最初のステップを除外）
        df_nucleon = df[(df['type'] == nucleon_type) & (df['step'] > 0)].groupby("step").first().reset_index()

        # 前のステップからのエネルギー変化量の絶対値を計算
        energy_change = df_nucleon['one_E'].diff().abs()

        # 対応するサブプロットの位置を計算
        subplot_row = i // 2
        subplot_col = i % 2

        # ログスケールでプロット
        axs[subplot_row, subplot_col].plot(df_nucleon['step'][1:], energy_change[1:], marker='o',markersize="1")
        axs[subplot_row, subplot_col].set_title(f'{nucleon_type.capitalize()} {nucleon_idx}',fontsize=20)
        axs[subplot_row, subplot_col].set_xlabel('Iteration Number',fontsize=14)
        axs[subplot_row, subplot_col].set_ylabel("$|E_{\\text{one}} - E_{\\text{one}_{\\text{prev}}}|$",fontsize=14)
        axs[subplot_row, subplot_col].set_yscale("log")
    fig.suptitle('Absolute Change in One Particle Energy',fontsize=30)
    plt.tight_layout()
    return fig

def plot_one_particle_energy(resfolder):
    filename = os.path.join(resfolder, "sort_information.txt")
    # filenameをパスのコンポーネントに分割
    path_components = filename.split('/')
    index = -2

    # 新しいパスを作成
    savefolder = '/'.join(path_components[:index+1] + ['plots'])+"/"
    if not os.path.exists(savefolder):
        os.makedirs(savefolder)
    with open(filename, 'r') as file:
        data = file.readlines()

    df = parse_data(data)
    # 核子のタイプとインデックスに基づいてプロット
    fig = plot_data_2x2_subplots(df)
    plt.savefig(savefolder+"one_particle_energy.png")
    plt.close(fig)
    fig = plot_energy_change_log_scale(df)
    plt.savefig(savefolder+"one_particle_energy_change.png")
    plt.close(fig)
    return

# メイン処理
def main():
    # python plot_nucleon_data.py ./../res/0122/125200/sort_information.txt
    filename=os.path.dirname(sys.argv[1])
    plot_one_particle_energy(filename)


if __name__ == "__main__":
    main()