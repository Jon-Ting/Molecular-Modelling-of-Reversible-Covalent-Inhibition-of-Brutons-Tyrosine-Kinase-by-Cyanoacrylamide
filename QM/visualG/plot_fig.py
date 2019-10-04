"""
Python script to plot figures for QM data analysis.
Written by Jonathan Yik Chang Ting (Student ID: 44254124) for UQ BAdSc(Hons) Honours Project 2019.
Project Title: Molecular Modelling of Reversible Covalent Inhibition of Bruton’s Tyrosine Kinase by Cyanoacrylamides
Free to be distributed, edited and utilized for purposes beneficial to the world. :)
Though recognition of my effort is appreciated! XD
"""

import seaborn as sns
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from QM.run_gaussian.settings import DATA_PATH
from QM.visualG.plot_config import combination_dict, charge_list, DI_list

sns.set(context='paper', font_scale=1.5)
# sns.despine()  # Remove "chartjunk"

SINGLE_PLOT_WIDTH, SINGLE_PLOT_HEIGHT = 6, 4
FOUR_PLOTS_WIDTH, FOUR_PLOTS_HEIGHT = 8, 6
EIGHT_PLOTS_WIDTH, EIGHT_PLOTS_HEIGHT = 10, 14


def lin_reg(m, c, r2, xlimit, leg_loc, font_size="x-small"):
    x = np.linspace(xlimit[0], xlimit[1], 2)
    y = m * x + c
    plt.plot(x, y, '--r', label="$y = {0:.1f}x + {1:.1f}$\n$R^2 = {2:.2f}$".format(m, c, r2))
    plt.legend(loc=leg_loc, fontsize=font_size)
    plt.locator_params(axis='x', nbins=6)


if __name__ == "__main__":
    combination = "CombinationI"
    csv_file = "{0}/QM/Conformational_Analysis/Most_stable_conformers/{1}_Properties_Correlation.csv".format(DATA_PATH, combination)
    df = pd.read_csv(csv_file, index_col=0)
    print(df)
    prop_dict = combination_dict[combination[-1]]

    # TS S-C Distance
    chosen_dict = prop_dict["SC"]
    m, c, r2, x_axis, y_axis, leg, fontsize = chosen_dict["M"], chosen_dict["C"], chosen_dict["R2"], chosen_dict[
        "X-AXIS"], chosen_dict["Y-AXIS"], chosen_dict["LEG"], chosen_dict["FONTSIZE"]
    fig = plt.figure(figsize=(SINGLE_PLOT_WIDTH, SINGLE_PLOT_HEIGHT))
    fig.subplots_adjust(top=0.9, bottom=0.1, left=0.1, right=0.95)
    ax = sns.scatterplot(x="TS S-C Distance (A)", y="Addition Barrier (kcal/mol)", data=df, hue="Ligand", legend=False, s=100)
    ax.set(xlabel=r"TS S-C$\beta$ Distance ($\AA$)")
    for line in range(1, df.shape[0] + 1):
        ax.text(df["TS S-C Distance (A)"][line] + x_axis, df["Addition Barrier (kcal/mol)"][line] + y_axis,
                df["Ligand"][line], horizontalalignment='center', size='small', color='black')
    lin_reg(m, c, r2, plt.xlim(), leg, fontsize)
    plt.ylim(None, None); plt.xlim(None, None)
    plt.tight_layout()
    plt.savefig(r"{0}/TS S-C Distance".format(combination))

    # Ligand LUMO Energy
    chosen_dict = prop_dict["LUMO"]
    m, c, r2, x_axis, y_axis, leg, fontsize = chosen_dict["M"], chosen_dict["C"], chosen_dict["R2"], chosen_dict[
        "X-AXIS"], chosen_dict["Y-AXIS"], chosen_dict["LEG"], chosen_dict["FONTSIZE"]
    fig = plt.figure(figsize=(SINGLE_PLOT_WIDTH, SINGLE_PLOT_HEIGHT))
    fig.subplots_adjust(top=0.9, bottom=0.1, left=0.1, right=0.95)
    ax = sns.scatterplot(x="Ligand LUMO Energy (kcal/mol)", y="Addition Barrier (kcal/mol)", data=df, hue="Ligand", legend=False, s=100)
    for line in range(1, df.shape[0] + 1):
        ax.text(df["Ligand LUMO Energy (kcal/mol)"][line] + x_axis, df["Addition Barrier (kcal/mol)"][line] + y_axis,
                df["Ligand"][line], horizontalalignment='center', size='small', color='black')
    lin_reg(m, c, r2, plt.xlim(), leg, fontsize)
    plt.ylim(None, None); plt.xlim(None, None)
    plt.tight_layout()
    plt.savefig("{0}/Ligand LUMO Energy".format(combination))

    # Beta-Carbon Charges
    fig = plt.figure(figsize=(EIGHT_PLOTS_WIDTH, EIGHT_PLOTS_HEIGHT))
    fig.subplots_adjust(top=0.9, bottom=0.1, left=0.1, right=0.95, wspace=0.3, hspace=0.45)
    fig.suptitle(r"Correlation between $\beta$-Carbon Charge and Addition Barrier",
                 horizontalalignment='center', fontsize=16, weight='bold')

    for i, model in enumerate(charge_list):
        chosen_dict = prop_dict["Charge"][model]
        m, c, r2, x_axis, y_axis, leg, txt_align, fontsize, x_name = chosen_dict["M"], chosen_dict["C"], chosen_dict["R2"], chosen_dict[
            "X-AXIS"], chosen_dict["Y-AXIS"], chosen_dict["LEG"], chosen_dict["ALIGN"], chosen_dict["FONTSIZE"], chosen_dict["X-NAME"]
        ax1 = fig.add_subplot(4, 2, i + 1)
        p1 = sns.scatterplot(x=x_name, y="Addition Barrier (kcal/mol)", data=df,
                             hue="Ligand", legend=False, s=100)
        for line in range(1, df.shape[0] + 1):
            ax1.text(df[x_name][line] + x_axis, df["Addition Barrier (kcal/mol)"][line] + y_axis,
                     df["Ligand"][line], horizontalalignment=txt_align, size='small', color='black')
        lin_reg(m, c, r2, plt.xlim(), leg, fontsize)
    plt.savefig("{0}/Beta-Carbon Charges".format(combination))
    # plt.show()

    # Distortion/Interaction Analysis
    fig = plt.figure(figsize=(FOUR_PLOTS_WIDTH, FOUR_PLOTS_HEIGHT))
    # fig, axes = plt.subplots(nrows=2, ncols=2, sharex=False, sharey=True)
    fig.subplots_adjust(top=0.9, bottom=0.1, left=0.1, right=0.95, wspace=0.3, hspace=0.45)
    fig.suptitle("Correlation between Distortion/Interaction Energies and Addition Barrier", horizontalalignment='center', fontsize=16, weight='bold')

    for i, aspect in enumerate(DI_list):
        chosen_dict = prop_dict["DI"][aspect]
        m, c, r2, x_axis, y_axis, leg, txt_align, fontsize, x_name = chosen_dict["M"], chosen_dict["C"], chosen_dict["R2"], chosen_dict[
            "X-AXIS"], chosen_dict["Y-AXIS"], chosen_dict["LEG"], chosen_dict["ALIGN"], chosen_dict["FONTSIZE"], chosen_dict["X-NAME"]
        ax1 = fig.add_subplot(2, 2, i + 1)
        p1 = sns.scatterplot(x=x_name, y="Addition Barrier (kcal/mol)", data=df,
                             hue="Ligand", legend=False, s=100)
        for line in range(1, df.shape[0] + 1):
            ax1.text(df[x_name][line] + x_axis, df["Addition Barrier (kcal/mol)"][line] + y_axis,
                     df["Ligand"][line], horizontalalignment=txt_align, size='small', color='black')
        lin_reg(m, c, r2, plt.xlim(), leg, fontsize)
    plt.ylim(None, None); plt.xlim(None, None)
    plt.savefig("{0}/Distortion-Interaction Analysis".format(combination))
    # plt.show()
