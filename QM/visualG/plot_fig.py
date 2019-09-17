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
from QM.visualG.plot_config import combination_dict, charge_list

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

    # Ligand LUMO Energy
    # sorted_df = df.sort_values(by=["Ligand LUMO Energy (kcal/mol)"])
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

    raise
    chosen_dict = prop_dict["Charge"]["Mulliken"]
    m, c, r2, x_axis, y_axis, leg, fontsize = chosen_dict["M"], chosen_dict["C"], chosen_dict["R2"], chosen_dict[
        "X-AXIS"], chosen_dict["Y-AXIS"], chosen_dict["LEG"], chosen_dict["FONTSIZE"]
    ax1 = fig.add_subplot(4, 2, 1)
    p1 = sns.scatterplot(x="Mulliken Charge (e)", y="Addition Barrier (kcal/mol)", data=df,
                         hue="Ligand", legend=False, s=100)
    for line in range(1, df.shape[0] + 1):
        ax1.text(df["Mulliken Charge (e)"][line] + x_axis, df["Addition Barrier (kcal/mol)"][line] + y_axis,
                 df["Ligand"][line], horizontalalignment='right', size='small', color='black')
    lin_reg(m, c, r2, plt.xlim(), leg, fontsize)

    chosen_dict = prop_dict["Charge"]["NBO"]
    m, c, r2, x_axis, y_axis, leg, fontsize = chosen_dict["M"], chosen_dict["C"], chosen_dict["R2"], chosen_dict[
        "X-AXIS"], chosen_dict["Y-AXIS"], chosen_dict["LEG"], chosen_dict["FONTSIZE"]
    ax2 = fig.add_subplot(4, 2, 2)
    p2 = sns.scatterplot(x="NBO Charge (e)", y="Addition Barrier (kcal/mol)", data=df,
                         hue="Ligand", legend=False, s=100)
    for line in range(1, df.shape[0] + 1):
        ax2.text(df["NBO Charge (e)"][line] + x_axis, df["Addition Barrier (kcal/mol)"][line] + y_axis,
                 df["Ligand"][line], horizontalalignment='right', size='small', color='black')
    lin_reg(m, c, r2, plt.xlim(), leg, fontsize)

    chosen_dict = prop_dict["Charge"]["MK"]
    m, c, r2, x_axis, y_axis, leg, fontsize = chosen_dict["M"], chosen_dict["C"], chosen_dict["R2"], chosen_dict[
        "X-AXIS"], chosen_dict["Y-AXIS"], chosen_dict["LEG"], chosen_dict["FONTSIZE"]
    ax3 = fig.add_subplot(4, 2, 3)
    p3 = sns.scatterplot(x="Merz-Kollman Charge (e)", y="Addition Barrier (kcal/mol)", data=df, hue="Ligand",
                         legend=False, s=100)
    for line in range(1, df.shape[0] + 1):
        ax3.text(df["Merz-Kollman Charge (e)"][line] + x_axis, df["Addition Barrier (kcal/mol)"][line] + y_axis,
                 df["Ligand"][line], horizontalalignment='left', size='small', color='black')
    lin_reg(m, c, r2, plt.xlim(), leg, fontsize)

    chosen_dict = prop_dict["Charge"]["Hirshfeld"]
    m, c, r2, x_axis, y_axis, leg, fontsize = chosen_dict["M"], chosen_dict["C"], chosen_dict["R2"], chosen_dict[
        "X-AXIS"], chosen_dict["Y-AXIS"], chosen_dict["LEG"], chosen_dict["FONTSIZE"]
    ax4 = fig.add_subplot(4, 2, 4)
    p4 = sns.scatterplot(x="Hirshfeld Charge (e)", y="Addition Barrier (kcal/mol)", data=df, hue="Ligand",
                         legend=False, s=100)
    for line in range(1, df.shape[0] + 1):
        ax4.text(df["Hirshfeld Charge (e)"][line] + x_axis, df["Addition Barrier (kcal/mol)"][line] + y_axis,
                 df["Ligand"][line], horizontalalignment='left', size='small', color='black')
    lin_reg(m, c, r2, plt.xlim(), leg, fontsize)

    chosen_dict = prop_dict["Charge"]["CM5"]
    m, c, r2, x_axis, y_axis, leg, fontsize = chosen_dict["M"], chosen_dict["C"], chosen_dict["R2"], chosen_dict[
        "X-AXIS"], chosen_dict["Y-AXIS"], chosen_dict["LEG"], chosen_dict["FONTSIZE"]
    ax5 = fig.add_subplot(4, 2, 5)
    p5 = sns.scatterplot(x="CM5 Charge (e)", y="Addition Barrier (kcal/mol)", data=df, hue="Ligand",
                         legend=False, s=100)
    for line in range(1, df.shape[0] + 1):
        ax5.text(df["CM5 Charge (e)"][line] + x_axis, df["Addition Barrier (kcal/mol)"][line] + y_axis,
                 df["Ligand"][line], horizontalalignment='right', size='small', color='black')
    lin_reg(m, c, r2, plt.xlim(), leg, fontsize)

    chosen_dict = prop_dict["Charge"]["AIM"]
    m, c, r2, x_axis, y_axis, leg, fontsize = chosen_dict["M"], chosen_dict["C"], chosen_dict["R2"], chosen_dict[
        "X-AXIS"], chosen_dict["Y-AXIS"], chosen_dict["LEG"], chosen_dict["FONTSIZE"]
    ax6 = fig.add_subplot(4, 2, 6)
    p6 = sns.scatterplot(x="AIM Charge (e)", y="Addition Barrier (kcal/mol)", data=df, hue="Ligand",
                         legend=False, s=100)
    for line in range(1, df.shape[0] + 1):
        ax6.text(df["AIM Charge (e)"][line] + x_axis, df["Addition Barrier (kcal/mol)"][line] + y_axis,
                 df["Ligand"][line], horizontalalignment='right', size='small', color='black')
    lin_reg(m, c, r2, plt.xlim(), leg, fontsize)

    chosen_dict = prop_dict["Charge"]["ChelpG"]
    m, c, r2, x_axis, y_axis, leg, fontsize = chosen_dict["M"], chosen_dict["C"], chosen_dict["R2"], chosen_dict[
        "X-AXIS"], chosen_dict["Y-AXIS"], chosen_dict["LEG"], chosen_dict["FONTSIZE"]
    ax7 = fig.add_subplot(4, 2, 7)
    p5 = sns.scatterplot(x="ChelpG Charge (e)", y="Addition Barrier (kcal/mol)", data=df, hue="Ligand",
                         legend=False, s=100)
    for line in range(1, df.shape[0] + 1):
        ax7.text(df["ChelpG Charge (e)"][line] + x_axis, df["Addition Barrier (kcal/mol)"][line] + y_axis,
                 df["Ligand"][line], horizontalalignment='right', size='small', color='black')
    lin_reg(m, c, r2, plt.xlim(), leg, fontsize)

    chosen_dict = prop_dict["Charge"]["Omega"]
    m, c, r2, x_axis, y_axis, leg, fontsize = chosen_dict["M"], chosen_dict["C"], chosen_dict["R2"], chosen_dict[
        "X-AXIS"], chosen_dict["Y-AXIS"], chosen_dict["LEG"], chosen_dict["FONTSIZE"]
    ax8 = fig.add_subplot(4, 2, 8)
    p6 = sns.scatterplot(x="Electrophilicity Index", y="Addition Barrier (kcal/mol)", data=df, hue="Ligand",
                         legend=False, s=100)
    for line in range(1, df.shape[0] + 1):
        ax8.text(df["Electrophilicity Index"][line] + x_axis, df["Addition Barrier (kcal/mol)"][line] + y_axis,
                 df["Ligand"][line], horizontalalignment='right', size='small', color='black')
    lin_reg(m, c, r2, plt.xlim(), leg, fontsize)

    # plt.legend(p4, loc="center right", borderaxespad=0.1, title="Ligands")
    # plt.ylim(3, 15); plt.xlim(None, None)
    # plt.tight_layout()
    plt.savefig("{0}/Beta-Carbon Charges".format(combination))
    # plt.show()

    # Distortion/Interaction Analysis
    fig = plt.figure(figsize=(FOUR_PLOTS_WIDTH, FOUR_PLOTS_HEIGHT))
    # fig, axes = plt.subplots(nrows=2, ncols=2, sharex=False, sharey=True)
    fig.subplots_adjust(top=0.9, bottom=0.1, left=0.1, right=0.95, wspace=0.3, hspace=0.45)
    fig.suptitle("Correlation between Distortion/Interaction Energies and Addition Barrier", horizontalalignment='center', fontsize=16, weight='bold')

    chosen_dict = prop_dict["DI"]["Thiolate"]
    m, c, r2, x_axis, y_axis, leg, fontsize = chosen_dict["M"], chosen_dict["C"], chosen_dict["R2"], chosen_dict[
        "X-AXIS"], chosen_dict["Y-AXIS"], chosen_dict["LEG"], chosen_dict["FONTSIZE"]
    ax1 = fig.add_subplot(2, 2, 1)
    p1 = sns.scatterplot(x="Thiolate Distortion Energy (kcal/mol)", y="Addition Barrier (kcal/mol)", data=df, hue="Ligand", legend=False, s=100)
    for line in range(1, df.shape[0] + 1):
        ax1.text(df["Thiolate Distortion Energy (kcal/mol)"][line] + x_axis, df["Addition Barrier (kcal/mol)"][line] + y_axis,
                df["Ligand"][line], horizontalalignment='right', size='small', color='black')
    # m, c, r2 = 17.241, 7.6018, 0.798
    lin_reg(m, c, r2, plt.xlim(), leg, fontsize)

    chosen_dict = prop_dict["DI"]["Ligand"]
    m, c, r2, x_axis, y_axis, leg, fontsize = chosen_dict["M"], chosen_dict["C"], chosen_dict["R2"], chosen_dict[
        "X-AXIS"], chosen_dict["Y-AXIS"], chosen_dict["LEG"], chosen_dict["FONTSIZE"]
    ax2 = fig.add_subplot(2, 2, 2)
    p2 = sns.scatterplot(x="Ligand Distortion Energy (kcal/mol)", y="Addition Barrier (kcal/mol)", data=df, hue="Ligand", legend=False, s=100)
    for line in range(1, df.shape[0] + 1):
        ax2.text(df["Ligand Distortion Energy (kcal/mol)"][line] + x_axis, df["Addition Barrier (kcal/mol)"][line] + y_axis,
                df["Ligand"][line], horizontalalignment='left', size='small', color='black')
    # m, c, r2 = 0.8082, 5.3355, 0.9760
    lin_reg(m, c, r2, plt.xlim(), "best")

    chosen_dict = prop_dict["DI"]["Activation"]
    m, c, r2, x_axis, y_axis, leg, fontsize = chosen_dict["M"], chosen_dict["C"], chosen_dict["R2"], chosen_dict[
        "X-AXIS"], chosen_dict["Y-AXIS"], chosen_dict["LEG"], chosen_dict["FONTSIZE"]
    ax3 = fig.add_subplot(2, 2, 3)
    p3 = sns.scatterplot(x="Activation Energy (kcal/mol)", y="Addition Barrier (kcal/mol)", data=df, hue="Ligand", legend=False, s=100)
    for line in range(1, df.shape[0] + 1):
        ax3.text(df["Activation Energy (kcal/mol)"][line] + x_axis, df["Addition Barrier (kcal/mol)"][line] + y_axis,
                df["Ligand"][line], horizontalalignment='right', size='small', color='black')
    # m, c, r2 = 0.9013, 10.439, 0.8943
    lin_reg(m, c, r2, plt.xlim(), "upper left")

    chosen_dict = prop_dict["DI"]["Interaction"]
    m, c, r2, x_axis, y_axis, leg, fontsize = chosen_dict["M"], chosen_dict["C"], chosen_dict["R2"], chosen_dict[
        "X-AXIS"], chosen_dict["Y-AXIS"], chosen_dict["LEG"], chosen_dict["FONTSIZE"]
    ax4 = fig.add_subplot(2, 2, 4)
    p4 = sns.scatterplot(x="Interaction Energy (kcal/mol)", y="Addition Barrier (kcal/mol)", data=df, hue="Ligand", legend=False, s=100)
    for line in range(1, df.shape[0] + 1):
        ax4.text(df["Interaction Energy (kcal/mol)"][line] + x_axis, df["Addition Barrier (kcal/mol)"][line] + y_axis,
                df["Ligand"][line], horizontalalignment='left', size='small', color='black')
    # m, c, r2 = -1.6397, 0.1818, 0.428
    lin_reg(m, c, r2, plt.xlim(), "lower left")

    # plt.legend(p4, loc="center right", borderaxespad=0.1, title="Ligands")
    plt.ylim(None, None); plt.xlim(None, None)
    # plt.tight_layout()
    plt.savefig("{0}/Distortion-Interaction Analysis".format(combination))
    # plt.show()



