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

sns.set(context='paper', font_scale=1.5)
# sns.despine()  # Remove "chartjunk"

SINGLE_PLOT_WIDTH, SINGLE_PLOT_HEIGHT = 6, 4
FOUR_PLOTS_WIDTH, FOUR_PLOTS_HEIGHT = 8, 6
SIX_PLOTS_WIDTH, SIX_PLOTS_HEIGHT = 10, 8


def lin_reg(m, c, r2, xlimit, leg_loc, font_size="x-small"):
    x = np.linspace(xlimit[0], xlimit[1], 2)
    y = m * x + c
    plt.plot(x, y, '--r', label="$y = {0}x + {1}$\n\t$R^2 = {2}$".format(m, c, r2))
    plt.legend(loc=leg_loc, fontsize=font_size)


if __name__ == "__main__":
    csv_file = "{0}/QM/Conformational_Analysis/Most_stable/Properties_Correlation.csv".format(DATA_PATH)
    df = pd.read_csv(csv_file, index_col=0)
    print(df)

    # Ligand LUMO Energy
    # sorted_df = df.sort_values(by=["Ligand LUMO Energy (kcal/mol)"])
    fig = plt.figure(figsize=(SINGLE_PLOT_WIDTH, SINGLE_PLOT_HEIGHT))
    fig.subplots_adjust(top=0.9, bottom=0.1, left=0.1, right=0.95)
    ax = sns.scatterplot(x="Ligand LUMO Energy (kcal/mol)", y="Addition Barrier (kcal/mol)", data=df, hue="Ligand", legend=False, s=100)
    for line in range(1, df.shape[0] + 1):
        ax.text(df["Ligand LUMO Energy (kcal/mol)"][line], df["Addition Barrier (kcal/mol)"][line] + 0.4,
                df["Ligand"][line], horizontalalignment='center', size='small', color='black')
    m, c, r2 = 0.6276, 19.017, 0.7402
    lin_reg(m, c, r2, plt.xlim(), "upper left", font_size="medium")
    plt.ylim(None, None); plt.xlim(None, None)
    plt.tight_layout()
    plt.savefig("Ligand LUMO Energy")

    # Beta-Carbon Charges
    fig = plt.figure(figsize=(SIX_PLOTS_WIDTH, SIX_PLOTS_HEIGHT))
    fig.subplots_adjust(top=0.9, bottom=0.1, left=0.1, right=0.95, wspace=0.3, hspace=0.45)
    fig.suptitle(r"Correlation between $\beta$-Carbon Charge and Addition Barrier",
                 horizontalalignment='center', fontsize=16, weight='bold')
    ax1 = fig.add_subplot(3, 2, 1)
    p1 = sns.scatterplot(x="Mulliken Charge (e)", y="Addition Barrier (kcal/mol)", data=df,
                         hue="Ligand", legend=False, s=100)
    for line in range(1, df.shape[0] + 1):
        ax1.text(df["Mulliken Charge (e)"][line]-0.03, df["Addition Barrier (kcal/mol)"][line],
                 df["Ligand"][line], horizontalalignment='right', size='small', color='black')
    m, c, r2 = -9.5378, 8.4412, 0.6933
    lin_reg(m, c, r2, plt.xlim(), "upper right")

    ax2 = fig.add_subplot(3, 2, 2)
    p2 = sns.scatterplot(x="NBO Charge (e)", y="Addition Barrier (kcal/mol)", data=df,
                         hue="Ligand", legend=False, s=100)
    for line in range(1, df.shape[0] + 1):
        ax2.text(df["NBO Charge (e)"][line]-0.01, df["Addition Barrier (kcal/mol)"][line],
                 df["Ligand"][line], horizontalalignment='right', size='small', color='black')
    m, c, r2 = -13.261, 7.3075, 0.2655
    lin_reg(m, c, r2, plt.xlim(), "upper right")

    ax3 = fig.add_subplot(3, 2, 3)
    p3 = sns.scatterplot(x="Merz-Kollman Charge (e)", y="Addition Barrier (kcal/mol)", data=df, hue="Ligand",
                         legend=False, s=100)
    for line in range(1, df.shape[0] + 1):
        ax3.text(df["Merz-Kollman Charge (e)"][line]+0.01, df["Addition Barrier (kcal/mol)"][line],
                 df["Ligand"][line], horizontalalignment='left', size='small', color='black')
    m, c, r2 = 7.8488, 4.5388, 0.1579
    lin_reg(m, c, r2, plt.xlim(), "lower right")

    ax4 = fig.add_subplot(3, 2, 4)
    p4 = sns.scatterplot(x="Hirshfeld CM5 Charge (e)", y="Addition Barrier (kcal/mol)", data=df, hue="Ligand",
                         legend=False, s=100)
    for line in range(1, df.shape[0] + 1):
        ax4.text(df["Hirshfeld CM5 Charge (e)"][line]+0.03, df["Addition Barrier (kcal/mol)"][line],
                 df["Ligand"][line], horizontalalignment='left', size='small', color='black')
    m, c, r2 = -5.7879, 9.3252, 0.4627
    lin_reg(m, c, r2, plt.xlim(), "upper right")

    ax5 = fig.add_subplot(3, 2, 5)
    p5 = sns.scatterplot(x="ChelpG Charge (e)", y="Addition Barrier (kcal/mol)", data=df, hue="Ligand",
                         legend=False, s=100)
    for line in range(1, df.shape[0] + 1):
        ax5.text(df["ChelpG Charge (e)"][line]-0.025, df["Addition Barrier (kcal/mol)"][line],
                 df["Ligand"][line], horizontalalignment='right', size='small', color='black')
    m, c, r2 = -20.086, 21.333, 0.3557
    lin_reg(m, c, r2, plt.xlim(), "upper left")

    ax6 = fig.add_subplot(3, 2, 6)
    p6 = sns.scatterplot(x="AIM Charge (e)", y="Addition Barrier (kcal/mol)", data=df, hue="Ligand",
                         legend=False, s=100)
    for line in range(1, df.shape[0] + 1):
        ax6.text(df["AIM Charge (e)"][line]-0.008, df["Addition Barrier (kcal/mol)"][line],
                 df["Ligand"][line], horizontalalignment='right', size='small', color='black')
    m, c, r2 = -63.398, 10.76, 0.1817
    lin_reg(m, c, r2, plt.xlim(), "upper right")

    # plt.legend(p4, loc="center right", borderaxespad=0.1, title="Ligands")
    # plt.ylim(3, 15); plt.xlim(None, None)
    # plt.tight_layout()
    plt.savefig("Beta-Carbon Charges")
    # plt.show()

    # Distortion/Interaction Analysis
    fig = plt.figure(figsize=(FOUR_PLOTS_WIDTH, FOUR_PLOTS_HEIGHT))
    # fig, axes = plt.subplots(nrows=2, ncols=2, sharex=False, sharey=True)
    fig.subplots_adjust(top=0.9, bottom=0.1, left=0.1, right=0.95, wspace=0.3, hspace=0.45)
    fig.suptitle("Correlation between Distortion/Interaction Energies and Addition Barrier", horizontalalignment='center', fontsize=16, weight='bold')
    ax1 = fig.add_subplot(2, 2, 1)
    p1 = sns.scatterplot(x="Thiolate Distortion Energy (kcal/mol)", y="Addition Barrier (kcal/mol)", data=df, hue="Ligand", legend=False, s=100)
    for line in range(1, df.shape[0] + 1):
        ax1.text(df["Thiolate Distortion Energy (kcal/mol)"][line]-0.03, df["Addition Barrier (kcal/mol)"][line],
                df["Ligand"][line], horizontalalignment='right', size='small', color='black')
    m, c, r2 = 18.279, 4.4489, 0.711
    lin_reg(m, c, r2, plt.xlim(), "upper left")

    ax2 = fig.add_subplot(2, 2, 2)
    p2 = sns.scatterplot(x="Ligand Distortion Energy (kcal/mol)", y="Addition Barrier (kcal/mol)", data=df, hue="Ligand", legend=False, s=100)
    for line in range(1, df.shape[0] + 1):
        ax2.text(df["Ligand Distortion Energy (kcal/mol)"][line]+5, df["Addition Barrier (kcal/mol)"][line],
                df["Ligand"][line], horizontalalignment='left', size='small', color='black')
    m, c, r2 = 0.0194, 8.4567, 0.0323
    lin_reg(m, c, r2, plt.xlim(), "best")

    ax3 = fig.add_subplot(2, 2, 3)
    p3 = sns.scatterplot(x="Activation Energy (kcal/mol)", y="Addition Barrier (kcal/mol)", data=df, hue="Ligand", legend=False, s=100)
    for line in range(1, df.shape[0] + 1):
        ax3.text(df["Activation Energy (kcal/mol)"][line]-0.5, df["Addition Barrier (kcal/mol)"][line],
                df["Ligand"][line], horizontalalignment='right', size='small', color='black')
    m, c, r2 = 0.9708, 5.5026, 0.9026
    lin_reg(m, c, r2, plt.xlim(), "upper left")

    ax4 = fig.add_subplot(2, 2, 4)
    p4 = sns.scatterplot(x="Interaction Energy (kcal/mol)", y="Addition Barrier (kcal/mol)", data=df, hue="Ligand", legend=False, s=100)
    for line in range(1, df.shape[0] + 1):
        ax4.text(df["Interaction Energy (kcal/mol)"][line]+5, df["Addition Barrier (kcal/mol)"][line],
                df["Ligand"][line], horizontalalignment='left', size='small', color='black')
    m, c, r2 = 0.0088, 8.7397, 0.0068
    lin_reg(m, c, r2, plt.xlim(), "lower left")

    # plt.legend(p4, loc="center right", borderaxespad=0.1, title="Ligands")
    plt.ylim(None, None); plt.xlim(None, None)
    # plt.tight_layout()
    plt.savefig("Distortion-Interaction Analysis")
    # plt.show()


