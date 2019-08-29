"""
Python script to plot figures for QM data analysis.
Written by Jonathan Yik Chang Ting (Student ID: 44254124) for UQ BAdSc(Hons) Honours Project 2019.
Project Title: Molecular Modelling of Reversible Covalent Inhibition of Bruton’s Tyrosine Kinase by Cyanoacrylamides
Free to be distributed, edited and utilized for purposes beneficial to the world. :)
Though recognition of my effort is appreciated! XD
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from run_gaussian.settings import DATA_PATH

sns.set(context='paper', font_scale=1.5)
# sns.despine()  # Remove "chartjunk"

SINGLE_PLOT_WIDTH, SINGLE_PLOT_HEIGHT = 6, 4
SUBPLOTS_WIDTH, SUBPLOTS_HEIGHT = 10, 6


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
    plt.ylim(None, None); plt.xlim(None, None)
    plt.tight_layout()
    plt.savefig("Ligand LUMO Energy")

    # Beta-Carbon Charges
    fig = plt.figure(figsize=(SINGLE_PLOT_WIDTH, SINGLE_PLOT_HEIGHT))
    fig.subplots_adjust(top=0.9, bottom=0.1, left=0.1, right=0.95)
    ax = sns.scatterplot(x="B-Carbon Charge (e)", y="Addition Barrier (kcal/mol)", data=df, hue="Ligand", legend=False, s=100)
    for line in range(1, df.shape[0] + 1):
        ax.text(df["B-Carbon Charge (e)"][line], df["Addition Barrier (kcal/mol)"][line] + 0.4,
                df["Ligand"][line], horizontalalignment='center', size='small', color='black')
    plt.ylim(None, None); plt.xlim(None, None)
    plt.tight_layout()
    plt.savefig("Beta-Carbon Charges")

    # Distortion/Interaction Analysis
    fig = plt.figure(figsize=(SUBPLOTS_WIDTH, SUBPLOTS_HEIGHT))
    # fig, axes = plt.subplots(nrows=2, ncols=2, sharex=False, sharey=True)
    fig.subplots_adjust(top=0.9, bottom=0.1, left=0.1, right=0.95, wspace=0.3, hspace=0.45)
    fig.suptitle("Correlation between Distortion/Interaction Energies and Addition Barrier", horizontalalignment='center', fontsize=16, weight='bold')
    ax1 = fig.add_subplot(2, 2, 1)
    p1 = sns.scatterplot(x="Thiolate Distortion Energy (kcal/mol)", y="Addition Barrier (kcal/mol)", data=df, hue="Ligand", legend=False, s=100)
    for line in range(1, df.shape[0] + 1):
        ax1.text(df["Thiolate Distortion Energy (kcal/mol)"][line]-0.03, df["Addition Barrier (kcal/mol)"][line],
                df["Ligand"][line], horizontalalignment='right', size='small', color='black')
    ax2 = fig.add_subplot(2, 2, 2)
    p2 = sns.scatterplot(x="Ligand Distortion Energy (kcal/mol)", y="Addition Barrier (kcal/mol)", data=df, hue="Ligand", legend=False, s=100)
    for line in range(1, df.shape[0] + 1):
        ax2.text(df["Ligand Distortion Energy (kcal/mol)"][line]+5, df["Addition Barrier (kcal/mol)"][line],
                df["Ligand"][line], horizontalalignment='left', size='small', color='black')
    ax3 = fig.add_subplot(2, 2, 3)
    p3 = sns.scatterplot(x="Activation Energy (kcal/mol)", y="Addition Barrier (kcal/mol)", data=df, hue="Ligand", legend=False, s=100)
    for line in range(1, df.shape[0] + 1):
        ax3.text(df["Activation Energy (kcal/mol)"][line]-0.5, df["Addition Barrier (kcal/mol)"][line],
                df["Ligand"][line], horizontalalignment='right', size='small', color='black')
    ax4 = fig.add_subplot(2, 2, 4)
    p4 = sns.scatterplot(x="Interaction Energy (kcal/mol)", y="Addition Barrier (kcal/mol)", data=df, hue="Ligand", legend=False, s=100)
    for line in range(1, df.shape[0] + 1):
        ax4.text(df["Interaction Energy (kcal/mol)"][line]+5, df["Addition Barrier (kcal/mol)"][line],
                df["Ligand"][line], horizontalalignment='left', size='small', color='black')
    # plt.legend(p4, loc="center right", borderaxespad=0.1, title="Ligands")
    plt.ylim(None, None); plt.xlim(None, None)
    # plt.tight_layout()
    plt.savefig("Distortion-Interaction Analysis")
    # plt.show()


