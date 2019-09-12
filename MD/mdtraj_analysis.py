'''
Python script for MD trajectory data analysis.
Written by Jonathan Yik Chang Ting (Student ID: 44254124) for UQ BAdSc(Hons) Honours Project 2019.
Project Title: Molecular Modelling of Reversible Covalent Inhibition of Bruton's Tyrosine Kinase by Cyanoacrylamides
Free to be distributed, edited and utilized for purposes beneficial to the world. :)
Though recognition of my effort is appreciated! XD
'''

import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

sns.set(context='paper', font_scale=1.5)
# sns.despine()  # Remove "chartjunk"

MD_PATH = "/home/s4425412/data/Honours/MD"
SINGLE_PLOT_WIDTH, SINGLE_PLOT_HEIGHT = 9, 7
SUBPLOTS_WIDTH, SUBPLOTS_HEIGHT = 10, 6


if __name__ == "__main__":

    dist_interest_thresh = 12

    # Base for elimination
    for i in [1, 2]:
        # Update DATA_PATH and reinitialise residue_dict for every replicate
        DATA_DIR = "{0}/{1}/MD/cov/Rep{2}/analysis/base".format(MD_PATH, "3", i)
        residue_dict = {"arg": {"filename_list": [], "interest_list": []},
                        "glu": {"filename_list": [], "interest_list": []},
                        "asp": {"filename_list": [], "interest_list": []},
                        "his": {"filename_list": [], "interest_list": []},
                        "lys": {"filename_list": [], "interest_list": []}}
        # Sort out all .dat files
        for j, filename in enumerate(os.listdir(DATA_DIR)):
            # Filter out files that are not .dat files
            if ".dat" not in filename:
                continue
            # assign all .dat files to respective residues dictionary
            for k, residue in enumerate(residue_dict.keys()):
                if residue in filename:
                    residue_dict[residue]["filename_list"].append(filename)
                else:
                    raise Exception("Residue not found!")
            with open(filename, "r+") as dist_data:
                for l, distance in enumerate(dist_data.readlines()):
                    if float(distance) <= dist_interest_thresh:
                        residue_dict[residue]["interest_list"].append(filename)
                        continue

        # Analyse each residue
        for j, residue in enumerate(residue_dict.keys()):
            combined_df = pd.DataFrame(list(range(0, 20001, 1)), columns=["Time (ps)"])  # Create empty dataframe
            for k, dat_file in enumerate(residue_dict[residue]["filename_list"]):
                resname = dat_file.split(".")[0].split("_")[-1].upper()
                df = pd.read_csv("{0}/{1}".format(DATA_DIR, dat_file), index_col=0, delim_whitespace=True)
                df.columns = [resname]
                combined_df = pd.concat([combined_df, df], axis=1, ignore_index=False)
            combined_df = combined_df.drop([0], axis=0)
            combined_df = combined_df.melt("Time (ps)", var_name="Residues", value_name=r"Distance ($\AA$)")
            fig = plt.figure(figsize=(SINGLE_PLOT_WIDTH, SINGLE_PLOT_HEIGHT))
            fig.subplots_adjust(top=0.9, bottom=0.1, left=0.1, right=0.9)
            ax = sns.lineplot(x="Time (ps)", y=r"Distance ($\AA$)", data=combined_df.head(100000), hue="Residues")
            leg = ax.legend(ncol=3, loc='best', fontsize='small')
            leg.set_title("")
            plt.ylim(None, None); plt.xlim(None, None)
            plt.tight_layout()
            plt.show()
            plt.savefig(r"Distance of C$_\alpha$ H from {0} Residues Replicate {1}".format(residue.upper(), i))

