'''
Python script for MD trajectory data analysis.
Written by Jonathan Yik Chang Ting (Student ID: 44254124) for UQ BAdSc(Hons) Honours Project 2019.
Project Title: Molecular Modelling of Reversible Covalent Inhibition of Bruton's Tyrosine Kinase by Cyanoacrylamides
Free to be distributed, edited and utilized for purposes beneficial to the world. :)
Though recognition of my effort is appreciated! XD
'''

import os
import json
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

sns.set(context='paper', font_scale=1.5)

MD_PATH = "/home/s4425412/data/Honours/MD"
SINGLE_PLOT_WIDTH, SINGLE_PLOT_HEIGHT = 9, 7
SUBPLOTS_WIDTH, SUBPLOTS_HEIGHT = 10, 6


if __name__ == "__main__":

    dist_interest_thresh = 12

    # Base for elimination
    for i in [1, 2]:

        # Update DATA_PATH and reinitialise residue_dict for every replicate
        DATA_DIR = "{0}/{1}/MD/old_cov/Rep{2}/analysis/base".format(MD_PATH, "3", i)
        '''
        residue_dict = {"arg": {"filename_list": [], "interest_list": []},
                        "glu": {"filename_list": [], "interest_list": []},
                        "asp": {"filename_list": [], "interest_list": []},
                        "his": {"filename_list": [], "interest_list": []},
                        "lys": {"filename_list": [], "interest_list": []}}
        # Empty the .dat files
        for j, residue in enumerate(residue_dict.keys()):
            open("{0}_list.txt".format(residue), "w").close()
        # Sort out all .dat files
        for j, filename in enumerate(os.listdir(DATA_DIR)):
            # Filter out files that are not .dat files
            if ".dat" not in filename:
                continue
            # assign all .dat files to respective residues dictionary
            for k, residue in enumerate(residue_dict.keys()):
                if residue in filename:
                    with open("{0}/{1}".format(os.getcwd(), "{0}_list.txt".format(residue)), "a+") as res_list_file:
                        res_list_file.write("{0}\n".format(filename))
                        residue_dict[residue]["filename_list"].append(filename)
                        with open("{0}/{1}".format(DATA_DIR, filename), "r") as dist_data:
                            with open("{0}/{1}".format(os.getcwd(), "close_res_list.txt"), "a") as close_res_file:
                                for l, line_entry in enumerate(dist_data.readlines()):
                                    distance = line_entry.split()[-1]
                                    try:
                                        if float(distance) <= dist_interest_thresh:
                                            residue_dict[residue]["interest_list"].append(filename)
                                            close_res_file.write("{0}\n".format(filename))
                                            break
                                    except ValueError:
                                        continue
        json.dump(residue_dict, open("residue_dict.txt", "w"))
        '''

        # Analyse each residue
        residue_dict = json.load(open("residue_dict.txt"))
        for j, residue in enumerate(residue_dict.keys()):
            combined_df = pd.DataFrame(list(range(0, 20001, 1)), columns=["Time (ps)"])  # Create empty dataframe
            for k, dat_file in enumerate(residue_dict[residue]["filename_list"]):
                resname = dat_file.split(".")[0].split("_")[-1].upper()
                df = pd.read_csv("{0}/{1}".format(DATA_DIR, dat_file), index_col=0, delim_whitespace=True)
                df.columns = [resname]
                combined_df = pd.concat([combined_df, df], axis=1, ignore_index=False)
            # combined_df = combined_df.drop([0], axis=0)
            combined_df = combined_df.melt("Time (ps)", var_name="Residues", value_name=r"Distance ($\AA$)")
            fig, ax = plt.subplots()
            # plt.figure(figsize=(SINGLE_PLOT_WIDTH, SINGLE_PLOT_HEIGHT))
            fig.subplots_adjust(top=0.9, bottom=0.1, left=0.1, right=0.9)
            g = sns.lineplot(x="Time (ps)", y=r"Distance ($\AA$)", data=combined_df.head(30000), legend=False)
            # sns.despine(ax=ax, left=True)  # Remove "chartjunk"
            leg = ax.legend(ncol=3, loc='best', fontsize='small')
            leg.set_title("Residues")
            label_list = combined_df["Residues"]
            for k, res_label in zip(leg.texts, label_list):
                k.set_text(res_label)
            plt.ylim(None, None); plt.xlim(None, None)
            plt.tight_layout()
            plt.savefig("{0}/Distance of CaH from {1} Residues Replicate {2}".format(os.getcwd(), residue.upper(), i))
            plt.show()
