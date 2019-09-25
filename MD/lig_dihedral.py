'''
Python script for MD trajectory data analysis.
Written by Jonathan Yik Chang Ting (Student ID: 44254124) for UQ BAdSc(Hons) Honours Project 2019.
Project Title: Molecular Modelling of Reversible Covalent Inhibition of Bruton's Tyrosine Kinase by Cyanoacrylamides
Free to be distributed, edited and utilized for purposes beneficial to the world. :)
Though recognition of my effort is appreciated! XD
'''


import os
import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from Honours.MD.config import MD_PATH, SINGLE_PLOT_WIDTH, SINGLE_PLOT_LEG_HEIGHT,SINGLE_PLOT_LEG_WIDTH, \
    SINGLE_PLOT_LEG_HEIGHT, FOUR_PLOTS_WIDTH, FOUR_PLOTS_HEIGHT, SIX_PLOTS_WIDTH, SIX_PLOTS_HEIGHT, inhibitor_list

sns.set(context='paper', font_scale=1.5)
CURR_DIR = os.getcwd()


if __name__ == "__main__":

    sort_dict, system_type = True, "noncov"

    if sort_dict:  # Sorting and storing the dihedral data
        dih_data_dict = {"1": {"Rep1A": None, "Rep1B": None, "Rep2A": None, "Rep2B": None},
                         "3": {"Rep1A": None, "Rep1B": None, "Rep2A": None, "Rep2B": None},
                         "4": {"Rep1A": None, "Rep1B": None, "Rep2A": None, "Rep2B": None},
                         "5": {"Rep1A": None, "Rep1B": None, "Rep2A": None, "Rep2B": None},
                         "7": {"Rep1A": None, "Rep1B": None, "Rep2A": None, "Rep2B": None},
                         "9": {"Rep1A": None, "Rep1B": None, "Rep2A": None, "Rep2B": None}}

        for i, inhibitor in enumerate(inhibitor_list):
            for j in [1, 2]:
                DATA_DIR = "{0}/{1}/MD/{2}/Rep{3}/analysis/dih".format(MD_PATH, str(inhibitor), system_type, j)
                SYSTEM_DIR = "{0}/{1}/{2}".format(CURR_DIR, str(inhibitor), system_type)

                # Sort out all .dat files
                for j, filename in enumerate(os.listdir(DATA_DIR)):
                    # assign all .dat files to respective residues dictionary
                    for k, residue in enumerate(residue_dict.keys()):
                        if residue in filename:
                            if not os.path.exists("{0}/residue_dict.txt".format(SYSTEM_DIR)):
                                with open("{0}/{1}".format(SYSTEM_DIR, "{0}_list.txt".format(residue)), "a+") as res_list_file:
                                    res_list_file.write("{0}\n".format(filename))
                                    residue_dict[residue].append(filename)

                json.dump(within_dist_dict, open("{0}/within_dist_dict.txt".format(STORE_DIR), "w"), indent=4)



    for i, inhibitor in enumerate(inhibitor_list):

        # Analyse each residue
        within_dist_dict = json.load(open("{0}/within_dist_dict.txt".format(STORE_DIR)))  # Residues of interest only
        residue_dict = json.load(open("{0}/residue_dict.txt".format(SYSTEM_DIR)))  # All residues

        chosen_dict = within_dist_dict  # Choose the dictionary with right data ***************************************

        # Charged Residues from CaH
        PLOTS_WIDTH, PLOTS_HEIGHT = SIX_PLOTS_WIDTH, SIX_PLOTS_HEIGHT
        num_row, subplot_index = int(num_subplots / 2), 0
        fig, axes = plt.subplots(nrows=num_row, ncols=2, sharex=True, sharey=True, figsize=(PLOTS_WIDTH, PLOTS_HEIGHT))
        fig.subplots_adjust(top=0.95, bottom=0.07, left=0.08, right=0.97, wspace=0.15, hspace=0.15)
        fig.suptitle(r"Distance of Charged Residues of Interest from C$_\alpha$ H", horizontalalignment='center',
                     fontsize=14, weight='bold')
        fig.text(0.5, 0.02, "Time (ns)", va='center', ha='center')
        fig.text(0.02, 0.5, r"Distance ($\AA$)", va='center', ha='center', rotation='vertical')
        for j, residue in enumerate(chosen_dict.keys()):
            label_list = []
            combined_df = pd.DataFrame(list(np.arange(0.000, 100.000, 0.005)), columns=["Time (ns)"])
            for k, dat_file in enumerate(chosen_dict[residue]):
                resname = dat_file.split(".")[0].split("_")[-1].upper()
                df = pd.read_csv("{0}/{1}".format(DATA_DIR, dat_file), index_col=0, delim_whitespace=True)
                df.columns = [resname]
                combined_df = pd.concat([combined_df, df], axis=1, ignore_index=False)
                label_list.append(resname)
                '''
                combined_df = combined_df.melt("Time (ps)", var_name="Residues", value_name=r"Distance ($\AA$)")
                # combined_df = combined_df.drop([0], axis=0)
                fig, ax = plt.subplots()
                fig.subplots_adjust(top=0.9, bottom=0.1, left=0.1, right=0.9)
                sns.lineplot(x="Time (ps)", y=r"Distance ($\AA$)", data=combined_df, legend="full")
                sns.despine(ax=ax, left=True)  # Remove "chartjunk"
                leg = ax.legend(ncol=3, loc='best', fontsize='small')
                leg.set_title("Residues")
                label_list = combined_df["Residues"]
                for k, res_label in zip(leg.texts, label_list):
                    k.set_text(res_label)
                '''  # Seaborn plots (unsuccessful)
            if combined_df.shape[1] <= 1:  # Don't plot if only time column exists
                continue
            subplot_index += 1
            ax1 = plt.subplot(int("{0}2{1}".format(num_row, subplot_index)))
            # plt.figure(figsize=(SINGLE_PLOT_WIDTH, SINGLE_PLOT_LEG_HEIGHT))
            combined_df.plot(x="Time (ns)", y=label_list, kind="line", ax=plt.gca(), lw=0.3, linestyle=":")
            ax1.tick_params(top=False, bottom=True, left=True, right=False, labelleft=True, labelbottom=True)
            ax1.xaxis.label.set_visible(False)
            leg = plt.legend(loc="upper right")
            # leg = plt.legend(loc="upper center", ncol=3, bbox_to_anchor=(0.5, -0.3))
            plt.setp(leg.get_lines(), linewidth=4)
            # plt.ylim(None, None); plt.xlim(None, None)
        fig.delaxes(axes[-1][-1])
        plt.savefig(r"{0}/Distance of CaH from Charged Residues Rep {1} within {2} A".format(STORE_DIR, i, dist_thresh))
        # plt.show()
