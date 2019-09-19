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
CURR_DIR = os.getcwd()
SINGLE_PLOT_WIDTH, SINGLE_PLOT_HEIGHT = 9, 7
SINGLE_PLOT_LEG_WIDTH, SINGLE_PLOT_LEG_HEIGHT = 9, 10
SUBPLOTS_WIDTH, SUBPLOTS_HEIGHT = 10, 6


if __name__ == "__main__":

    sort_dict, dist_thresh, inhibitor, system_type = True, 10, 3, "old_cov"

    # Base for elimination
    for i in [1, 2]:

        # Update DATA_PATH and reinitialise residue_dict for every replicate
        DATA_DIR = "{0}/{1}/MD/{2}/Rep{3}/analysis/base".format(MD_PATH, str(inhibitor), system_type, i)
        SYSTEM_DIR = "{0}/{1}/{2}".format(CURR_DIR, str(inhibitor), system_type)
        STORE_DIR = "{0}/within_{1}A".format(SYSTEM_DIR, dist_thresh)

        if sort_dict:
            # Sorting and storing the residue data
            within_dist_dict = {"arg": [], "glu": [], "asp": [], "his": [], "lys": []}
            if not os.path.exists("{0}/residue_dict.txt".format(SYSTEM_DIR)):
                residue_dict = {"arg": [], "glu": [], "asp": [], "his": [], "lys": [], "ligand": []}
            else:
                residue_dict = json.load(open("{0}/residue_dict.txt".format(SYSTEM_DIR)))  # All residues
                # Empty the .dat files
                for j, residue in enumerate(residue_dict.keys()):
                    open("{0}/{1}_list.txt".format(SYSTEM_DIR, residue), "w").close()
            # Sort out all .dat files
            for j, filename in enumerate(os.listdir(DATA_DIR)):
                # Filter out files that are not .dat files
                if ".dat" not in filename: continue
                # assign all .dat files to respective residues dictionary
                for k, residue in enumerate(residue_dict.keys()):
                    if residue in filename:
                        if not os.path.exists("{0}/residue_dict.txt".format(SYSTEM_DIR)):
                            with open("{0}/{1}".format(SYSTEM_DIR, "{0}_list.txt".format(residue)), "a+") as res_list_file:
                                res_list_file.write("{0}\n".format(filename))
                                residue_dict[residue].append(filename)
                        with open("{0}/{1}".format(DATA_DIR, filename), "r") as dist_data:

                            # Create necessary directories and files if they don't exist already
                            if not os.path.exists(STORE_DIR):
                                os.mkdir(STORE_DIR)
                            if not os.path.exists("{0}/close_res_list.txt".format(STORE_DIR)):
                                open("{0}/close_res_list.txt".format(STORE_DIR), "w").close()

                            with open("{0}/{1}".format(STORE_DIR, "close_res_list.txt"), "a") as close_res_file:
                                for l, line_entry in enumerate(dist_data.readlines()):
                                    distance = line_entry.split()[-1]
                                    try:
                                        if float(distance) <= dist_thresh:
                                            within_dist_dict[residue].append(filename)
                                            close_res_file.write("{0}\n".format(filename))
                                            break
                                    except (ValueError, TypeError) as error:  # If distance is NaN or dist_thresh=="All"
                                        continue
                    elif residue == "ligand":
                        try:
                            resname = filename.split(".")[0].split("_")[-1].upper()
                            int(resname)
                        except TypeError:
                            continue
                    else:
                        continue
            json.dump(within_dist_dict, open("{0}/within_dist_dict.txt".format(STORE_DIR), "w"), indent=4)
            if not os.path.exists("{0}/residue_dict.txt".format(SYSTEM_DIR)):
                json.dump(residue_dict, open("{0}/residue_dict.txt".format(SYSTEM_DIR), "w"), indent=4)

        # Analyse each residue
        within_dist_dict = json.load(open("{0}/within_dist_dict.txt".format(STORE_DIR)))  # Residues of interest only
        residue_dict = json.load(open("{0}/residue_dict.txt".format(SYSTEM_DIR)))  # All residues

        chosen_dict = within_dist_dict  # Choose the dictionary with right data

        for j, residue in enumerate(chosen_dict.keys()):
            label_list = []
            combined_df = pd.DataFrame(list(range(0, 20001, 1)), columns=["Time (ps)"])  # Create empty dataframe
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
            plt.figure(figsize=(SINGLE_PLOT_WIDTH, SINGLE_PLOT_LEG_HEIGHT))
            combined_df.plot(x="Time (ps)", y=label_list, kind="line", ax=plt.gca(), lw=0.3, linestyle=":")
            plt.xlabel("Time (ps)")
            plt.ylabel(r"Distance ($\AA$)")
            leg = plt.legend(loc="upper center", ncol=3, bbox_to_anchor=(0.5, -0.12))
            plt.setp(leg.get_lines(), linewidth=4)
            plt.ylim(None, None); plt.xlim(None, None)
            plt.tight_layout()
            plt.savefig(r"{0}/Distance of CaH from {1} Residues Replicate {2} within {3} A".format(STORE_DIR,
                        residue.upper(), i, dist_thresh))
            # plt.show()
