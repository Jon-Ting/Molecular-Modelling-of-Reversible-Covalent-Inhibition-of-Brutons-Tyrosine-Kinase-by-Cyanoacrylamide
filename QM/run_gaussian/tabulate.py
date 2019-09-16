"""
Python script to tabulate data obtained from Gaussian (version 16 at the moment)
Written by Jonathan Yik Chang Ting (Student ID: 44254124) for UQ BAdSc(Hons) Honours Project 2019.
Project Title: Molecular Modelling of Reversible Covalent Inhibition of Bruton’s Tyrosine Kinase by Cyanoacrylamides
Free to be distributed, edited and utilized for purposes beneficial to the world. :)
Though recognition of my effort is appreciated! XD
"""

import os
import pandas as pd
import re
from os.path import isdir
from QM.run_gaussian.settings import DATA_PATH


def sort_human(list):
    convert = lambda text: float(text) if text.isdigit() else text
    alphanum = lambda key: [convert(c) for c in re.split('([-+]?[0-9]*\.?[0-9]*)', key)]
    list.sort(key=alphanum)
    return list


def replaceMultiple(main_str, toBeReplaced, new_str):
    for elem in toBeReplaced:
        if elem in main_str:
            main_str = main_str.replace(elem, new_str)
    return main_str


def find_val(line_list, target_str):
    val, isEnergy, isMethod = None, 'Energies' in target_str[0] or 'Enthalpies' in target_str[0], '%chk' in target_str[0]
    for i, string in enumerate(target_str):
            for j, line in enumerate(line_list):
                if string in line:
                    if isMethod:
                        value_inc = line_list[j - 2]
                        val = value_inc.split(' ')[2]; break
                    else:
                        value_inc = line.split(string)[-1].strip()
                        if isEnergy:
                            val = float(value_inc.split('=')[-1].strip().replace(' ', '')); break
                        else:
                            if '\\' in value_inc:
                                val = value_inc.split('\\')[0].split('=')[-1]
                            else:
                                val = value_inc + line_list[j - 1].split('\\')[0]
                                val = val.split('=')[-1]
                            try:
                                val = val.replace(' ', '')
                                float(val); break
                            except ValueError:
                                continue
    if val is None:
        raise Exception('Target string {0} not found!'.format(target_str))
    # print(target_str, val, line_list[-3])  # For debugging
    return val


if __name__ == "__main__":
    class_type = 'Second_Products'
    input_dir = '{0}/Cross_Conjugation/Conformational_Search/mystery/Chloroform/{1}'.format(DATA_PATH, class_type)
    groups = [f for f in os.listdir(input_dir) if isdir('{0}/{1}'.format(input_dir, f)) and 'exclude' not in f]  # and 'Reactant' not in f
    print("\n# Tabulating values of interest from Gaussian .out files to an Excel sheet...")
    print("\n# Input directory:\n", input_dir, "\n\n# Groups:\n", groups, "\n")
    method_list, title_list, molecule_list, conformer_list, NImag_list, Z_list, E_list, H_list, G_list = [], [], [], [], [], [], [], [], []
    var_fill_list = [method_list, NImag_list, Z_list, E_list, H_list, G_list]
    keyword_list = [['%chk'], ['NI', 'Im'], ['zero-point Energies'], ['HF'], ['thermal Enthalpies'], ['thermal Free Energies']]
    for i, group in enumerate(groups):  # Potentially R, TSS, P
        group_dir = '{0}/{1}'.format(input_dir, group)
        conformers = [g for g in os.listdir(group_dir) if isdir('{0}/{1}'.format(group_dir, g))]
        for j, title in enumerate(conformers):
            conformer_dir = '{0}/{1}'.format(group_dir, title)
            print("Conformer:", title)
            with open('{0}/{1}.out'.format(conformer_dir, title), 'r') as f:
                line_list = f.readlines()
                line_list.reverse()
                for i, var_list in enumerate(var_fill_list):
                    val = find_val(line_list, keyword_list[i])
                    var_list.append(val)
                title_list.append(title)
                molecule_list.append(replaceMultiple(title.split('c')[0].replace('_', ''), ['TR', 'TSS', 'TP'], ''))
                conformer_list.append('c' + title.split('c')[-1])
    data = {'Method': method_list, 'Title': title_list, 'Molecule': molecule_list, 'Conformer': conformer_list, 'NImag': NImag_list,
            'Z (Hartree)': Z_list, 'E (Hartree)': E_list, 'H (Hartree)': H_list, 'G (Hartree)': G_list}
    df = pd.DataFrame(data, columns=['Method', 'Title', 'Molecule', 'Conformer', 'NImag', 'Z (Hartree)', 'E (Hartree)', 'H (Hartree)', 'G (Hartree)'])
    title_list = sort_human(title_list)
    sorted_df = df.set_index('Title').reindex(title_list).reset_index()
    print("# Sorted data frame:\n", sorted_df)

    print("# Writing to Excel sheet...")
    # sorted_df.to_excel('{0}/Energies.xlsx'.format(input_dir), sheet_name='Sheet1')
    writer = pd.ExcelWriter('{0}/Energies_{1}.xlsx'.format(input_dir, class_type), engine='xlsxwriter')
    sorted_df.to_excel(writer, startrow=1, sheet_name='Sheet1', index=False)
    workbook = writer.book
    worksheet = writer.sheets['Sheet1']
    for i, col in enumerate(sorted_df.columns):
        column_len = sorted_df[col].astype(str).str.len().max()
        column_len = max(column_len, len(col)) + 2
        worksheet.set_column(i, i, column_len)
    writer.save()

