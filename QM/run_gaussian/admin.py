"""
Python script to do administrative work on data directory
Python script to tabulate data obtained from Gaussian (version 16 at the moment)
Written by Jonathan Yik Chang Ting (Student ID: 44254124) for UQ BAdSc(Hons) Honours Project 2019.
Project Title: Molecular Modelling of Reversible Covalent Inhibition of Bruton’s Tyrosine Kinase by Cyanoacrylamides
Free to be distributed, edited and utilized for purposes beneficial to the world. :)
Though recognition of my effort is appreciated! XD
"""

import os
from os.path import isdir
from QM.run_gaussian.settings import DATA_PATH


def change_names(dir_path, new_str, old_str=None):
    file_list = [f for f in os.listdir(dir_path)]
    for j, title in enumerate(file_list):
        file_path = '{0}/{1}'.format(dir_path, title)
        if old_str:
            os.rename(file_path, file_path.replace(old_str, new_str))
        else:  # If no old_str specified, append
            os.rename(file_path, file_path + new_str)


def make_dir_for_files(dir_path):
    file_list = [f for f in os.listdir(dir_path) if not isdir('{0}/{1}'.format(dir_path, f)) and '.xyz' in f]
    for i, file in enumerate(file_list):
        new_dir_path = '{0}/{1}'.format(dir_path, rid_ext(file))
        os.mkdir(new_dir_path)


def group_files(dir_path):
    group_list = [g for g in os.listdir(dir_path) if isdir('{0}/{1}'.format(dir_path, g))]
    if len(group_list) == 0:
        make_dir_for_files(dir_path)
        group_list = [g for g in os.listdir(dir_path) if isdir('{0}/{1}'.format(dir_path, g))]
    file_list = [f for f in os.listdir(dir_path) if not isdir('{0}/{1}'.format(dir_path, f))]
    for i, file in enumerate(file_list):
        file_path = '{0}/{1}'.format(dir_path, file)
        for j, group in enumerate(group_list):
            if rid_ext(file) == group:
                old_str, new_str = '/{0}'.format(file), '/{0}/{1}'.format(group, file)
                os.rename(file_path, file_path.replace(old_str, new_str))


def rid_ext(file_name):
    if '.' in file_name:
        split_names = file_name.split('.')
        new_file_name = ''.join(split_names[:-1])
    else:
        new_file_name = file_name
    return new_file_name


def split_file(file_path, break_str, ext):
    break_index_list, prev_index = [], None
    with open(file_path, 'r') as f:
        line_list = f.readlines()
        for i, line in enumerate(line_list):
            if break_str in line:
                break_index_list.append(i)
        break_index_list.append(len(line_list)-1)  # Include the last conformer!
        for j, index in enumerate(break_index_list):
            if j == 0:
                prev_index = index
                continue
            title = rid_ext(line_list[prev_index + 1].split('\\')[-1].split('\n')[0])
            dir_path = '/'.join(file_path.split('/')[:-1])
            with open('{0}/{1}_c{2}.{3}'.format(dir_path, title, j, ext), 'w') as g:
                g.write(''.join(line_list[prev_index:index]))
            prev_index = index
    os.remove(file_path)


if __name__ == "__main__":
    operation_list = ['split', 'rename', 'group']
    input_dir = '{0}/Cross_Conjugation/Conformational_Search/mystery/Chloroform/Second_Products'.format(DATA_PATH)
    dircs = [f for f in os.listdir(input_dir) if isdir('{0}/{1}'.format(input_dir, f))]  # and 'Reactant' not in f
    operation_index = 2
    for i, dirc in enumerate(dircs):
        dir_path = '{0}/{1}'.format(input_dir, dirc)
        print("\n# Directory:\n", dir_path, "\n# Operation:", operation_list[operation_index])
        if operation_list[operation_index] == 'split':
            print("Splitting concatenated molecules (all.mol2) to individual conformers...")
            big_file_path = '{0}/{1}/all.mol2'.format(input_dir, dirc)
            if not os.path.exists(big_file_path):
                print("all.mol2 file not found!")
                continue
            split_file(big_file_path, '@<TRIPOS>MOLECULE', 'mol2')
        elif operation_list[operation_index] == 'rename':
            change_names(dir_path, '', 'mw')
        elif operation_list[operation_index] == 'group':
            print("Grouping conformers into individual dir...")
            group_files(dir_path)
        else:
            raise Exception("Unknown operation!")
    print("Operation success!")

