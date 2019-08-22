"""
Python script to generate Gaussian input files and job submission files on HPC clusters
Written by Jonathan Yik Chang Ting (Student ID: 44254124) for UQ BAdSc(Hons) Honours Project 2019.
Project Title: Molecular Modelling of Reversible Covalent Inhibition of Bruton’s Tyrosine Kinase by Cyanoacrylamides
Free to be distributed, edited and utilized for purposes beneficial to the world. :)
Though recognition of my effort is appreciated! XD
"""

import sys
import os
from os.path import isdir
from run_gaussian.settings import theory_lvl_list, keyword_dict, DATA_PATH

NCPUS, SOFTWARE, VERSION, CLUSTER = int(8), 'g16', 'b01', 'Raijin'  # Unlikely to change
WALLTIME, VMEM, JOBFS = '12:00:00', int(70000), int(3000)  # Common variables
MEM, CHARGE, SPIN = VMEM, int(0), int(1)
METHOD, BASIS_SET, CALC_TYPE, FREQ, SCRF = 'M062X', '6-31+G(d)', ' opt', ' freq=noraman', ' scrf=(cpcm,solvent=water)'


def gen_from_dir(input_dir, run_calc_type, cluster=CLUSTER, combination=None, solvent='water'):
    assert run_calc_type in keyword_dict.keys(), 'Calculation type not known!'
    groups = [f for f in os.listdir(input_dir) if isdir('{0}/{1}'.format(input_dir, f))]  # and 'Reactant' not in f
    for i, group in enumerate(groups):  # Potentially R, TSS, P
        group_dir = '{0}/{1}'.format(input_dir, group)  # Absolute path to the directories
        conformers = [g for g in os.listdir(group_dir) if isdir('{0}/{1}'.format(group_dir, g))]  # and 'TR' not in g
        for j, title in enumerate(conformers):
            conformer_dir = '{0}/{1}'.format(group_dir, title)
            edit_charge = True if run_calc_type == 'FBRGO' or run_calc_type == 'TSGOVF' or 'TS' in title or '-' in title else False
            SCRF = ' scrf=(cpcm,solvent={0})'.format(solvent)
            gauss_inp_geom_opt(title, conformer_dir, keyword_dict[run_calc_type]['mem'], NCPUS, combination, keyword_dict[run_calc_type]['freq'],
                               SCRF, keyword_dict[run_calc_type]['type'], CHARGE, SPIN, edit_charge)
            write_job_script(title, conformer_dir, cluster, NCPUS, SOFTWARE, keyword_dict[run_calc_type]['time'], VMEM, keyword_dict[run_calc_type]['jobfs'], VERSION)


def gauss_inp_geom_opt(title, dirc, mem=MEM, ncpus=NCPUS, combination='{0}/{1}'.format(METHOD, BASIS_SET), freq=FREQ,
                       scrf=SCRF, calc_type=CALC_TYPE, charge=CHARGE, spin=SPIN, edit_charge=False):
    charge += keyword_dict['FBRGO']['edit_charge'] if edit_charge else int(0)
    with open('{0}/{1}.inp'.format(dirc, title), 'w') as f:
        f.write('%mem={0}mb\n%nproc={1}\n%chk={2}.chk\n'.format(mem, ncpus, title))
        f.write('# {0}{1}{2}{3}'.format(combination, scrf, calc_type, freq))
        f.write('\n\n{0}\n\n{1} {2}\n'.format(title, charge, spin))
        with open('{0}/{1}.xyz'.format(dirc, title), 'r') as g:
            line_list = g.readlines()
            has_energy_name = 'Energy' in line_list[1] or 'w' in line_list[1] or 'Name' in line_list[1]  # Potential bug-breeding spots
            has_path = ':' in line_list[1]
            for i, line in enumerate(line_list):
                if has_energy_name and i > 1:
                    f.write('{0}'.format(line))
                elif has_path and i > 1:
                    f.write('\n{0}'.format(line)) if i == 2 else f.write('{0}'.format(line))
                elif not(has_energy_name) and not(has_path) and i > 0:
                    f.write('{0}'.format(line))
                else:
                    pass
            f.write('\n')  # Necessary


def gauss_inp_freeze_opt_gzmat(title, dirc, mem=MEM, ncpus=NCPUS, method=METHOD, basis_set=BASIS_SET, freq=FREQ,
                               scrf=SCRF, calc_type=CALC_TYPE, charge=CHARGE, spin=SPIN, edit_charge=True):
    charge_spin_line = int(4)  # Potential bug-breeding spot
    charge += keyword_dict['FBRGO']['edit_charge'] if edit_charge else int(0)
    target_bond, bond_num = find_bond(dirc, title, 'S')
    with open('{0}/{1}.inp'.format(dirc, title), 'w') as f:
        f.write('%mem={0}mb\n%nproc={1}\n%chk={2}.chk\n'.format(mem, ncpus, title))
        f.write('# {0}/{1}{2}{3}{4}'.format(method, basis_set, scrf, calc_type, freq))
        f.write('\n\n{0}\n\n{1} {2}\n'.format(title, charge, spin))
        with open('{0}/{1}.gzmat'.format(dirc, title), 'r') as g:
            line_list = g.readlines()
            for i, line in enumerate(line_list):
                if i > charge_spin_line:
                    if '=' in line:
                        line = line.replace('=', '')
                    if 'Variables' in line:
                        line = '\n'
                    elif target_bond in line and i > bond_num:
                        line = line.replace('\n', ' f\n')
                    f.write('{0}'.format(line))
            f.write('\n')


def find_bond(dirc, title, element):
    target_bond, bond_num = None, None
    with open('{0}/{1}.gzmat'.format(dirc, title), 'r') as f:
        line_list = f.readlines()
        for i, line in enumerate(line_list):
            if element in line:
                target_bond, bond_num = line.split('  ')[2], i
    if not target_bond or not bond_num:
        raise Exception('Element not found!')
    return target_bond, bond_num


def write_job_script(title, dirc, cluster=CLUSTER, ncpus=NCPUS, software=SOFTWARE, walltime=WALLTIME, vmem=VMEM, jobfs=JOBFS, version=VERSION):
    with open('{0}/{1}.sh'.format(dirc, title), 'w') as f:
        if cluster == 'Raijin':
            f.write('#!/bin/bash\n#PBS -l wd\n#PBS -q normal\n')
            f.write('#PBS -l walltime={0},mem={1}mb,ncpus={2},software={3},jobfs={4}mb'.format(walltime, vmem, ncpus, software, jobfs))
            f.write('\n\nmodule load gaussian/{0}{1}'.format(software, version))
            f.write('\n{0} < {1}.inp > {1}.out 2>&1'.format(software, title))
        elif cluster == 'RCC':
            f.write('#!/bin/bash\n#PBS -S /bin/bash\n#PBS -l walltime={0}\n#PBS -A UQ-SCI-SCMB\n'.format(walltime))
            f.write('#PBS -l select=1:ncpus={0}:mem={1}MB'.format(ncpus, vmem))
            f.write('\n\ncd $PBS_O_WORKDIR')
            f.write('\n\nmodule load gaussian/{0}-{1}-bash'.format(software, version.upper()))
            # f.write('\nulimit -c 0; ulimit -d hard; ulimit -f hard; ulimit -l hard; ulimit -m hard; ulimit -n hard; '
            #         'ulimit -s hard; ulimit -t hard; ulimit -u hard;')
            f.write('\n{0} < {1}.inp > {1}.out'.format(software, title))
        else:
            raise Exception('Cluster not recognized!')


def SCS_corr(line_list):  # Spin-component scaled correction (to MP2 methods only)
    param_start_line, E_SCF, E_SCS_MP2 = "Spin components of T(2) and E(2)", None, None
    for i, line in enumerate(line_list):
        if 'HF' in line or '\\MP2' in line:
            value_inc = line.split('HF')[-1].strip() if 'HF' in line else line.split('MP2')[0].strip()
            if '\\' in value_inc:
                E_SCF = value_inc.split('\\')[0].split('=')[-1]
            else:
                E_SCF = value_inc + line_list[i-1].split('\\')[0]
                E_SCF = E_SCF.split('=')[-1]
            try:
                E_SCF = E_SCF.replace(' ', '')
                E_SCF = float(E_SCF);
            except ValueError:
                print("E_SCF not number!")
                continue
        if param_start_line in line:  # Found the line where parameters for SCS can be located
            aaE2, abE2, bbE2 = float(line_list[i-1].split('     ')[-1].strip().replace('D', 'e')), float(line_list[i-2].split('     ')[-1].strip().replace('D', 'e')), float(line_list[i-3].split('     ')[-1].strip().replace('D', 'e'))
            E_2_S, E_2_T = abE2, aaE2 + bbE2
            E_2_SCS = 6/5*E_2_S + E_2_T/3
            E_SCS_MP2 = E_SCF + E_2_SCS
            print("E_2_SCS:", E_2_SCS, "; E_SCF", E_SCF)
            break
    return E_SCS_MP2


if __name__ == "__main__":
    if len(sys.argv) == 2:
        input_dir = os.path.abspath(sys.argv[1])
        cluster = 'Raijin'
    elif len(sys.argv) > 2:
        input_dir = os.path.abspath(sys.argv[1])
        cluster = sys.argv[2]
    else:
        input_dir = '{0}/QM_Calculations/Benchmarking'.format(DATA_PATH)
    job_list = ['geom_opt', 'spe', 'SCS']
    job = job_list[1]
    if job == 'geom_opt':
        input_dir = '{0}/Cross_Conjugation/Conformational_Search/alpha_santonin/Water/Reactants'.format(DATA_PATH)
        combination = '{0}/{1}'.format(METHOD, BASIS_SET)
        gen_from_dir(input_dir, run_calc_type='GOVF', cluster='RCC', combination=combination, solvent='water')
    elif job == 'spe':
        input_dir = '{0}/Cross_Conjugation/Single_Point_Energy'.format(DATA_PATH)
        for i, combination in enumerate(theory_lvl_list):
            # theory_dir = '{0}/Combination{1}/Diethyl_ether'.format(input_dir, i + 1)
            # theory_dir = '{0}/Combination8/Water'.format(input_dir)
            theory_dir = '{0}/alpha_santonin/Water'.format(input_dir)
            gen_from_dir(theory_dir, run_calc_type='SPEiS', cluster='RCC', combination=combination, solvent='water')
    elif job == 'SCS':
        input_dir = '{0}/QM_Calculations/Benchmarking/Combination3/Diethyl_ether'.format(DATA_PATH)
        groups = [f for f in os.listdir(input_dir) if
                  isdir('{0}/{1}'.format(input_dir, f)) and 'exclude' not in f]  # and 'Reactant' not in f
        print("\n# Correcting SPE from MP2 methods by scaling spin components...")
        print("\n# Input directory:\n", input_dir, "\n\n# Groups:\n", groups, "\n")
        for i, group in enumerate(groups):  # Potentially R, TSS, PE_SCF
            group_dir = '{0}/{1}'.format(input_dir, group)
            conformers = [g for g in os.listdir(group_dir) if isdir('{0}/{1}'.format(group_dir, g))]
            for j, title in enumerate(conformers):
                conformer_dir = '{0}/{1}'.format(group_dir, title)
                print("Conformer:", title)
                with open('{0}/{1}.out'.format(conformer_dir, title), 'r') as f:
                    line_list = f.readlines()
                    line_list.reverse()
                    E_SCS_MP2 = SCS_corr(line_list)
                    print("SCS-MP2 Potential Energy:", E_SCS_MP2)
    else:
        raise Exception("Job specified not identified!")


