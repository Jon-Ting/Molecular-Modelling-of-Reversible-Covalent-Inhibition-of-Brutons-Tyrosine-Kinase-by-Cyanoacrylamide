# Molecular Modelling of Covalent Inhibition of Bruton’s Tyrosine Kinase by Cyanoacrylamides

## Project Details
* Author: Jonathan Yik Chang Ting
* Student ID: 44254124
* Supervisors: Associate Professor Elizabeth Krenske, Professor Alan Mark, Dr Martin Stroet
* Institution: University of Queensland
* Organisation: School of Chemistry and Molecular Biosciences
* Programme: Bachelor of Advanced Science (Honours)
* Course: Honours in Chemistry (CHEM6511)
* Year: 2019

## Project Phases 

### Quantum Mechanical Calculations
* QM directory contains the Python files written for various purposes related to the QM calculation phase of the project:
    > Files in run_gaussian directory are for the automation of generation of QM calculation input files and job submission files, grouping of files with same names, and QM data tabulation.
    > Files in Csearch directory are some functions I used while trying to come up with ways to explain the conformaitonal searching phase of the project for my Honours seminar.
    > Files in visualG directory are either for plotting of energy profiles for Cross-Conjugation project (a side project) or the linear regression analysis of certain molecular properties. The data are stored in plot_config.py while the plotting codes are mainly written in plot_fig.py. Figures generated from the analysis on data generated using Combination G and I were stored in the respective directories.
    > calculation.py was used for interconversion between kinetics and thermodynamics parameters during the early phase of the project.
* QM_Gaussian directory contains the same files as QM/run_gaussian. It is generalised for the usage of others who wish to automate some tasks involved in the conformational searching phase of the QM calculations of flexible molecules using Gaussian.

 ### Molecular Dynamics Simulations
* MD directory contains the Python files written for the analysis of simulation trajectories of BTK inhibited by different inhibitors.

