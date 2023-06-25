This directory contains several Python files and Bash file that I have written during my Honours project for the automation of certain tasks.

Some descriptions on the files:
settings.py - Contains variables for different type of Gaussian jobs that is used by prepGaussian.py.
admin.py - Group files with the same names (before extension, e.g. abc.inp is the same as abc.xyz) into individual directories.
prepGaussian.py - Batch generation of Gaussian input files and job submission files on HPCs with PBS Scheduler.
tabulate.py - Batch tabulation of interested values from Gaussian output files into an Excel file.
gsub.sh - Batch submission of Gaussian QM calculation jobs on HPCs.

Typical work flow for a mechanism-based project where flexible molecules are involved would be:
1) Conduct conformational searches on the species along the reaction coordinate (using MacroModel).
2) Export all conformers within 3 kcal/mol of the lowest energy structure to a directory in .xyz format. The naming convention is very important, be consistent, but make sure each conformer has a different name (I do this by adding numbers at the end of their names, signifying their ranks from the conformational searches).
3) Change your input directory in admin.py to where you store the coordinate files and run it. This will group all of the conformers into individual directories.
4) Change your input directory in prepGaussian.py to the same location and run it. This will generate the Gaussian input files and job submission files in the corresponding directory. Make sure you check at least a few of the files generated to see that you got the charges, spacing at the end of file, solvents, resources requested, etc right. I named all of the Gaussian job submission jobs with *.sh, feel free to change it according to your preference (Line 52 in prepGaussian.py).
5) Copy the directories across to the HPCs (Raijin/Gadi/Tinaroo/Awoonga, NOTE: The details for job submissions on Raijin and RCC HPCs are different, specify the cluster before running prepGaussian.py in Step 4).
6) Change directory to the directory that contains all the conformers subdirectories and run 'gsub.sh' (Make sure it's an executable, if you can't run it, use chmod to change it). This will submit all of the Gaussian submission jobs to the HPC. Note that prior to this you need to adjust Line 6 in the gsub.sh file to the naming convention you give to your Gaussian submit files if you have changed them in the prepGaussian.py.
7) After the jobs are done, copy them over to your local machine.
8) Change your input directory in tabulate.py to the directory that contains all the conformers and run it. Note that you need to have the Python package pandas installed for it to work.

Detailed description could be found in the files.

