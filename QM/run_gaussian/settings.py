"""
Python script to store the setting variables for other scripts
Written by Jonathan Yik Chang Ting (Student ID: 44254124) for UQ BAdSc(Hons) Honours Project 2019.
Project Title: Molecular Modelling of Reversible Covalent Inhibition of Bruton’s Tyrosine Kinase by Cyanoacrylamides
Free to be distributed, edited and utilized for purposes beneficial to the world. :)
Though recognition of my effort is appreciated! XD
"""

DATA_PATH = 'C:/Users/s4425412/Documents/Honours/Data'
keyword_dict = {
    'SPEXGO': {'type': '', 'time': '48:00:00', 'vmem': int(10000), 'ncpus': int(16)},  # single-point energy calculation without geometry optimization
    'GOVF': {'type': ' opt int(grid=ultrafine) scf=tight', 'freq': ' freq=noraman', 'mem': int(4000), 'jobfs': int(2400), 'time': '10:00:00'},  # geometry optimization with vibrational frequency calculation
    'TSGOVF': {'type': ' opt=(ts,calcfc,noeigentest,maxcyc=200) int(grid=ultrafine) scf=tight', 'freq': ' freq=noraman', 'mem': int(4000), 'jobfs': int(4000), 'time': '10:00:00'},  # transition state GOVF
    'IRC': {'type': ' irc=(maxpoints=6,stepsize=10,calcfc,maxcyc=200)', 'time': '10:00:00'},  # intrinsic reaction coordinate calculation
    'SPEiS': {'type': ' int(grid=ultrafine) scf=tight', 'freq': '', 'time': '07:00:00', 'mem': int(70000), 'jobfs': int(100000)},  # single-point energy calculation in solvent
    'MO': {'type': ' gfinput iop(6/7=3)', 'time': '10:00:00'},  # print molecular orbitals for Molden display
    'FBRGO': {'type': ' opt=modredundant int(grid=ultrafine) scf=tight', 'freq': ' freq=noraman', 'edit_charge': int(-1), 'time': '18:00:00', 'mem': int(4000), 'jobfs': int(4000)},  # freeze a bond, optimize the rest
    'SBGO': {'type': ' opt=(z-matrix)'},  # stretch a bond, optimize for each point
    }
# theory_lvl_list = ['B2PLYPD/6-311G(2d,p)', 'M062X/6-311G(2d,p)', 'MP2/6-311G(2d,p)',
#                   'wB97XD/6-311G(2d,p)', 'wB97XD/AUG-cc-pVTZ', 'wB97XD/6-311G(d,p)', 'wB97XD/6-311+G(d,p)']
theory_lvl_list = ['B2PLYPD/6-31+G(d)', 'M062X/6-311+G(d,p)', 'M062X/6-311G(d,p)', 'wB97XD/6-311+G(d,p)']
# Correction needed for SCS-MP2 afterwards, function in gaussian.py
