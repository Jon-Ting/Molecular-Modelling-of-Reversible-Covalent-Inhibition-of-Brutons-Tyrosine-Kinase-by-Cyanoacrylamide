"""
Python script to store the setting variables for other scripts
Written by Jonathan Yik Chang Ting (Student ID: 44254124) for UQ BAdSc(Hons) Honours Project 2019.
Project Title: Molecular Modelling of Reversible Covalent Inhibition of Bruton’s Tyrosine Kinase by Cyanoacrylamides
Free to be distributed, edited and utilized for purposes beneficial to the world. :)
Though recognition of my effort is appreciated! XD
"""

keyword_dict = {
    'GOVF': {'type': ' opt int(grid=ultrafine) scf=tight', 'freq': ' freq=noraman', 'mem': int(8000), 'jobfs': int(2400), 'time': '10:00:00'},  # geometry optimization with vibrational frequency calculation
    'TSGOVF': {'type': ' opt=(ts,calcfc,noeigentest,maxcyc=200) int(grid=ultrafine) scf=tight', 'freq': ' freq=noraman', 'mem': int(4000), 'jobfs': int(4000), 'time': '10:00:00'},  # transition state GOVF

    'SPEXGO': {'type': '', 'time': '48:00:00', 'vmem': int(10000), 'ncpus': int(16)},  # single-point energy calculation without geometry optimization
    'IRC': {'type': ' irc=(maxpoints=6,stepsize=10,calcfc,maxcyc=200)', 'time': '10:00:00'},  # intrinsic reaction coordinate calculation
    'SPEiS': {'type': ' int(grid=ultrafine) scf=tight', 'freq': '', 'time': '07:00:00', 'mem': int(70000), 'jobfs': int(100000)},  # single-point energy calculation in solvent
    'MO': {'type': ' gfinput iop(6/7=3)', 'time': '10:00:00'},  # print molecular orbitals for Molden display
    'FBRGO': {'type': ' opt=modredundant int(grid=ultrafine) scf=tight', 'freq': ' freq=noraman', 'edit_charge': int(-1), 'time': '18:00:00', 'mem': int(4000), 'jobfs': int(4000)},  # freeze a bond, optimize the rest
    'SBGO': {'type': ' opt=(z-matrix)'},  # stretch a bond, optimize for each point
    }

