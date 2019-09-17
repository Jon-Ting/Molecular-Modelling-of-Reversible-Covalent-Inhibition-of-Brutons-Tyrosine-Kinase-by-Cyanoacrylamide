"""
Python script to store the configuration for plotting functions in plot_fig.py.
Written by Jonathan Yik Chang Ting (Student ID: 44254124) for UQ BAdSc(Hons) Honours Project 2019.
Project Title: Molecular Modelling of Reversible Covalent Inhibition of Bruton’s Tyrosine Kinase by Cyanoacrylamides
Free to be distributed, edited and utilized for purposes beneficial to the world. :)
Though recognition of my effort is appreciated! XD
"""


charge_list = ["Mulliken", "NBO", "MK", "Hirshfeld", "CM5", "AIM", "ChelpG", "Omega"]

combination_dict = {
    "I": {
        "LUMO": {
            "M": 0.6121,
            "C": 22.143,
            "R2": 0.8354,
            "X-AXIS": 0,
            "Y-AXIS": 0.4,
            "LEG": "upper left",
            "ALIGN": "center",
            "FONTSIZE": "medium",
            "X-NAME": "Ligand LUMO Energy (kcal/mol)"
        },

        "Charge": {
            "Mulliken": {
                "M": -3.5523,
                "C": 10.639,
                "R2": 0.2152,
                "X-AXIS": -0.03,
                "Y-AXIS": 0,
                "LEG": "upper right",
                "ALIGN": "right",
                "FONTSIZE": "x-small",
                "X-NAME": "Mulliken Charge (e)"
            },
            "NBO": {
                "M": -24.802,
                "C": 9.7637,
                "R2": 0.7967,
                "X-AXIS": -0.01,
                "Y-AXIS": 0,
                "LEG": "upper right",
                "ALIGN": "right",
                "FONTSIZE": "x-small",
                "X-NAME": "NBO Charge (e)"
            },
            "MK": {
                "M": -7.1748,
                "C": 11.747,
                "R2": 0.3866,
                "X-AXIS": 0.02,
                "Y-AXIS": 0,
                "LEG": "upper right",
                "ALIGN": "left",
                "FONTSIZE": "x-small",
                "X-NAME": "Merz-Kollman Charge (e)"
            },
            "Hirshfeld": {
                "M": -89.348,
                "C": 12.659,
                "R2": 0.8491,
                "X-AXIS": 0.005,
                "Y-AXIS": 0,
                "LEG": "upper right",
                "ALIGN": "left",
                "FONTSIZE": "x-small",
                "X-NAME": "Hirshfeld Charge (e)"
            },
            "CM5": {
                "M": -54.454,
                "C": 8.8508,
                "R2": 0.7741,
                "X-AXIS": -0.005,
                "Y-AXIS": 0,
                "LEG": "lower left",
                "ALIGN": "right",
                "FONTSIZE": "x-small",
                "X-NAME": "CM5 Charge (e)"
            },
            "AIM": {
                "M": -92.849,
                "C": 14.411,
                "R2": 0.8573,
                "X-AXIS": -0.003,
                "Y-AXIS": 0,
                "LEG": "upper right",
                "ALIGN": "right",
                "FONTSIZE": "x-small",
                "X-NAME": "AIM Charge (e)"
            },
            "ChelpG": {
                "M": -19.394,
                "C": 9.9205,
                "R2": 0.6864,
                "X-AXIS": -0.004,
                "Y-AXIS": 0.2,
                "LEG": "upper right",
                "ALIGN": "right",
                "FONTSIZE": "x-small",
                "X-NAME": "ChelpG Charge (e)"
            },
            "Omega": {
                "M": -0.7225,
                "C": 56.999,
                "R2": 0.9376,
                "X-AXIS": -0.4,
                "Y-AXIS": 0,
                "LEG": "upper right",
                "ALIGN": "right",
                "FONTSIZE": "x-small",
                "X-NAME": "Electrophilicity Index"
            },
        },
        "DI": {
            "Thiolate": {
                "M": 17.241,
                "C": 7.6018,
                "R2": 0.798,
                "X-AXIS": -0.02,
                "Y-AXIS": 0,
                "LEG": "upper left",
                "ALIGN": "right",
                "FONTSIZE": "x-small",
                "X-NAME": "Thiolate Distortion Energy (kcal/mol)"
            },
            "Ligand": {
                "M": 0.8082,
                "C": 5.3355,
                "R2": 0.9760,
                "X-AXIS": -1.0,
                "Y-AXIS": 0.7,
                "LEG": "best",
                "ALIGN": "left",
                "FONTSIZE": "x-small",
                "X-NAME": "Ligand Distortion Energy (kcal/mol)"
            },
            "Activation": {
                "M": 0.9013,
                "C": 10.439,
                "R2": 0.8943,
                "X-AXIS": -0.3,
                "Y-AXIS": 0,
                "LEG": "upper left",
                "ALIGN": "right",
                "FONTSIZE": "x-small",
                "X-NAME": "Activation Energy (kcal/mol)"
            },
            "Interaction": {
                "M": -1.6397,
                "C": 0.1818,
                "R2": 0.428,
                "X-AXIS": 0.18,
                "Y-AXIS": 0,
                "LEG": "upper right",
                "ALIGN": "left",
                "FONTSIZE": "x-small",
                "X-NAME": "Interaction Energy (kcal/mol)"
            },
        }
    },
    "G": {
        "LUMO": {
                "M": -3.5523,
                "C": 10.639,
                "R2": 0.2152,
                "X-AXIS": -0.03,
                "Y-AXIS": 0,
                "LEG": "upper right",
                "FONTSIZE": "x-small"
            },
        "Charge": {
            "Mulliken": {
                "M": -3.5523,
                "C": 10.639,
                "R2": 0.2152,
                "X-AXIS": -0.03,
                "Y-AXIS": 0,
                "LEG": "upper right",
                "FONTSIZE": "x-small"
            },
            "NBO": {
                "M": -3.5523,
                "C": 10.639,
                "R2": 0.2152,
                "X-AXIS": -0.03,
                "Y-AXIS": 0,
                "LEG": "upper right",
                "FONTSIZE": "x-small"
            },
            "MK": {
                "M": -3.5523,
                "C": 10.639,
                "R2": 0.2152,
                "X-AXIS": -0.03,
                "Y-AXIS": 0,
                "LEG": "upper right",
                "FONTSIZE": "x-small"
            },
            "Hirshfeld": {
                "M": -3.5523,
                "C": 10.639,
                "R2": 0.2152,
                "X-AXIS": -0.03,
                "Y-AXIS": 0,
                "LEG": "upper right",
                "FONTSIZE": "x-small"
            },
            "CM5": {
                "M": -3.5523,
                "C": 10.639,
                "R2": 0.2152,
                "X-AXIS": -0.03,
                "Y-AXIS": 0,
                "LEG": "upper right",
                "FONTSIZE": "x-small"
            },
            "ChelpG": {
                "M": -3.5523,
                "C": 10.639,
                "R2": 0.2152,
                "X-AXIS": -0.03,
                "Y-AXIS": 0,
                "LEG": "upper right",
                "FONTSIZE": "x-small"
            },
            "AIM": {
                "M": -3.5523,
                "C": 10.639,
                "R2": 0.2152,
                "X-AXIS": -0.03,
                "Y-AXIS": 0,
                "LEG": "upper right",
                "FONTSIZE": "x-small"
            },
            "Omega": {
                "M": -3.5523,
                "C": 10.639,
                "R2": 0.2152,
                "X-AXIS": -0.03,
                "Y-AXIS": 0,
                "LEG": "upper right",
                "FONTSIZE": "x-small"
            },
        },
        "DI": {
            "Thiolate": {
                "M": -3.5523,
                "C": 10.639,
                "R2": 0.2152,
                "X-AXIS": -0.03,
                "Y-AXIS": 0,
                "LEG": "upper right",
                "FONTSIZE": "x-small"
            },
            "Ligand": {
                "M": -3.5523,
                "C": 10.639,
                "R2": 0.2152,
                "X-AXIS": -0.03,
                "Y-AXIS": 0,
                "LEG": "upper right",
                "FONTSIZE": "x-small"
            },
            "Activation": {
                "M": -3.5523,
                "C": 10.639,
                "R2": 0.2152,
                "X-AXIS": -0.03,
                "Y-AXIS": 0,
                "LEG": "upper right",
                "FONTSIZE": "x-small"
            },
            "Interaction": {
                "M": -3.5523,
                "C": 10.639,
                "R2": 0.2152,
                "X-AXIS": -0.03,
                "Y-AXIS": 0,
                "LEG": "upper right",
                "FONTSIZE": "x-small"
            },
        }
    }
}