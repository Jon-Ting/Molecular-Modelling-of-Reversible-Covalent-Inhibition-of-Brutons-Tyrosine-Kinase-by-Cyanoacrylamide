"""
Python script to store the configuration for plotting functions in plot_fig.py.
Written by Jonathan Yik Chang Ting (Student ID: 44254124) for UQ BAdSc(Hons) Honours Project 2019.
Project Title: Molecular Modelling of Reversible Covalent Inhibition of Bruton’s Tyrosine Kinase by Cyanoacrylamides
Free to be distributed, edited and utilized for purposes beneficial to the world. :)
Though recognition of my effort is appreciated! XD
"""


charge_list = ["Mulliken", "NBO", "MK", "Hirshfeld", "CM5", "AIM", "ChelpG", "Omega"]
DI_list = ["Thiolate", "Ligand", "Activation", "Interaction"]

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
            "M": 0.6341,
            "C": 11.665,
            "R2": 0.8038,
            "X-AXIS": 0,
            "Y-AXIS": 0.4,
            "LEG": "upper left",
            "ALIGN": "center",
            "FONTSIZE": "medium",
            "X-NAME": "Ligand LUMO Energy (kcal/mol)"
        },

        "Charge": {
            "Mulliken": {
                "M": -3.8031,
                "C": 13.26,
                "R2": 0.2089,
                "X-AXIS": -0.04,
                "Y-AXIS": 0,
                "LEG": "upper right",
                "ALIGN": "right",
                "FONTSIZE": "x-small",
                "X-NAME": "Mulliken Charge (e)"
            },
            "NBO": {
                "M": -24.677,
                "C": 12.381,
                "R2": 0.678,
                "X-AXIS": -0.01,
                "Y-AXIS": 0,
                "LEG": "upper right",
                "ALIGN": "right",
                "FONTSIZE": "x-small",
                "X-NAME": "NBO Charge (e)"
            },
            "MK": {
                "M": -15.228,
                "C": 12.617,
                "R2": 0.7188,
                "X-AXIS": 0.02,
                "Y-AXIS": 0,
                "LEG": "upper right",
                "ALIGN": "left",
                "FONTSIZE": "x-small",
                "X-NAME": "Merz-Kollman Charge (e)"
            },
            "Hirshfeld": {
                "M": -85.733,
                "C": 15.179,
                "R2": 0.7221,
                "X-AXIS": 0.004,
                "Y-AXIS": 0,
                "LEG": "lower left",
                "ALIGN": "left",
                "FONTSIZE": "x-small",
                "X-NAME": "Hirshfeld Charge (e)"
            },
            "CM5": {
                "M": -52.002,
                "C": 11.53,
                "R2": 0.6279,
                "X-AXIS": -0.005,
                "Y-AXIS": 0,
                "LEG": "lower left",
                "ALIGN": "right",
                "FONTSIZE": "x-small",
                "X-NAME": "CM5 Charge (e)"
            },
            "AIM": {
                "M": -95.9,
                "C": 15.637,
                "R2": 0.7617,
                "X-AXIS": -0.004,
                "Y-AXIS": 0,
                "LEG": "upper right",
                "ALIGN": "right",
                "FONTSIZE": "x-small",
                "X-NAME": "AIM Charge (e)"
            },
            "ChelpG": {
                "M": -23.509,
                "C": 12.432,
                "R2": 0.8589,
                "X-AXIS": -0.007,
                "Y-AXIS": 0.18,
                "LEG": "upper right",
                "ALIGN": "right",
                "FONTSIZE": "x-small",
                "X-NAME": "ChelpG Charge (e)"
            },
            "Omega": {
                "M": -0.9655,
                "C": 61.066,
                "R2": 0.9274,
                "X-AXIS": -0.3,
                "Y-AXIS": 0,
                "LEG": "upper right",
                "ALIGN": "right",
                "FONTSIZE": "x-small",
                "X-NAME": "Electrophilicity Index"
            },
        },
        "DI": {
            "Thiolate": {
                "M": 16.642,
                "C": 10.236,
                "R2": 0.6431,
                "X-AXIS": -0.02,
                "Y-AXIS": 0,
                "LEG": "upper left",
                "ALIGN": "right",
                "FONTSIZE": "x-small",
                "X-NAME": "Thiolate Distortion Energy (kcal/mol)"
            },
            "Ligand": {
                "M": 0.9119,
                "C": 6.9095,
                "R2": 0.9469,
                "X-AXIS": -1.4,
                "Y-AXIS": 0.55,
                "LEG": "best",
                "ALIGN": "left",
                "FONTSIZE": "x-small",
                "X-NAME": "Ligand Distortion Energy (kcal/mol)"
            },
            "Activation": {
                "M": 0.9652,
                "C": 10.666,
                "R2": 0.9044,
                "X-AXIS": -0.35,
                "Y-AXIS": 0.24,
                "LEG": "upper left",
                "ALIGN": "right",
                "FONTSIZE": "x-small",
                "X-NAME": "Activation Energy (kcal/mol)"
            },
            "Interaction": {
                "M": -2.1245,
                "C": 4.387,
                "R2": 0.2975,
                "X-AXIS": 0.1,
                "Y-AXIS": 0,
                "LEG": "lower left",
                "ALIGN": "left",
                "FONTSIZE": "x-small",
                "X-NAME": "Interaction Energy (kcal/mol)"
            },
        }
    },
}
