# This repository contains all files (structures, code, simulation input details etc.), required for re-producing water vapor adsorption studies #

![TOC_new](https://github.com/user-attachments/assets/eab96bcd-b949-4d1d-9e3d-c6c0f62b124d)


## All data is provided along with visualization codes in the folders ##
Example --- For electrostatic potential analysis, i) Calculate the potential using ewald_potential.py, and ii) plot the potential_data using the visualization code provided there.
Also, please look at repo_structure.txt for folder-file structure tree (example below for adsorption_data/10.78A/desorption_q_0.7830/):
.
├── README.md \
├── adsorption_data \
│   ├── 10.78A              \
│   │   ├── desorption_q_0.7830
│   │   │   ├── AA.csv
│   │   │   ├── AAR.csv
│   │   │   ├── AR.csv
│   │   │   ├── Inner.csv
│   │   │   ├── Outer.csv
│   │   │   ├── Visualize_Isotherm.ipynb
│   │   │   └── files.zip
│   │   ├── q_0.3915
│   │   │   ├── AA.csv
│   │   │   ├── AAR.csv
│   │   │   ├── AR.csv
│   │   │   ├── Inner.csv
│   │   │   ├── Outer.csv
│   │   │   └── Visualize_Isotherm.ipynb
│   │   ├── q_0.7830
│   │   │   ├── AA.csv
│   │   │   ├── AAR.csv
│   │   │   ├── AR.csv
│   │   │   ├── Inner.csv
│   │   │   ├── Outer.csv
│   │   │   └── Visualize_Isotherm.ipynb

NOTE: In energy analysis, files with extension _Ads.csv are for adsorbate-adsorbate potential energy while files with _E.csv extension are the host-adsorbate energies
