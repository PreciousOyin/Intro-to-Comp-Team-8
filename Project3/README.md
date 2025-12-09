# README.md for all four tasks
This directory contains four MATLAB tasks that reproduce and analyze results from the Stefanini paper. Each task runs independently and generates the associated figures or population simulations.

## Requirements
- For task 1-4, preinstall parallel computing toolbox
- To run the code each task, just click on green triangle run button in matlab, no further input is needed

## Stefanini_driver_task1and2.m
- Computes and provides visual for Figs 1A, S1A, S2A, and Figs 1B, S1B, S2B from the Stefanini paper
- Computes and provides visual for Figs 1C, S1C, S2C and Figs 1D,S1D, S2D from the Stefanini paper
- Figures are generated automatically and displayed in MATLAB.

## Stefanini_driver_task3.m
- PopPK
- simulates the impact of variability in kcl,A (clearance rate constant of the antibody) and in Vr (‘rest of body’ volume) with 100 virtual patients.
- This assummed a normal distribution for both parameters, coefficients of variation of 50% for kcl,A and 25% for Vr
- This visualized the concentrations of VEGF, VEGF-VEGFR2 Complex, VEGFR2, VEGF-Antibody Complex for both tumor and tissue compartments
- Figures are generated automatically and displayed in MATLAB.
- Caluclates the Average AUC for 100 patients


## Stefanini_driver_task1and4.m
- PopPD
- simulated the impact of VEGF, VEGFR2 expression variability for 100 specific patients in VirtualPatientsTask4.csv
- This visualized the concentrations of VEGF, VEGF-VEGFR2 Complex, VEGFR2, VEGF-Antibody Complex for both tumor and tissue compartments
- Figures are generated automatically and displayed in MATLAB.
- Caluclates the Average AUC for 100 specific patients
