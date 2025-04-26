# GIMIC Automation

A suite of Bash scripts designed to automate the initialization of Python virtual environments and streamline the execution of Gauge-Including Magnetically Induced Current (GIMIC) calculations. These scripts are specifically tailored to workflows utilizing Gaussian-generated input files.

> ⚠️ **Note:** Currently, only Gaussian output files (`.fchk` or `.log`) are supported.

## Overview

This repository facilitates the setup and deployment of [GIMIC](https://github.com/qmcurrents/gimic.git), a computational tool for analyzing magnetically induced current densities in molecular systems. The provided scripts handle virtual environment creation, dependency installation, and the orchestration of GIMIC execution steps.

## Contents

- `GIMIC_virtual-env_automation.sh`  (use this if you have installed GIMIC in virtual environment)
  
  Automates the creation of a Python virtual environment and installs all required dependencies for GIMIC.  
  > 🔧 *Make sure to modify the `ACTIVATE_SCRIPT` variable in the script to reflect the correct path to your virtual environment, e.g.:*  
  `ACTIVATE_SCRIPT="/home/path/to/your/gimic_env/bin/activate"`

- `GIMIC_automation.sh`  (use this if you have installed GIMIC system-wide)
  
  Manages the execution of GIMIC computations using Gaussian output files and a predefined `gimic.inp` configuration file. This script supports automated and batch-mode workflows.

## Usage


1. Ensure either desired scripts are executable: 

       
        chmod +x GIMIC_virtual-env_automation.sh
        chmod +x GIMIC_automation.sh

3. Place the following files in the working directory:

   > A Gaussian output file (*.fchk or *.log)

   > A valid gimic.inp configuration file

   > One of the automation scripts

3.Run the desired script:
 
  # To set up the GIMIC cakculation in virtual environment
    bash GIMIC_virtual-env_automation.sh

  # To execute the GIMIC calculation
    bash GIMIC_automation.sh


4.Requirements

    Bash shell (Linux or Unix-based system)

    GIMIC source code and compatible toolchain (e.g., CMake, compilers)

    Gaussian output files (.fchk or .log)

License

This project is distributed under the MIT License.
