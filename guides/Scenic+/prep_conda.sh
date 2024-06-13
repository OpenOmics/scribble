#! /bin/bash

# This follows the instructions for setting up a conda environment on Biowulf as of 2/16/24.
# Make sure to uses the package from the github site and only that as it ensures you have the correct versions of all required packages.

# You don't need this if using the singularity object.

python3 -m venv .scenicPlus_env
source .scenicPlus_env/bin/activate
pip install --upgrade pip

git clone https://github.com/aertslab/scenicplus
pip install -e scenicplus

deactivate
