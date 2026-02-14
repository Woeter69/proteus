#!/bin/bash
source /opt/anaconda/etc/profile.d/conda.sh
conda env remove -n proteus_env -y
conda env create -f environment.yml
