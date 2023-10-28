#!/bin/bash

source /opt/conda/bashrc && micromamba activate \
    && micromamba activate pygplates-tutorials\
    && jupyter notebook --allow-root --ip=0.0.0.0 --no-browser