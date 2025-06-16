#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate spec_wrapper

unset HDF5 HDF5_ROOT HDF5_HOME FFTW FFTW_DIR

#pip install  -v -e . 2>&1 | tee compile.log

rm -rf build

cmake -S . -B build

cd build

make -j4
