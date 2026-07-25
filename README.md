# Enhancing Presolve in Mixed Integer Programming by Combining Probing and Dual Fixing
## Authors: Zhao-Wei Wang, Wei-Kun Chen, Yu-Hong Dai


## Project Overview
This project provides a version of the HiGHS solver (https://github.com/ERGO-Code/HiGHS) containing an implementation of the two presolve techniques proposed in the paper [Enhancing Presolve in Mixed Integer Programming by Combining Probing and Dual Fixing](https://arxiv.org/abs/2607.10767) by Zhao-Wei Wang, Wei-Kun Chen, and Yu-Hong Dai.


## Repository Organization
* ```highs-1.14.0/```: the HiGHS solver with the two proposed presolve techniques. 
    We modified the code based on the release [HiGHS 1.14.0](https://github.com/ERGO-Code/HiGHS/releases/tag/v1.14.0).
    The new presolve techniques are implemented in ```highs-1.14.0/highs/presolve/HPresolve.cpp```, ```highs-1.14.0/highs/mip/HighsDomain.cpp```, and ```highs-1.14.0/highs/mip/HighsImplication.cpp```
* ```paper-settings/:```: settings files containing parameters to enable the two proposed presolve techniques in the paper.
* ```results.csv```: detailed instance-wise results of the experiments in the paper.

## Installation
* Simply install HiGHS as usual:
```shell
cd highs-1.14.0/
cmake -S . -B build
cmake --build build
```
* More details about the installation of HiGHS can be found at the [official installation page](https://github.com/ERGO-Code/HiGHS#installation).


## Testset
 * The MIPLIB 2017 benchmark is used as the testset.
 * We use 5 random seeds for each of the 240 problems, and each problem and seed combination is treated as an individual observation, referred to as an "instance". 

## Running experiments
* To solve an instance with different settings, use
```shell
# run Default
./build/bin/highs /path/to/your/instance --options_file ../paper-settings/Default.setting
# run DF-Probing
./build/bin/highs /path/to/your/instance --options_file ../paper-settings/DF-Probing.setting
# run GDF
./build/bin/highs /path/to/your/instance --options_file ../paper-settings/GDF.setting
# run All
./build/bin/highs /path/to/your/instance --options_file ../paper-settings/All.setting
```
* Detailed instance-wise results are given in the csv file ```results.csv```
