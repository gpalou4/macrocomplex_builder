# Macrocomplex builder
*Guillermo Palou Márquez and Javier Sánchez Utgés*

## **TABLE OF CONTENTS**

<!-- TOC depthFrom:1 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->
- [Description](#description)
- [SBI-Python project](#sbi-python-project)
- [Installation](#installation)
  - [via Git](#1-via-git)
  - [via Setup](#2-via-setup)
  - [via Pip](#3-via-pip)
- [Requirements](#Requirements)
- [Information](#information)
  - [Input files](#input-files)
  - [Limitations](#Limitations)
- [Biological background](#background)
- [Algorithm](#algorithm)
- [Tutorial](#tutorial)
  - [Command-line arguments](#command-line-arguments)
  - [Example 1](#example-1-6ezm)
  - [Example 2](#example-2-1g65)
  - [Example 3](#example-3-5vox)
  - [Example 4](#example-4-5oom)
  - [Example 5](#example-5-3kuy)
<!-- /TOC -->

## Description
This program is able to reconstruct biological macrocomplexes of protein-protein interactions as well as protein-DNA/RNA interactions given a set of binary interactions and the desired number of chains of the target complex.

## SBI-Python project

## Installation

### **1) via Git**

The package can be downloaded using Git:

 ```bash
  git clone https://github.com/gpalou4/macrocomplex_builder
  cd macrocomplex_builder
 ```

### **2) via Setup**

### **3) via Pip**

## Tutorial

### **Command-line arguments**

  - -h, --help: this flag will show the usage of the program as well as a description of what it does as well as an explanation of all the parameters it has and can modify or offer some information when executing the program.
  - -i, --input: this argument is **required** can either be an absolute or relative path of the input folder containing all the binary-interaction PDB files that are going to be used to build the complex.
  - -o, --output: this argument is **optional** and if set, all the output files will be saved in this folder. If not set, by default, the output files will be saved in a folder named: _input_foldername_output_.
  - -v, --verbose: this argument is **optional** and will print the progression log in the standard error if set.
  - -pi, --pdb_iterations: this argument is **optional** and will save a new PDB file every time a chain is added to the complex if set.
  - -nc, --number_chains: this argument is **required** and indicates the number of chains that the final complex must have.
  - -rmsd, --rmsd_threshold: this argument is **optional** and if set, the RMSD threshold will take its value. If not, it will take a value of 0.3 by default.
  - -cl, --clashes_theshold: this argument is **optional** and if set, the clashes threshold will take its value. If not, it will take a value of 30 by default.
  - -it, --iterations: this argument is **optional** and if set, the maximum number of iterations will run for, if the number of chains is not reached will take its value. If not set, it will take a value of 100 by default.

### Example 1, 6EZM

The first example is the Imidazoleglycerol-phosphate dehydratase from Saccharomyces cerevisiae. It is a Homo 24-mer (stoichiometry: A24). Based on the input provided files, the following command will recover the complete complex:

```bash
python3 macrocomplex_builder.py -i 6ezm -nc 24
```
 Where:
 
 -6ezm: is the input folder containing all the binary-interaction PDB files
 
 -24: indicates that the final complex must have 24 chains
 
 The computation time is around 2-5 seconds, and the RMDS between the reconstructed complex and the original PDB file is 0.639 Â. The input folder contains 23 files and 24 chains. In this particular example, there is always a common chain (chain A) between any two binary interactions. For that reason, in each iteration occurs a superimposition between both equal chains A with an RMSD of 0, giving a total of 24 iterations and 24 chains.
 
### Example 2, 1G65
  
The second example corresponds to the 20S proteosome from Saccharomyces cerevisiae. It is a Hetero 28-mer (stoichiometry: A2B2C2D2E2F2G2H2I2J2K2L2M2N2). Based on the input provided files, the following command will recover the complete complex:

```bash
python3 macrocomplex_builder.py -i 1g65 -nc 28 -rmsd 0.5 -cl 45 
```
 Where:
 
 -0.5: is the new RMSD threshold (default is 0.3)
 
 -45: is the new clashes threshold (default is 30)
 
The computation time is around 10-15 seconds and the RMDS between the reconstructed complex and the original PDB file is 0.975 Â. In this example, we must change the RMDS and clashes thresholds because it can only recover 27 chains with the default values. The input folder contains 6 files and 8 different chains, and only some of them have a common chain. All superimpositions below the set thresholds are used to add the rotated chains to the complex. Notice that even though there are less number of files/chains than the total number of chains of the original complex, it is able to reconstruct the complete macrocomplex in 27 iterations.

| <img src="images/1G65_ORIGINAL.png" width="200" height="200"> | <img src="images/1G65_BUILT.png" width="200" height="200"> | <img src="images/1G65_SUPERIMPOSED.png" width="200" height="200"> |
| :---: | :---: | :---: |
| *Original* | *Built* | *Superimposition*|

### Example 3, 5VOX

The third example corresponds to the Yeast V-ATPase in complex with Legionella pneumophila effector SidK. It is a Hetero 33-mer  (stoichiometry: A8B3C3D3E3F3GHIJKLMNOP). Based on the input provided files, the following command will recover the complete complex:

```bash
python3 macrocomplex_builder.py -i 5vox -nc 33 
```

The computation time is around 5-10 seconds and the RMDS between the reconstructed complex and the original PDB file is surprisingly 0 Â. The input folder contains 51 files and 33 different chains
 
 
### Example 4, 5OOM
  
```bash
python3 macrocomplex_builder.py -i 5oom -nc 53 
```
 
### Example 5, 3KUY
  
```bash
python3 macrocomplex_builder.py -i 6kuy -nc 10
```
