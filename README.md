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
  - -rm, --rmsd: this argument is **optional** and if set, the RMSD threshold will take its value. If not, it will take a value of 0.3 by default.
  - -cl, --clashes: this argument is **optional** and if set, the clashes threshold will take its value. If not, it will take a value of 30 by default.
  - -it, --iterations: this argument is **optional** and if set, the maximum number of iterations will run for, if the number of chains is not reached will take its value. If not set, it will take a value of 100 by default.

### Example 1, 6EZM

The first example is the Imidazoleglycerol-phosphate dehydratase from Saccharomyces cerevisiae. It is a Homo 24-mer (stoichiometry: A24). Based on the input provided files, the following command will recover the complete complex:
```bash
python3 macrocomplex_builder -i 6ezm -nc 24
```
 Where:
 -6ezm: is the input folder containing all the binary-interaction PDB files
 -24: indicates that the final complex must have 24 chains
 The computation time is around 3 seconds, and the RMDS between the reconstructed complex and the original PDB file is 0.639 Â.
 
### Example 2, 1G65
  
```bash
python3 macrocomplex_builder -i 1g65 -nc 28 -rmsd 0.5 -cl 45 
```
 
### Example 3, 5VOX
  
```bash
python3 macrocomplex_builder -i 5vox -nc 33 
```
 
### Example 4, 5OOM
  
```bash
python3 macrocomplex_builder -i 5oom -nc 53 
```
 
### Example 5, 3KUY
  
```bash
python3 macrocomplex_builder -i 6kuy -nc 10
```
