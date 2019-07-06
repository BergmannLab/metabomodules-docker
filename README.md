# METABOMODULES
This repository holds the source code for **metabomodules**, a Linux/macOS command-line toolbox to for automatic metabolite identification in 1D NMR data using unsupervised clustering.

## SOURCE CODE
To get a copy of the source code run:
```$ git clone --recursive https://github.com/BergmannLab/metabomodules-docker```

## PREREQUISITES
Either ```docker``` or ```singularity``` must be installed. Please visit https://www.docker.com or http://singularity.lbl.gov

The tool was tested on *Ubuntu Linux 18.04*, *CentOS Linux 7.5* and *macOS Sierra* Version 10.12.

## INSTALLATION
To install: ```./install```

To uninstall: ```./uninstall```

## RUNNING
To get help with running the tool, install then invoke, from any location: ```metabomodules --help```
For example: **metabomodules --input=/tmp/input.csv --container=docker**

## INPUT
In input, provide pre-processed (aligned) 1D NMR data in a tabular form (peak list)
* firt row: PPM axis
* first column: sample id
* each cell i,j, contains the area under the peak found at PPM i for sample j

## Methods
The following methods will be applied to the input data in order to automatically identify metabolites
* **ISA**: Iterative Signature Algorithm
* **ACP**: Average Correlation Profiles Method
* **PCA**: Principal Component Analysis

## PUBLICATION
Khalili, Bita & Tomasoni, Mattia & Mattei, Mirjam & Mallol Parera, Roger & Sonmez, Reyhan & Krefl, Daniel & Rueedi, Rico & Bergmann, Sven. (2019). Automated analysis of large-scale NMR data generates metabolomic signatures and links them to candidate metabolites. 10.1101/613935.
