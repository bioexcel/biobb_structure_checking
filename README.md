# biobb_structure_checking
### Introduction
Biobb_structure_checking performs a checking of the quality of a 3D structure intended to facilitate the setup of molecular dynamics simulation of protein or nucleic acids systems.

Biobb_structure_checking package allows to configure the system (selection of model/chains,alternative location, addition of disulfide bonds and hydrogen atoms, side chain mutations), to detect and fix structure errors (missing side chain atoms, backbone breaks, amide assignments, incorrect chirality). It works with structures obtained from the Protein Data Bank or user provided.

The Biobb_structure_checking package provides a command line utility ([check_structure](https://biobb-structure-checking.readthedocs.io/en/latest/command_line_usage.html)) and a python [API](https://biobb-structure-checking.readthedocs.io/en/latest/biobb_structure_checking.html).

The latest documentation of this package can be found in our readthedocs site:
[latest package documentation](http://biobb_structure_checking.readthedocs.io/en/latest/).

### Version
v3.9.7 December 2021

### Requirements

* Biopython
### Optional requirements
* psutil (required for --debug, included in conda pkg.)
* Modeller (required for some functionalities, not included in conda pkg.)
* jupyter & nglview (required for demonstration notebooks, not included in conda pkg.)

### Installation
Using PIP:

> **Important:** PIP only installs the package. All the dependencies must be installed separately. To perform a complete installation, please use ANACONDA.

* Installation:

        pip install "biobb_structure_checking>=3.9.7"

* Usage: [Python API documentation](https://biobb_structure_checking.readthedocs.io/en/latest/modules.html).

Using ANACONDA:

* Installation:

        conda install -c bioconda "biobb_structure_checking>=3.9.7"

* Usage: With conda installation BioBBs can be used with the [Python API documentation](https://biobb_structure_checking.readthedocs.io/en/latest/modules.html) and the  [Command Line documentation](https://biobb_structure_checking.readthedocs.io/en/latest/command_line.html)

### Copyright & Licensing
This software has been developed in the MMB group (http://mmb.irbbarcelona.org) at the
BSC (https://www.bsc.es/) & IRB (https://www.irbbarcelona.org/) for the European BioExcel (https://bioexcel.eu/), funded by the European Commission
(EU H2020 [675728](https://cordis.europa.eu/projects/675728)).

* (c) 2015-2021 [Barcelona Supercomputing Center](https://www.bsc.es/)
* (c) 2015-2021 [Institute for Research in Biomedicine](https://www.irbbarcelona.org/)

Licensed under the
[GNU Lesser General Public License v2.1](https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html).

![](https://bioexcel.eu/wp-content/uploads/2015/12/Bioexcell_logo_1080px_transp.png "Bioexcel")
