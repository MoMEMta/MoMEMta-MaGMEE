# MoMEMta-MaGMEE
MoMEMta - MadGraph Matrix Element Exporter

**This is the standalone version of the exporter**
Use it to export C++ matrix elements to be used **without** MoMEMta.

## Requirements

- Python >= 2.7
- MG5_aMC@NLO >= 2.5.1 (see [here](https://launchpad.net/~maddevelopers))
- A C++-11 capable compiler

## Install

Once your copy of MG5 is setup, go the the `PLUGIN` subfolder. There, download the MoMEMta-MaGMEE plugin (either by cloning the git repository, or by downloading and extracting the archive). 

You are now ready to export matrix elements in C++!

## Usage

Example:

- Go to the MG5 folder
- Launch MG5: `bin/mg5_aMC`
- Generate your favourite process: `generate p p > t t~, (t > w+ b, w+ > e+ ve), (t~ > w- b~, w- > e- ve~)`
- Output the process: `output MoMEMta_standalone myHappyME`. `myHappyME` will be used to define the namespace in which the matrix element code is enclosed, so it has to be a valid C++ variable name.

You now have a folder called `myHappyME` containing the model's code, a parameter card and the matrix element's code itself. 

To build the code and use it in your project, do:
```
cd myHappyME
mkdir build
cd build
cmake ..
make -j 4
```
This generates a shared library that can be linked with your code.
