# MoMEMta-MaGMEE
MoMEMta - MadGraph Matrix Element Exporter

**Note**: this version of the exporter is only compatible with MoMEMta v1.0.X. For newer versions of MoMEMta, select the appropriate branch of this repository.

## Requirements

- Python >= 2.7
- MG5_aMC@NLO >= 2.5.0 (see [here](https://launchpad.net/~maddevelopers))
- MoMEMta v1.0.X
- A C++-11 capable compiler

**Note**: MoMEMta needs to be installed on the system (locally or globally), cf. MoMEMta documentation.

## Install

Once your copy of MG5 is setup, go the the `PLUGIN` subfolder. There, download the MoMEMta-MaGMEE plugin (either by cloning the git repository, or by downloading and extracting the archive). 

You are now ready to export matrix elements to use in MoMEMta!

## Usage

Example:

- Go to the MG5 folder
- Launch MG5: `bin/mg5_aMC`
- Generate your favourite process: `generate p p > t t~, (t > w+ b, w+ > e+ ve), (t~ > w- b~, w- > e- ve~)`
- Output the process: `output MoMEMta myHappyME`. `myHappyME` will be used to define the namespace in which the matrix element code is enclosed, so it has to be a valid C++ variable name.

You now have a folder called `myHappyME` containing the model's code, a parameter card and the matrix element's code itself. 

To build the code and use it in your project, you need to have MoMEMta installed on your system (see the MoMEMta documentation for more details). Then, do:
```
cd myHappyME
mkdir build
cd build
cmake ..
make -j 4
```
This generates a shared library that can be dynamically loaded by MoMEMta (using the `load_modules()` function in the Lua config file).

The following options are available when configuring the the build (when running `cmake ..`):
- `-DCMAKE_PREFIX_PATH=(path)`: Path to the installation of MoMEMta. **Must be specified** if your version of MoMEMta is installed locally and not in your system directories

The matrix element has a name assigned to it, of type `myHappyME_P1_Sigma_pp_ttx_...`: 
this is the name to use when defining the matrix element module in your Lua config file. 
If you wish to modify it, simply change the name string at the end of
the .cc file in the corresponding subfolder of the `SubProcesses` folder before building the code.
