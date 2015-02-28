# Model of germinal center B-cell affinity maturation

## System Requirements

Requires Python 2.7.x, SQLite, and a C++11-compliant compiler.

The build script (`make.py) assumes a Unix system and have been tested on Mac OS X and Linux.
It should be possible to compile the code on Windows, but not with the supplied build script.

The build script will check for the Intel C++ compiler, LLVM/Clang, and GCC, preferentially in that order.

## Getting and building the model

Run these Git commands to get the repository:

```sh
git clone https://github.com/cobeylab/gcdynamics.git
cd gcdynamics
```

To download all necessary libraries and build the 

```sh
./make.py
```

and an executable named `gcdynamics` will appear in the `bin/` directory.

## Performing a single model run

The model requires a JSON-formatted parameters file containing values for all parameters.
To prevent mistakes, no defaults are provided in the executable, and the code will fail if any parameters are missing from the parameters file.
An example parameters file, with descriptions of all parameters, is provided in `example/constant_parameters.json`.

To perform a run:
1. Create a working directory for the run
2. Copy the `example/constant_parameters.json` file to the working directory and modify the parameter values.
3. Run the code inside the working directory.

The whole thing would look something like this at the command line, along with output redirection for stdout and stderr:

```
#!sh
mkdir <run-dir>
cp <gcdynamics-dir>/example/parameters.json <run-dir>/parameters.json
cd <run-dir>
# [edit <run-dir>/parameters.json]
[path-to-gcdynamics]/gcdynamics parameters.json 2> stderr.txt 1> stdout.txt
```
