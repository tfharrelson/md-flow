# MD_FLOW
-----------
## Summary
md_flow is a python package that utilizes metaflow to run convenient molecular dynamics workflows. It can be used as a 
standard program (e.g. via `python -m ...`), or it may be used as an SDK. Standard end-to-end flows are provided, and 
each of the step functions are modularly defined. Thus, they may be used as part of user-defined custom flows.

The main workflow developed is an end-to-end molecular dynamics flow that takes a UniProt accession ID, creates a solvated 
molecular dynamics box, runs an MD simulation

## Installation Instructions
I don't know exactly how to write down the install instructions yet b/c they're pretty complicated.
I plan to make a docker file to automate the construction of an image that will always work. So, the 
instructions for installing from sources are as follows:

- Install gromacs from source
    - Make sure you have mpich installed (e.g. using `sudo apt install mpich` in ubuntu); I could not get the python bindings to work without mpi4py, which requires mpich.
    - Download the tar ball, and untar
    - Follow the "dirty" install instructions on the [gromacs website](https://manual.gromacs.org/documentation/current/install-guide/index.html).
    - NOTE: you may add extra flags to the cmake that are specific to your system/preferences (e.g. I enabled cuda).
- Create a new python virtual environment (using conda, venv, etc.).
- Activate the environment (e.g. `conda activate my_env`).
- Run `pip install -r requirements.txt` from the top level directory.
    - This will install md_flow, gmxapi, and their dependencies.

