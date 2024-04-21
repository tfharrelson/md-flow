# md_flow
-----------
## Summary
`md_flow` is a python package that utilizes dask as a distributed runtime environment to run convenient molecular dynamics workflows. It can be used as a 
standard program (e.g. via `python -m ...`), or it may be used as an SDK. Standard end-to-end flows are provided, and 
each of the step functions are modularly defined. Thus, they may be used as part of user-defined custom flows. Each provided step function is under a
dask `delayed` decorator, meaning that each step will return extremely quickly, and provide a graph of the workflow without evaluating it directly.
These graphs can be composed together or extended to create more complex delayed workflows. Since the flows haven't been evaluated yet, we can
run the computations on distributed resources, when advantageous.

The main workflow developed is an end-to-end molecular dynamics flow that takes a UniProt accession ID, gets a structure using the AlphaFold endpoint, then it creates a solvated 
molecular dynamics box, runs an MD simulation (after doing proper equilibration and setup).

## Installation Instructions
I'm currently working on a docker image and build script to expedite the installation b/c right now you have to install gromacs by hand, which is tricky to get right. 
However, currently the instructions for installing from sources are as follows:

- Install gromacs from source
    - Make sure you have mpich installed (e.g. using `sudo apt install mpich` in ubuntu); I could not get the python bindings to work without mpi4py, which requires mpich.
    - Download the tar ball, and untar
    - Follow the "dirty" install instructions on the [gromacs website](https://manual.gromacs.org/documentation/current/install-guide/index.html).
    - NOTE: you may add extra flags to the cmake that are specific to your system/preferences (e.g. I enabled cuda).
- Create a new python virtual environment (using conda, venv, etc.).
- Activate the environment (e.g. `conda activate my_env`).
- Run `pip install -e .` from the top level directory.
    - This will install `md_flow`, `gmxapi`, and their dependencies.
 
## Getting started
In a jupyter notebook, you can import the control structure that abstracts away the cluster information, and run a flow. See the following:

```python
from md_flow import MDCluster

cluster = MDCluster()
result = cluster.npt_md("P00250")   # a Uniprot ID for a protein

result
```
This should run a bunch of molecular dynamics simulations culminating in a "production" run using the NPT ensemble (particle <u>N</u>umber, <u>P</u>ressure, and <u>T</u>emperature held constant).
The result contains all the relevant metadata contained in the molecular dynamics run along with the results, which include the trajectory file, and the energy file. Each of these can be downstream
processed by additional analysis scripts (I recommend using MDAnalysis).

