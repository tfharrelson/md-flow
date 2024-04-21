import requests
import logging
import os
import gmxapi
from md_flow import md_inputs
from typing import Any
from dask import delayed
from .models import MDRunInput, MDRun, ProteinInput


logger = logging.getLogger(__file__)


@delayed
def get_alphafold_pdb(uniprot_id: str) -> str:
    """
    Helper function that takes a uniprot ID and returns a path to a PDB file.

    Parameters:
        uniprot_id (str): A string representing a unique ID of a protein
            existingin the uniprot database

    Returns:
        pdb_string (str): the string (stored in RAM) to the PDB file
            representing a structural instance of the protein encoded by the
            uniprot ID
    """
    url = "https://alphafold.ebi.ac.uk/api"
    route = f"/prediction/{uniprot_id}"

    response = requests.get(url + route)
    response_info = response.json()
    logger.info(f"response json from alphafold = {response_info}")
    # TODO: make this a unique exception with a better error message
    if not response_info:
        raise Exception
    if len(response_info) == 0:
        raise Exception

    # extract the pdb url from the response and get the pdb from the url
    # TODO: figure out why we get a list back
    pdb_url = response_info[0].get("pdbUrl")
    # TODO: make this a unique exception with a better error message
    if not pdb_url:
        raise Exception
    pdb_response = requests.get(pdb_url)

    # get the pwd to prepend to relative file paths to be made
    cwd = os.getcwd()

    # create the temp pdb file from the string input
    temp_pdb = open(os.path.join(cwd, "temp_input.pdb"), "w")
    temp_pdb.write(pdb_response.text)
    logger.info(f"writing files to {temp_pdb.name} directory")

    return temp_pdb.name


@delayed
def pdb2gmx(
    pdb_file: str,
    gro_name="conf",
    top_name="topol",
    itp_name="posre",
) -> ProteinInput:
    """
    Function that invokes the pdb2gmx helper function in Gromacs. This will
    convert a PDB file format containing the structure of a protein (say pulled
    from alphafold) to a gro and topology file for use in gromacs simulations.

    Parameters:
        pdb_file (str): contents of the PDB file in string form.
        gro_name (str): the name of the .gro file that will be created
        top_name (str): the name of the topology file to be created
        itp_name (str): the name of the .itp file to be created
    Returns:
        OutputDataProxy: output object of gmx cli api from python.
    """

    # get the pwd to prepend to relative file paths to be made
    cwd = os.getcwd()

    # build the gro file name, itp, and top file name
    if len(gro_name) < 4 or gro_name[-4:] != ".gro":
        gro_file = gro_name + ".gro"
    else:
        gro_file = gro_name
    if len(top_name) < 4 or top_name[-4:] != ".top":
        top_file = top_name + ".top"
    else:
        top_file = top_name
    if len(itp_name) < 4 or itp_name[-4:] != ".itp":
        itp_file = itp_name + ".itp"
    else:
        itp_file = itp_name
    # add directory prefix to everybody
    gro_file = os.path.join(cwd, gro_file)
    top_file = os.path.join(cwd, top_file)
    itp_file = os.path.join(cwd, itp_file)
    logger.info(f"going to store the pdb2gmx outputs at {gro_file} and {top_file}")

    args = ["pdb2gmx", "-ff", "amber03", "-water", "tip3p"]
    input_files = {"-f": pdb_file}
    output_files = {"-p": top_file, "-i": itp_file, "-o": gro_file}
    make_top = gmxapi.commandline_operation("gmx", args, input_files, output_files)

    # return the resolved file paths
    return ProteinInput.from_pdb2gmx(make_top)


# NOTE: since the gromacs output structs are heavily mixed with C++ types,
#       there's no way to provide effective type annotations here :(
@delayed
def hydrate_simulation_box(protein_gro: ProteinInput) -> ProteinInput:
    """
    Goal of this step is to hydrate the simulation box with an appropriate
    amount of water molecules. A concentration may be provided to ensure the
    simulation box is a certain size to get the right protein conc.

    Parameters:
        protein_gro (Any): The output data structure representing the protein.
    Returns:
        solvated_gro (Any): The solvated output data structure of the protein.
    """
    # get pwd to prepend to filepaths below
    cwd = os.getcwd()

    # increase the size of the protein bounding box and center it
    edit_args = ["editconf", "-c", "-d", "1.5"]
    edit_inputs = {"-f": protein_gro.gro_file}
    edit_outputs = {"-o": os.path.join(cwd, "empty_protein.gro")}
    make_edits = gmxapi.commandline_operation(
        "gmx", edit_args, edit_inputs, edit_outputs
    )

    # solvate the protein structure file with water molecules
    solv_args = ["solvate"]
    solv_inputs = {"-cs": "spc216", "-cp": make_edits.output.file["-o"].result()}
    solv_outputs = {
        "-o": os.path.join(cwd, "solvated_protein.gro"),
        "-p": os.path.join(cwd, "topol.top"),
    }
    solv_process = gmxapi.commandline_operation(
        "gmx", solv_args, solv_inputs, solv_outputs
    )

    # generate enough ions to neutralize the system
    mdp_file = os.path.join(os.path.dirname(md_inputs.__file__), "ions.mdp")
    logger.info(f"found mdp_file here: {mdp_file}")
    grompp_input_files = {
        "-f": mdp_file,
        "-c": solv_process.output.file["-o"],
        "-p": solv_process.output.file["-p"],
    }

    if not os.path.isfile(mdp_file):
        logger.error("mdp file is not found!")
        raise Exception
    grompp = gmxapi.commandline_operation(
        "gmx",
        ["grompp"],
        input_files=grompp_input_files,
        output_files={"-o": os.path.join(cwd, "ions.tpr")},
    )
    logger.info(f"grompp output = {grompp.output}")
    tpr_input_file = grompp.output.file["-o"].result()

    # call the genion command
    genion = gmxapi.commandline_operation(
        "gmx",
        ["genion", "-neutral"],
        input_files={"-s": tpr_input_file},
        output_files={
            "-o": os.path.join(cwd, "neutral.gro"),
            "-p": os.path.join(cwd, "topol.top"),
        },
        # this is hopefully always the solvent group, i wish i could reference
        # by name
        stdin="13",
    )

    # return just the output construct b/c i'm not sure if this lazy evaluates
    # it's nice to keep it lazy until it needs to happen
    return ProteinInput.from_genion(genion)


@delayed
def optimize_configuration(solv_output: ProteinInput) -> MDRun:
    """
    This step performs a steepest descent optimization to the local minimum of
    the molecular potential. The goal is to remove Any anomalously large forces
    that would cause the equilibration steps to explode.
    """

    # start by getting the tpr file using grompp
    # need to find the nvt equilibration mdp file first
    mdp_file = get_mdp_path("steep.mdp")
    tpr_file = md_grompp(mdp_file=mdp_file, input=solv_output, tpr_file_name="em.tpr")

    # read the tpr file into an input
    return standard_md_run(tpr_file, file_prefix="em")


@delayed
def md_temp_equilibrate(input: MDRun | MDRunInput, settings="nvt_eq.mdp") -> MDRun:
    """
    Equilibrate the recently minimized configuration with a short NVT MD
    simulation.
    """

    # start by getting the tpr file using grompp
    # need to find the nvt equilibration mdp file first
    mdp_file = get_mdp_path(settings)
    tpr_name = "nvt_eq.tpr"
    tpr_file = md_grompp(
        input=input,
        mdp_file=mdp_file,
        tpr_file_name=tpr_name,
        posres=True,
    )

    # read the tpr file into an input
    return standard_md_run(tpr_file, file_prefix=tpr_name.split(".")[0])


@delayed
def md_pressure_equilibrate(input: MDRun | MDRunInput, settings="npt_eq.mdp") -> MDRun:
    # def md_pressure_equilibrate(nvt_conf: str, top_file: str) -> MDRun:
    # def md_pressure_equilibrate(nvt_conf: str, top_file: str) -> tuple[Any, str, str]:
    """
    Equilibrate the recently temperature equilibrated configuration with a
    short NPT MD simulation.
    """
    # start by getting the tpr file using grompp
    # need to find the nvt equilibration mdp file first
    mdp_file = get_mdp_path(settings)
    tpr_name = "npt_eq.tpr"
    tpr_file = md_grompp(
        input=input,
        mdp_file=mdp_file,
        tpr_file_name=tpr_name,
        posres=True,
    )

    # read the tpr file into an input
    return standard_md_run(tpr_file, file_prefix=tpr_name.split(".")[0])


@delayed
def md_run(input: MDRun, settings="prod.mdp", nsteps: int | None = None) -> MDRun:
    # def md_run(npt_conf: str, top_file: str, nsteps: int | None = None) -> tuple[Any, str, str]:
    """
    Function that runs a molecular dynamics simulation for a given set of input
    files.
    """
    # start by getting the tpr file using grompp
    # need to find the nvt equilibration mdp file first
    mdp_file = get_mdp_path(settings)
    tpr_name = settings.split(".")[0] + ".tpr"
    tpr_file = md_grompp(
        input=input, mdp_file=mdp_file, tpr_file_name=tpr_name, posres=False
    )
    logger.info("editted params; ready to run MD simulation...")

    # read the tpr file into an input
    return standard_md_run(tpr_file, file_prefix="prod", nsteps=nsteps)


# -------------------------------------------------------------------
# ------------------------ Helper functions -------------------------
# -------------------------------------------------------------------
def md_grompp(
    input: ProteinInput | MDRun | MDRunInput,
    mdp_file: str,
    tpr_file_name: str = "topol.tpr",
    posres: bool = False,
) -> MDRunInput:
    """
    Helper function that runs gmx grompp in a basic way on the standard
    set of inputs, and returns the tpr output
    """
    if isinstance(input, MDRun):
        top = input.md_input.top_file
    elif isinstance(input, MDRunInput):
        top = input.top_file
    elif isinstance(input, ProteinInput):
        top = input.top_file
    else:
        raise Exception("Unknown input type.")
    grompp_input_files = {"-f": mdp_file, "-c": input.gro_file, "-p": top}
    if posres:
        grompp_input_files["-r"] = input.gro_file

    if not os.path.isfile(mdp_file):
        logger.error("mdp file is not found!")
        raise Exception
    grompp = gmxapi.commandline_operation(
        "gmx",
        ["grompp"],
        input_files=grompp_input_files,
        output_files={"-o": tpr_file_name},
    )
    logger.info(f"grompp output = {grompp.output.file.result()}")
    return MDRunInput.from_grompp(grompp_input_files, grompp)


def get_mdrun_out_files(
    md: Any,
    gro_file_name: str,
    energy_file_name: str,
) -> tuple[str, str]:
    """
    Helper that finds the full paths to the output configuration and energy
    files from a molecular dynamics run.
    """
    # get the output gro file
    logger.info(f"md dir object = {md.output.directory}")
    md_dir = md.output.directory.result()
    out_gro = os.path.join(md_dir, gro_file_name)
    out_energy = os.path.join(md_dir, energy_file_name)

    return (out_gro, out_energy)


def get_mdp_path(file_name: str) -> str:
    """
    Helper that finds the full path to the mdp file specified by the input.
    """
    return os.path.join(os.path.dirname(md_inputs.__file__), file_name)


def standard_md_run(
    md_input: MDRunInput,
    file_prefix: str = "run",
    nsteps: int | None = 10_000,
) -> MDRun:
    """
    Runs a standard MD run given a tpr, and returns the md object along with
    output gro and edr files.
    """

    tpr_input = gmxapi.read_tpr(md_input.tpr_file)
    if nsteps:
        tpr_input = gmxapi.modify_input(tpr_input, parameters={"nsteps": nsteps})
    logger.info(f"tpr input = {tpr_input}")
    gro_output = file_prefix + ".gro"
    edr_output = file_prefix + ".edr"
    traj_output = file_prefix + ".trr"
    rargs = {
        "-o": traj_output,
        "-e": edr_output,
        "-c": gro_output,
    }
    md = gmxapi.mdrun(input=tpr_input, runtime_args=rargs)
    md.run()

    # get the output gro file and edr file
    gro_path, edr_path = get_mdrun_out_files(
        md,
        gro_file_name=gro_output,
        energy_file_name=edr_output,
    )
    output = MDRun(
        md_input=md_input,
        gro_file=gro_path,
        energy=edr_path,
        md_object=md,
        trajectory=md.output.trajectory.result(),
    )

    return output
