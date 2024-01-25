# from metaflow import FlowSpec
import requests
import logging
import subprocess
import os
import gmxapi

# TODO: initialize logger here for now, but encapsulate it into a class later
logger = logging.getLogger("gunicorn.error")


def get_alphafold_pdb(uniprot_id: str, logger: logging.Logger) -> str:
    """
    Helper function that takes a uniprot ID and returns a path to a PDB file.

    Parameters:
        uniprot_id (str): A string representing a unique ID of a protein
            existingin the uniprot database
        logger (logging.Logger): A logger that manages critical errors

    Returns:
        pdb_string (str): the string (stored in RAM) to the PDB file
            representing a structural instance of the protein encoded by the
            uniprot ID
    """
    url = "https://alphafold.ebi.ac.uk/api"
    route = f"/prediction/{uniprot_id}"

    response = requests.get(url+route)
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

    return pdb_response.text


def pdb2gmx(
    file_string: str,
    logger: logging.Logger,
    gro_name="conf",
    top_name="topol",
    itp_name="posre",
    temp_dir=""
) -> None:
    """
    Function that invokes the pdb2gmx helper function in Gromacs. This will
    convert a PDB file format containing the structure of a protein (say pulled
    from alphafold) to a gro and topology file for use in gromacs simulations.

    Parameters:
        file_string (str): contents of the PDB file in string form.
        logger (logging.Logger): class that sends messages to stdin/err.
        gro_name (str): the name of the .gro file that will be created
        top_name (str): the name of the topology file to be created
        temp_dir (str): the name of the directory to store the outputs

    Returns:
        Nothing
    """
    # create the temp pdb file from the string input
    temp_pdb = open(os.path.join(temp_dir, "temp_input.pdb"), 'w')
    temp_pdb.write(file_string)
    logger.info(f"writing files to {temp_pdb.name} directory")

    # build the gro file name, itp, and top file name
    if gro_name[-4:] != ".gro":
        gro_file = gro_name + ".gro"
    else:
        gro_file = gro_name
    if top_name[-4:] != ".top":
        top_file = top_name + ".top"
    else:
        top_file = top_name
    if itp_name[-4:] != ".itp":
        itp_file = itp_name + ".itp"
    else:
        itp_file = itp_name

    # add the temp dir to each of the output file names
    gro_file = os.path.join(temp_dir, gro_file)
    top_file = os.path.join(temp_dir, top_file)
    itp_file = os.path.join(temp_dir, itp_file)
    logger.info(
        f"going to store the pdb2gmx outputs at {gro_file} and {top_file}"
    )

    args = ['pdb2gmx', '-ff', 'amber03', '-water', 'tip3p']
    input_files = {'-f': temp_pdb.name}
    output_files = {
        '-p': top_file,
        '-i': itp_file,
        '-o': gro_file
    }
    make_top = gmxapi.commandline_operation(
        'gmx', args, input_files, output_files
    )

    # return the resolved file paths
    return make_top.output


def hydrate_simulation_box():
    pass


def optimize_configuration():
    pass


def md_equilibrate():
    pass


def md_run():
    """
    Function that runs a molecular dynamics simulation for a given set of input
    files.
    """
    # input =
    pass
