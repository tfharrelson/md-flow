from md_flow.flow import get_alphafold_pdb, pdb2gmx, hydrate_simulation_box
import logging
import pytest
import os
import tempfile


@pytest.fixture
def logger() -> logging.Logger:
    return logging.getLogger("gunicorn.error")


def test_get_pdb(caplog, logger):
    caplog.set_level(logging.INFO)
    dir = get_alphafold_pdb("P00250", logger)

    assert isinstance(dir, str)


def test_convert_to_gmx(caplog, logger):
    caplog.set_level(logging.INFO)

    pdb_string = get_alphafold_pdb("P00250", logger)
    with tempfile.TemporaryDirectory() as temp_dir:
        output = pdb2gmx(pdb_string, logger, temp_dir=temp_dir)
        gro_file = output.file['-o'].result()
        top_file = output.file['-p'].result()
        assert gro_file[-3:] == "gro"
        assert top_file[-3:] == "top"

        # check the size of the file to make sure something was made?
        assert os.path.isfile(gro_file)
        assert os.path.isfile(top_file)


def test_hydrate_box(caplog, logger):
    caplog.set_level(logging.INFO)

    # TODO: get these dumb files from alphafold so i don't have to keep
    #       peppering their system for every test
    pdb_string = get_alphafold_pdb("P00250", logger)
    with tempfile.TemporaryDirectory() as temp_dir:
        output = pdb2gmx(pdb_string, logger, temp_dir=temp_dir)

        solvated_output = hydrate_simulation_box(output, logger, temp_dir)
        solvated_file = solvated_output.file["-o"].result()
        logger.info(f"the solvated output is {"".join(open(solvated_file, "r").readlines())}")
        assert os.path.isfile(solvated_file)

        # assert that the SOL tag is in the name of the second-to-last line
        # this is a good test b/c solvent mols are added at the end of a gro
        # and the last line contains the box info
        lines = open(solvated_file, "r").readlines()
        logger.info(f"2nd to last line: {lines[-2]}")
        assert "SOL" in lines[-2].split()[0]

        # TODO: make an assertion to check for the correct bounding box size
        box_info = [float(token.strip()) for token in lines[-1].split()]
        assert all(b > 6. for b in box_info)
