from md_flow.flow import get_alphafold_pdb, pdb2gmx
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
