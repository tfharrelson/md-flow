from md_flow.flow import (
    get_alphafold_pdb, pdb2gmx, hydrate_simulation_box, optimize_configuration,
    md_temp_equilibrate, md_pressure_equilibrate, md_run
)
import logging
import pytest
import os
import tempfile


@pytest.fixture
def logger() -> logging.Logger:
    return logging.getLogger("gunicorn.error")


@pytest.fixture
def setup_pdb():
    with tempfile.TemporaryDirectory() as dir:
        cwd = os.getcwd()
        # check that we are in tmp_tests
        folder = os.path.split(cwd)[-1]

        # if we're not in the right folder raise an exception
        if folder != "tmp_tests":
            raise Exception
        os.chdir(dir)
        yield "P00250"

        # cleanup by going back to previous dir
        os.chdir(cwd)
    """
    cwd = os.getcwd()
    # check that we are in tmp_tests
    folder = os.path.split(cwd)[-1]

    # if we're not in the right folder raise an exception
    if folder != "tmp_tests":
        raise Exception
    # os.chdir(dir)
    yield "P00250"

    # cleanup by going back to previous dir
    os.chdir(cwd)
    """


def test_get_pdb(caplog, logger, setup_pdb):
    caplog.set_level(logging.INFO)
    dir = get_alphafold_pdb(setup_pdb, logger)

    assert isinstance(dir, str)


def test_convert_to_gmx(caplog, logger, setup_pdb):
    caplog.set_level(logging.INFO)

    pdb_string = get_alphafold_pdb(setup_pdb, logger)
    output = pdb2gmx(pdb_string, logger)
    gro_file = output.file['-o'].result()
    top_file = output.file['-p'].result()
    assert gro_file[-3:] == "gro"
    assert top_file[-3:] == "top"

    # check the size of the file to make sure something was made?
    assert os.path.isfile(gro_file)
    assert os.path.isfile(top_file)


def test_hydrate_box(caplog, logger, setup_pdb):
    caplog.set_level(logging.INFO)

    # TODO: get these dumb files from alphafold so i don't have to keep
    #       peppering their system for every test
    pdb_string = get_alphafold_pdb(setup_pdb, logger)
    output = pdb2gmx(pdb_string, logger)

    solvated_output = hydrate_simulation_box(output, logger)
    solvated_file = solvated_output.file["-o"].result()
    run_gro_tests(solvated_file, logger)


def run_gro_tests(file, logger):
    assert os.path.isfile(file)

    # assert that the SOL tag is in the name of the second-to-last line
    # this is a good test b/c solvent mols are added at the end of a gro
    # and the last line contains the box info
    lines = open(file, "r").readlines()
    logger.info(f"2nd to last line: {lines[-2]}")
    assert "NA" in lines[-2].split()[0]
    assert "SOL" in lines[-200].split()[0]

    # TODO: make an assertion to check for the correct bounding box size
    box_info = [float(token.strip()) for token in lines[-1].split()]
    assert all(b > 6. for b in box_info)


def test_optimize(caplog, logger, setup_pdb):
    caplog.set_level(logging.INFO)

    pdb_string = get_alphafold_pdb(setup_pdb, logger)
    output = pdb2gmx(pdb_string, logger)

    solvated_output = hydrate_simulation_box(output, logger)
    solvated_file = solvated_output.file["-o"].result()
    assert os.path.isfile(solvated_file)

    # from the solvated output run the optimization
    mdrun, min_gro, _ = optimize_configuration(solvated_output, output, logger)
    logger.info(f"mdrun = {mdrun}")

    # extract the traj file and check that it exists
    traj_file = mdrun.output.trajectory.result()
    logger.info(f"mdrun = {traj_file}")
    assert os.path.exists(traj_file)

    # check the gro file for bugs
    logger.info(f"opt gro file = {min_gro}")
    run_gro_tests(min_gro, logger)


def test_temp_eq(caplog, logger, setup_pdb):
    caplog.set_level(logging.INFO)

    pdb_string = get_alphafold_pdb(setup_pdb, logger)
    output = pdb2gmx(pdb_string, logger)

    solvated_output = hydrate_simulation_box(output, logger)
    solvated_file = solvated_output.file["-o"].result()
    assert os.path.isfile(solvated_file)

    mdrun, min_gro, _ = optimize_configuration(solvated_output, output, logger)
    # from the solvated output run the optimization
    conf_file = min_gro
    top_file = solvated_output.file['-p'].result()

    # run the sim
    mdrun, eq_gro, eq_edr = md_temp_equilibrate(conf_file, top_file, logger)
    logger.info(f"mdrun = {mdrun}")

    # extract the trajectory, make sure it exists
    traj_file = mdrun.output.trajectory.result()
    logger.info(f"mdrun = {traj_file}")

    assert os.path.exists(traj_file)

    # check the gro file for bugs
    logger.info(f"opt gro file = {eq_gro}")
    run_gro_tests(eq_gro, logger)


def test_pressure_eq(caplog, logger, setup_pdb):
    caplog.set_level(logging.INFO)

    pdb_string = get_alphafold_pdb(setup_pdb, logger)
    output = pdb2gmx(pdb_string, logger)

    solvated_output = hydrate_simulation_box(output, logger)
    solvated_file = solvated_output.file["-o"].result()
    assert os.path.isfile(solvated_file)

    # from the solvated output run the optimization
    # from the solvated output run the optimization
    conf_file = os.path.join(os.path.dirname(__file__), "nvt_eq.gro")
    top_file = os.path.join(os.path.dirname(__file__), "topol.top")
    top_file = solvated_output.file['-p'].result()

    # run the sim
    mdrun, eq_gro, eq_edr = md_pressure_equilibrate(conf_file, top_file, logger)
    logger.info(f"mdrun = {mdrun}")

    # extract the trajectory, make sure it exists
    traj_file = mdrun.output.trajectory.result()
    logger.info(f"mdrun = {traj_file}")

    assert os.path.exists(traj_file)

    # check the gro file for bugs
    logger.info(f"opt gro file = {eq_gro}")
    run_gro_tests(eq_gro, logger)


def test_prod_run(caplog, logger, setup_pdb):
    caplog.set_level(logging.INFO)

    pdb_string = get_alphafold_pdb(setup_pdb, logger)
    output = pdb2gmx(pdb_string, logger)

    solvated_output = hydrate_simulation_box(output, logger)
    solvated_file = solvated_output.file["-o"].result()
    assert os.path.isfile(solvated_file)

    # from the solvated output run the optimization
    conf_file = os.path.join(os.path.dirname(__file__), "npt_eq.gro")
    top_file = os.path.join(os.path.dirname(__file__), "topol.top")

    # run the sim
    mdrun, eq_gro, eq_edr = md_run(conf_file, top_file, logger, nsteps=10000)
    logger.info(f"mdrun = {mdrun}")

    # extract the trajectory, make sure it exists
    traj_file = mdrun.output.trajectory.result()
    logger.info(f"mdrun = {traj_file}")

    assert os.path.exists(traj_file)

    # check the gro file for bugs
    logger.info(f"opt gro file = {eq_gro}")
    run_gro_tests(eq_gro, logger)
