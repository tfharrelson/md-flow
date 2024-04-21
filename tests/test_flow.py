from md_flow.steps import (
    get_alphafold_pdb,
    pdb2gmx,
    hydrate_simulation_box,
    optimize_configuration,
    md_temp_equilibrate,
    md_pressure_equilibrate,
    md_run,
    md_grompp,
    get_mdp_path,
)
from md_flow.models import ProteinInput
import logging
import pytest
import os
import tempfile


@pytest.fixture
def logger() -> logging.Logger:
    return logging.getLogger(__file__)


@pytest.fixture(scope="function")
def setup_pdb():
    with tempfile.TemporaryDirectory() as dir:
        cwd = os.getcwd()
        os.chdir(dir)
        yield "P00250"

        # cleanup by going back to previous dir
        os.chdir(cwd)


def test_get_pdb(caplog, logger, setup_pdb):
    caplog.set_level(logging.INFO)
    dir = get_alphafold_pdb(setup_pdb)
    out = dir.compute()

    assert isinstance(out, str)


def test_convert_to_gmx(caplog, setup_pdb, logger):
    caplog.set_level(logging.INFO)

    pdb_string = get_alphafold_pdb(setup_pdb)
    output = pdb2gmx(pdb_string).compute()
    logger.info(f"output = {output}")
    gro_file = output.gro_file
    top_file = output.top_file
    assert gro_file[-3:] == "gro"
    assert top_file[-3:] == "top"

    # check the size of the file to make sure something was made?
    assert os.path.isfile(gro_file)
    assert os.path.isfile(top_file)


def test_hydrate_box(caplog, logger, setup_pdb):
    caplog.set_level(logging.INFO)

    # TODO: get these dumb files from alphafold so i don't have to keep
    #       peppering their system for every test
    pdb_string = get_alphafold_pdb(setup_pdb)
    output = pdb2gmx(pdb_string)

    solvated_output = hydrate_simulation_box(output)
    solvated_file = solvated_output.compute().gro_file
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
    assert all(b > 6.0 for b in box_info)


def test_optimize(caplog, logger, setup_pdb):
    caplog.set_level(logging.INFO)

    pdb_string = get_alphafold_pdb(setup_pdb)
    output = pdb2gmx(pdb_string)

    solvated_output = hydrate_simulation_box(output)

    # from the solvated output run the optimization
    mdrun = optimize_configuration(solvated_output)
    logger.info(f"mdrun = {mdrun}")
    solvated_file = solvated_output.compute().gro_file
    assert os.path.isfile(solvated_file)

    # extract the traj file and check that it exists
    md_out = mdrun.compute()
    traj_file = md_out.trajectory
    min_gro = md_out.gro_file
    logger.info(f"mdrun = {traj_file}")
    assert os.path.exists(traj_file)

    # check the gro file for bugs
    logger.info(f"opt gro file = {min_gro}")
    run_gro_tests(min_gro, logger)


def test_pressure_eq(caplog, logger, setup_pdb):
    caplog.set_level(logging.INFO)

    pdb_string = get_alphafold_pdb(setup_pdb)
    output = pdb2gmx(pdb_string)

    solvated_output = hydrate_simulation_box(output)

    opt = optimize_configuration(solvated_output)
    opt = opt.compute()
    print(f"finished opt = {opt}")
    # run the sim
    t_mdrun = md_temp_equilibrate(opt)
    t_mdrun = t_mdrun.compute()
    print(f"finished t eq = {t_mdrun}")
    mdrun = md_pressure_equilibrate(t_mdrun)
    logger.info(f"mdrun = {mdrun}")

    solvated_file = solvated_output.compute().gro_file
    assert os.path.isfile(solvated_file)

    # extract the trajectory, make sure it exists
    md_out = mdrun.compute()
    traj_file = md_out.trajectory
    logger.info(f"mdrun = {traj_file}")

    assert os.path.exists(traj_file)

    # check the gro file for bugs
    eq_gro = md_out.gro_file
    logger.info(f"opt gro file = {eq_gro}")
    run_gro_tests(eq_gro, logger)


def test_temp_eq(caplog, logger, setup_pdb):
    caplog.set_level(logging.INFO)

    pdb_string = get_alphafold_pdb(setup_pdb)
    output = pdb2gmx(pdb_string)

    solvated_output = hydrate_simulation_box(output)

    mdrun = optimize_configuration(solvated_output)

    # run the sim
    mdrun = md_temp_equilibrate(mdrun)
    logger.info(f"mdrun = {mdrun}")

    solvated_file = solvated_output.compute().gro_file
    assert os.path.isfile(solvated_file)

    # extract the trajectory, make sure it exists
    md_out = mdrun.compute()
    traj_file = md_out.trajectory
    logger.info(f"mdrun = {traj_file}")
    assert os.path.exists(traj_file)

    # check the gro file for bugs
    logger.info(f"opt gro file = {md_out.gro_file}")
    run_gro_tests(md_out.gro_file, logger)


def test_prod_run(caplog, logger, setup_pdb):
    caplog.set_level(logging.INFO)

    pdb_string = get_alphafold_pdb(setup_pdb)
    output = pdb2gmx(pdb_string)

    solvated_output = hydrate_simulation_box(output)

    # from the solvated output run the optimization
    conf_file = os.path.join(os.path.dirname(__file__), "npt_eq.gro")
    top_file = os.path.join(os.path.dirname(__file__), "topol.top")

    prot_input = ProteinInput(gro_file=conf_file, top_file=top_file)
    md_input = md_grompp(prot_input, mdp_file=get_mdp_path("prod.mdp"))

    # run the sim
    mdrun = md_run(md_input, nsteps=10000)
    logger.info(f"mdrun = {mdrun}")

    # check the solvated output
    solvated_file = solvated_output.compute().gro_file
    assert os.path.isfile(solvated_file)

    # extract the trajectory, make sure it exists
    md_out = mdrun.compute()
    traj_file = md_out.trajectory
    logger.info(f"mdrun = {traj_file}")

    assert os.path.exists(traj_file)

    # check the gro file for bugs
    eq_gro = md_out.gro_file
    logger.info(f"opt gro file = {eq_gro}")
    run_gro_tests(eq_gro, logger)
