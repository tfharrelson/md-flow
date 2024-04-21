import logging
from dask.distributed import Client
from dask.delayed import Delayed
from dataclasses import dataclass
from md_flow.steps import (
    get_alphafold_pdb,
    pdb2gmx,
    hydrate_simulation_box,
    optimize_configuration,
    md_temp_equilibrate,
    md_pressure_equilibrate,
    md_run,
)
from md_flow.models import MDRun


logger = logging.getLogger(__file__)


@dataclass
class MDCluster:
    threads_per_worker: int
    n_workers: int

    def __post_init__(self):
        self.client = Client(
            threads_per_worker=self.threads_per_worker, n_workers=self.n_workers
        )

    def optimize_structure(self, uniprot_id: str) -> MDRun:
        """
        Run the optimize structure flow.
        """
        return self.client.compute(structure_opt_flow(uniprot_id))

    def npt(self, uniprot_id: str) -> MDRun:
        return self.client.compute(npt_md_flow(uniprot_id))

    def run_flow(self, flow: Delayed) -> MDRun:
        """
        Run a custom molecular dynamics flow
        """
        return self.client.compute(flow)


def structure_opt_flow(uniprot_id: str) -> Delayed:
    protein_id = get_alphafold_pdb(uniprot_id)
    protein = pdb2gmx(protein_id)
    hydrated_protein = hydrate_simulation_box(protein)
    return optimize_configuration(hydrated_protein)


def npt_md_flow(uniprot_id: str) -> Delayed:
    opt_struct = structure_opt_flow(uniprot_id)
    t_equil = md_temp_equilibrate(opt_struct)
    p_equil = md_pressure_equilibrate(t_equil)
    return md_run(p_equil)
