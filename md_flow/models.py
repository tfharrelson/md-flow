from __future__ import annotations

from dataclasses import dataclass
from typing import Any
import logging


logger = logging.getLogger(__file__)


@dataclass
class ProteinInput:
    gro_file: str
    top_file: str

    @staticmethod
    def from_pdb2gmx(gmx_top: Any) -> ProteinInput:
        gro = gmx_top.output.file["-o"].result()
        top = gmx_top.output.file["-p"].result()
        return ProteinInput(gro_file=gro, top_file=top)

    @staticmethod
    def from_genion(genion: Any) -> ProteinInput:
        gro = genion.output.file["-o"].result()
        top = genion.output.file["-p"].result()
        return ProteinInput(gro, top)


@dataclass
class MDRunInput:
    tpr_file: str
    gro_file: str
    top_file: str
    settings_file: str | None = None
    itp_file: str | None = None
    nsteps: int = 10000
    grompp: Any | None = None

    @staticmethod
    def from_grompp(input_files: dict[str, str], gmp: Any) -> MDRunInput:
        logger.info(f"creating MD input from grompp input = {gmp}")
        gro = input_files["-c"]
        top = input_files["-p"]
        itp = input_files.get("-r")
        settings = input_files["-f"]
        tpr = gmp.output.file["-o"].result()
        return MDRunInput(tpr, gro, top, settings, itp, grompp=gmp)


@dataclass
class SteepInput:
    gro_file: str
    top_file: str
    idp_file: str | None = None
    settings_file: str | None = "steep.mdp"


@dataclass
class NVTEquilibrationInput:
    gro_file: str
    top_file: str
    idp_file: str | None = None
    settings_file: str | None = "nvt_eq.mdp"


@dataclass
class NPTEquilibrationInput:
    gro_file: str
    top_file: str
    idp_file: str | None = None
    settings_file: str | None = "npt_eq.mdp"


@dataclass
class MDRun:
    md_input: MDRunInput
    gro_file: str
    md_object: Any
    energy: str
    trajectory: str
