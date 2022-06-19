"""Functions for loading PDB files."""

import io
from pathlib import Path

from Bio.PDB import PDBParser
from Bio.PDB.Structure import Structure


def convert_to_path(path: str) -> Path:
    if not isinstance(path, Path):
        return Path(path)
    else:
        return path


def read_pdb_from_bytes(id: str, pdb_bytes: bytes) -> Structure:
    """
    Reads pdb file into a Structure object.

    Parameters:
    -----------
    id (str): id of structure object.
    pdb_bytes (bytes): byte object containing file data.

    Returns:
    --------
    Structure: Structure object containing data from the PDB file.
    """
    pdb_stream = io.StringIO(pdb_bytes.decode('utf-8').replace('\r', '\n'))
    return PDBParser().get_structure(id, pdb_stream)


def read_pdb_from_path(pdb_path: Path | str) -> Structure:
    """
    Reads a pdb file into a Structure object.

    Parameters:
    -----------
    pdb_path (Path | str): path to pdb file.

    Returns:
    --------
    Structure: Structure object containing data from the PDB file.
    """
    pdb_path = convert_to_path(pdb_path)
    id = pdb_path.stem
    return PDBParser().get_structure(id, pdb_path)
