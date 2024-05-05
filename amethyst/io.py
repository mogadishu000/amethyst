import os
import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Union

from loguru import logger
from rdkit.Chem.AllChem import Mol, MolFromSmiles, MolToSmiles
from rdkit.Chem.rdRGroupDecomposition import RelabelMappedDummies, RGroupLabelling

from amethyst.utils import mols_to_str

logger.add(sink="io.log", level=10)


@dataclass
class Substituents:
    """Dataclass used for storing a list of R-groups with their respective R#.

    Attributes:
        r_num (int): R#
        subs (List[Mol]): List of substituents in Mol object

    """

    r_num: int
    subs: List[Mol]


# FIXME - Set the atom map to appropriate R#
def parse_file_input(
    filepath: str, r_num: int, delimiter: Optional[str] = None
) -> Substituents:
    """Parses provided file to a Substituents dataclass. There is planned support for multiple R-groups in one file. For now please supply one R-group per file.

    Args:
        filepath (str): Path to file with R-groups. Required.
        delimiter (Optional[str], optional): Delimiter for one-line files. Defaults to newline.
        r_num (Optional[int], optional): Number of R-group attached to the Substituents dataclass. Can be passed as R# (case insensitive) in a filename.  Defaults to R1.

    Raises:
        FileNotFoundError: Raised if supplied path isn't a file.

    Returns:
        Substituents: Parsed file saved to dataclass.
    """
    if os.path.isfile(filepath):
        logger.debug("File path is good.")
        pass
    else:
        logger.error(f"{filepath} is not a file!")
        raise FileNotFoundError(f"{filepath} is not a file!")

    if r_num is None:
        path_split = re.split(r"(\\\\)|(/)|(\\)", filepath)
        m = re.match("[rR][0-9]", path_split[-1])
        if m is not None:
            r_num = m.group(0)
        else:
            raise ValueError("R# is missing.")

    logger.debug(f"File given for R{r_num}.")

    subs_list: List[Mol] = []

    with open(filepath, "r") as file:
        if delimiter is None:
            for i in file:
                subs_list.append(MolFromSmiles(i))
                logger.debug(f"SMILES added: {i}")
        else:
            # FIXME - change to readlines(). maybe another flag for parsing multiple r_nums in one file?
            line: str = file.readline()
            logger.debug(line)
            split_lines: List[str] = line.split(delimiter)
            [subs_list.append(MolFromSmiles(x)) for x in split_lines]
            logger.debug(f"SMILES added: {*split_lines,}")

    logger.debug(f"Final sub list: {mols_to_str(subs_list)}")

    for i in subs_list:
        RelabelMappedDummies(i, outputLabels=RGroupLabelling.AtomMap)
        logger.debug(f"Sub modified: {MolToSmiles(i)}")

    subs = Substituents(r_num, subs_list)

    return subs


def parse_string_input():
    pass


def parse_mol_input(mols: List[List[Mol]]) -> List[Substituents]:
    """Parses list of molecules to a Substituents class. R# is handled via the list index (n+1) e.g., first list of Mol's in the list passed will have R1 number and so on.

    Args:
        mols (List[List[Mol]]): List containing another list of R-groups.

    Returns:
        List[Substituents]: Returns subs parsed into a list of Substituents dataclass.
    """
    r_num = 1
    substituents_list = []
    for i in mols:
        substituents_list.append(Substituents(r_num, i))
        r_num = r_num + 1

    return substituents_list
