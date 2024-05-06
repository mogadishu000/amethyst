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
    """Parses provided file to a Substituents dataclass. There is planned support for multiple R-groups in one file. For now please supply one R-group per file. Accepts R-groups marked as either isotope labels or atom maps. Newline separated file can only be for one R#.

    Args:
        filepath (str): Path to file with R-groups. Required.
        r_num (Optional[int]): Number of R-group attached to the Substituents dataclass. Can be passed as R# (case insensitive) in a filename.  Defaults to R1.
        delimiter (Optional[str], optional): Delimiter for one-line files. Defaults to newline.

    Raises:
        FileNotFoundError: Raised if supplied path isn't a file.
        ValueError: Raised if no r_num was provided.

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
    r = f"[{r_num}*]"

    with open(filepath, "r") as file:
        if delimiter is None:
            for i in file:
                m = re.sub(r"\[[0-9]\*\]|\[.\:[0-9]\]", r, i)
                subs_list.append(MolFromSmiles(m))
                logger.debug(f"SMILES added: {m}")
        else:
            # FIXME - maybe another flag for parsing multiple r_nums in one file?
            lines: List[str] = file.readlines()
            logger.debug(*lines)
            for line in lines:
                m = re.sub(r"\[[0-9]\*\]|\[.\:[0-9]\]", r, line)
                split_lines: List[str] = m.split(delimiter)
                [subs_list.append(MolFromSmiles(x)) for x in split_lines]
                logger.debug(f"SMILES added: {*split_lines,}")

    logger.debug(f"Final sub list: {mols_to_str(subs_list)}")

    for i in subs_list:
        RelabelMappedDummies(i, outputLabels=RGroupLabelling.AtomMap)
        logger.debug(f"Sub modified: {MolToSmiles(i)}")

    subs = Substituents(r_num, subs_list)

    return subs


def parse_mol_input(mols: List[List[Union[Mol, str]]]) -> List[Substituents]:
    """Parses list of molecules to a Substituents class. R# is handled via the list index (n+1) e.g., first list of Mol's in the list passed will have R1 number and so on.

    Args:
        mols (List[List[Mol]]): List containing another list of R-groups.

    Raises:
        ValueError: Raised when input isn't Mol or str.

    Returns:
        List[Substituents]: Returns subs parsed into a list of Substituents dataclass.
    """
    r_num = 1
    substituents_list = []
    for i in mols:
        if type(i[0]) is Mol:
            logger.debug("Mol input detected.")
            mols_smi = [MolToSmiles(x) for x in i]
        elif type(i[0]) is str:
            logger.debug("String input detected.")
            mols_smi = mols[r_num - 1]
        else:
            raise ValueError("Wrong input type.")

        smis_relabelled = []
        for j in i:
            r = f"[*:{r_num}]"
            logger.debug(f"Inner: {j}")
            if type(j) is str:
                m = re.sub(r"\[[0-9]\*\]|\[\*\:[0-9]\]", r, j)
            elif type(j) is Mol:
                smi = MolToSmiles(j)
                m = re.sub(r"\[[0-9]\*\]|\[\*\:[0-9]\]", r, smi)
            else:
                raise ValueError("Wrong input type.")
            logger.debug(f"Relabelled SMILES: {m}")
            smis_relabelled.append(m)
            logger.debug(f"SMILES relabelled: {*smis_relabelled, }")

        mols_relabelled = [MolFromSmiles(x) for x in smis_relabelled]
        substituents_list.append(Substituents(r_num, mols_relabelled))
        logger.debug(f"R{r_num} SMILES: {mols_to_str(mols_relabelled)}")

        r_num = r_num + 1

    return substituents_list
