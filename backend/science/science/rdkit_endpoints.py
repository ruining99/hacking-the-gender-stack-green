from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from typing import List

ATOM_PROP_ATOM_LABEL = "atomLabel"


def generate_image(mol_smi: str, width: int = 400, height: int = 400) -> str:
    """
    Generates an image of an rdkit mol represented by the given smiles.

    *THIS WILL BE GIVEN TO PARTICIPANTS*

    :param mol_smi: SMILES of mol to display
    :param width: width of the image
    :param height: height of the image
    :return: generated image as an SVG string
    """
    mol = Chem.MolFromSmiles(mol_smi)
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText().encode()


def get_rgroup_labels(core_smi: str) -> List[str]:
    """
    Retrieves rlabels from core/scaffold (such as ["R1", "R2", "R3"]). Front
    end should use this when asking for R1s, R2s, etc

    :param core_smis: SMILES that represents core
    :return: List of R labels
    """
    core = Chem.MolFromSmiles(core_smi)
    core_rlabels = []
    for at in core.GetAtoms():
        if at.HasProp(ATOM_PROP_ATOM_LABEL):
            core_rlabels.append(at.GetProp(ATOM_PROP_ATOM_LABEL))
    return sorted(core_rlabels)


"""
def valid(smiles):
    return True
"""


def valid(mol_smi: str):
    mol = Chem.MolFromSmiles(mol_smi)
    if type(mol) == rdkit.Chem.rdchem.Mol:
        return True
    else:
        return False


def get_props(s):
    mol = Chem.MolFromSmiles(s)
    props = QED.properties(mol)
    props_str = str(props)
    split_props = props_str.split()

    # Hydrogen bond donors (No more than 5 hydrogen bond donors)
    HBD0 = split_props[3]
    HBD1 = HBD0.split("=")
    HBD = HBD1[1]

    # Hydrogen bond acceptors (No more than 10 hydrogen bond acceptors)
    HBA0 = split_props[2]
    HBA1 = HBA0.split("=")
    HBA = HBA1[1]

    # Molecular mass (A molecular mass less than 500)
    MW0 = split_props[0]
    MW1 = MW0.split("=")
    MW = MW1[1]

    # Octanol-water partition coefficient (does not exceed 5)
    ALOGP0 = split_props[1]
    ALOGP1 = ALOGP0.split("=")
    ALOGP = ALOGP1[1]

    return f"This molecule has {HBD[:-1]} hydrogen bond donors, {HBA[:-1]} hydrogen bond acceptors. The molecular mass is {MW[:-1]} and the ALOGP is {ALOGP[:-1]}"


def lipinski_violations(smile):
    violations = 0

    output = get_props(smile)
    output_split = output.split()
    HBD = int(output_split[3])
    HBA = int(output_split[7])
    MW = int(float(output_split[15]))
    ALOGP = int(float(output_split[20]))

    if HBD > 5:
        violations += 1
    elif HBA > 10:
        violations += 1
    elif MW > 500:
        violations += 1
    elif ALOGP > 5:
        violations += 1
    return violations


def lipinski_test_result(smile):
    violations = lipinski_violations(smile)
    if violations > 1:
        return "this molecule does not pass the lipinski test"
    else:
        return "this molecular passes the lipinski test"
