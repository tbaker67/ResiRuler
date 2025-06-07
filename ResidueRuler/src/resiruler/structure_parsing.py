from Bio.PDB import MMCIFParser
import numpy as np

def load_structure(file_path):
    """
    Load the structure from CIF file
    """
    if file_path.endswith('.cif'):
        parser = MMCIFParser(QUIET=True)
    else:
        raise ValueError("Unsupported file format. Only CIF files are supported.")
    structure = parser.get_structure('structure', file_path)
    return structure

def extract_CA_coords(structure):
    """
    Get the corrdinates of CA atoms and create a mapping from the chain and residue number to where they are stored in the list
    """

    #Error handling
    if not structure: 
        raise ValueError("Structure is None. Please load a valid structure.")
    if len(structure) == 0:
        raise ValueError("Structure is empty. Please load a valid structure.")
    if len(structure) > 1:
        raise ValueError("Multiple models found in the structure. Please provide a single model.")
    

    coords_list = []
    index_map = {}
    for atom in structure.get_atoms():
        if atom.get_name() == 'CA':
            res_id = atom.get_parent().get_id()
            resnum = res_id[1]
            chain_id = atom.get_parent().get_parent().get_id()
            coords_list.append(atom.get_coord())
            index_map[(chain_id, resnum)] = len(coords_list) - 1
    coords = np.array(coords_list)
    return index_map, coords


def get_coords_from_id(index_map, coords, chain, residue):
    """
    Get the coordinates for a given chain and residue from the index map.
    """
    try:
        return coords[index_map[(chain, int(residue))]].tolist()
    except KeyError:
        print(f"[WARNING] Missing coordinate for chain {chain}, residue {residue}")
        return None
    except Exception as e:
        print(f"[ERROR] Unexpected error for chain {chain}, residue {residue}: {e}")
        return None

