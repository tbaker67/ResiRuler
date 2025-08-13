from Bio.PDB import MMCIFParser
from Bio.SeqUtils import seq1
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

def check_input_structure(structure):
    if not structure: 
        raise ValueError("Structure is None. Please load a valid structure.")
    if len(structure) == 0:
        raise ValueError("Structure is empty. Please load a valid structure.")
    if len(structure) > 1:
        raise ValueError("Multiple models found in the structure. Please provide a single model.")
    

def extract_res_from_chain(chain):
    res_list = [res for res in chain.get_residues() if res.id[0] == ' ']
    return res_list

def extract_seq_from_chain(chain):
    """
    Extracts sequence and residue list from a Biopython chain,
    including support for insertion codes and unknown residues.
    """
    seq = ""

    for res in chain.get_residues():
        if res.id[0] != ' ':  # skip heteroatoms
            continue
        try:
            seq += seq1(res.get_resname())
        except KeyError:
            seq += "X"
    
    return seq


def extract_residues_from_structure(structure, chains=None):
    """
    Get a list of all residues in a strucutre.
    Optionally use a subset of chains from which to pull residues from
    """

    #Error handling
    check_input_structure(structure)
    if chains is None:
        residues = [res for res in structure.get_residues() if res.id[0] == ' ']
        return residues
    else:
        residues = [res for res in structure.get_residues() if res.id[0] == ' ' and res.get_parent().get_id() in chains]
        return residues

def extract_residues_by_type(structure, residue_type, chains):
    '''
    Extract a specific kind of residue type (ie LYS, CYS)
    Make sure to use the three letter code as the type
    '''
    #Error handling
    check_input_structure(structure)

    if residue_type not in set("ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
                               "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR",
                               "TRP", "TYR", "VAL"):
        raise ValueError("Unknown Residue Type. Please enter a valid 3 letter type (ie LYS)")
    
    if chains is None:
        residues = [res for res in structure.get_residues() if res.id[0] == ' ' and res.get_resname() == residue_type]
        return residues
    
    else:
        residues = [res for res in structure.get_residues() if res.id[0] == ' ' and res.get_parent().get_id() in chains
                                                                                and res.get_resname() == residue_type]
        return residues
    

def convert_to_CA_coords_list(res_list):
    """
    Convert a list of residues into a list full of the coordinates for each CA atom
    Also create an index map the maps (Chain ID, Resnum) -> index in corrdinate list
    """
    coords = []
    id_to_indices = {}
    index = 0
    for res in res_list:
        key = (res.get_parent().get_id(), res.get_id())
        coords.append(res["CA"].get_coord())
        id_to_indices[key] = index
        index += 1
    return coords, id_to_indices



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

