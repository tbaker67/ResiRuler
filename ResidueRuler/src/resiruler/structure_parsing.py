from Bio.PDB import MMCIFParser
from Bio.SeqUtils import seq1
import numpy as np

#This dictionary only works for cif files that are properly labeled
#ie DNA residue names will always begin with D
RESIDUE_TYPE = {
    #Standard amino acids
    "ALA": "protein", "ARG": "protein", "ASN": "protein", "ASP": "protein",
    "CYS": "protein", "GLN": "protein", "GLU": "protein", "GLY": "protein",
    "HIS": "protein", "ILE": "protein", "LEU": "protein", "LYS": "protein",
    "MET": "protein", "PHE": "protein", "PRO": "protein", "SER": "protein",
    "THR": "protein", "TRP": "protein", "TYR": "protein", "VAL": "protein",
    #Non-standard amino acids
    "MSE": "protein",  
    "SEC": "protein",  
    "PYL": "protein",  
   
    "DA": "dna", "DC": "dna", "DG": "dna", "DT": "dna",
    
    
    "A": "rna", "C": "rna", "G": "rna", "U": "rna",
    
    #Common RNA modifications
    "I": "rna",     
    "PSU": "rna",   
    "PSE": "rna",  
    "OMC": "rna",  
    "OMU": "rna",   
    "M2G": "rna",  
    "1MA": "rna",   
    "2MG": "rna",  
    "5MC": "rna",   
    "5MU": "rna",  
}

RESIDUE_TO_ONE_LETTER = {
    # --- Protein ---
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C",
    "GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
    "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
    "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
    # Common protein mods
    "MSE":"M","SEC":"C","PYL":"K",

    # --- DNA ---
    "DA":"A","DC":"C","DG":"G","DT":"T",
    # --- RNA ---
    "A":"A","C":"C","G":"G","U":"U","I":"I",
    # RNA modifications
    "PSU":"U","PSE":"U","OMC":"C","OMU":"U","1MA":"A",
    "5MC":"C","5MU":"U","M2G":"G","7MG":"G","H2U":"U","OMG":"G",
}

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
    Extracts sequence from a Biopython chain,
    """
    seq = ""

    for res in chain.get_residues():
        if res.id[0] != ' ':  # skip heteroatoms
            continue
        try:
            seq += RESIDUE_TO_ONE_LETTER[res.get_resname()]
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
    
def get_CA_from_residue(res):
    """
    Assumes the input residue is a biopython residue object
    """
    
    if "CA" in res:
        return res["CA"].get_coord()

    return None

def get_CB_from_residue(res):
    """
    Assumes the input residue is a biopython residue object
    """
    
    if "CB" in res:
        return res["CB"].get_coord()
    elif res.get_resname() == "GLY" and "CA" in res:
        return res["CA"].get_coord()
    
    return None

def get_SC_from_residue(res):
    """
    gets sidechain coordinate for a given residue based on the centroid of the heavy carbon atoms (not including CA)
    returns None if it can't find a coordinate
    Assumes the input residue is a biopython residue object
    """
    BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}

    for atom in res:

        sidechain_atom_coords = [atom.get_coord() for atom in res if atom.get_id() not in BACKBONE_ATOMS]

    if sidechain_atom_coords:
        return np.mean(sidechain_atom_coords, axis=0)
    
    ##TODO: Think about how we want to handle this
    elif res.get_resname() == "GLY" and "CA" in res:
        return res["CA"].get_coord()
    return None

def get_C1prime_from_residue(res):
    """
    Use to get the C1' coordinate for a given residue
    Assumes input is a valid biopython residue object
    """

    if "C1'" in res:
        return res["C1'"].get_coord()
    print(f"[WARNING] no C1' found for {res.get_parent().id}-{res.id[1]}")
    return None 

class ChainInfo:
    def __init__(self, chain):
        self.chain = chain
        self.seq = extract_seq_from_chain(chain)
        self.res_list = extract_res_from_chain(chain)
        self.type = self.detect_chain_type(self.res_list)
    
    def detect_chain_type(self, res_list):
        res_types = set()
        for res in res_list:
            res_type = RESIDUE_TYPE.get(res.resname, None)

            if res_type is None:
                print(f"[WARNING] Unknown residue '{res.resname}' in chain {self.chain.id}")
                continue
            res_types.add(res_type)
        
        if not res_types:
            print(f"[Warning] chain {self.chain.id} contains no known residues")
            return "unknown"

        elif len(res_types) > 1:
            print(f"[Warning] chain {self.chain.id} contains a mix of RNA/DNA/PROTEIN residues")
            return "unknown"

        else:
            return next(iter(res_types))

            
class ChainCollection:
    def __init__(self, structure):
        self.chains = {}
        for chain in structure.get_chains():
            self.chains[chain.id] = ChainInfo(chain)

    def valid_pairs(self, other_collection):
        """
        Creates a generate to iterate through valid pairings (matching chain types) one-by-one
        Intended for use with auto-aligning when a structure has a mix of DNA, RNA, and Protein chains
        """
        for ref_id, ref_info in self.chains.items():
            for tgt_id, tgt_info in other_collection.chains.items():
                if ref_info.type != tgt_info.type:
                    continue
                yield ref_id, tgt_id, ref_info, tgt_info
        
