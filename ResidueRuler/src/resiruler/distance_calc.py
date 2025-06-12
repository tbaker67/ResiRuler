import pandas as pd
import numpy as np
from .structure_parsing import get_coords_from_id, extract_CA_coords


def get_header_indices(df):
    """
    Find where the chain headers are in the DataFrae
    """
    header_indices = []
    first_col = df.iloc[:, 0].astype(str)

    #Look for non-numeric values in the first column
    is_header = ~first_col.str.isnumeric()
    header_indices = np.where(is_header)[0].tolist()
    header_indices.append(len(df))
    return header_indices

def compute_distance(coord1, coord2):
    if coord1 is None or coord2 is None:
        return np.nan
    #print(f"[DEBUG] Computing distance between {coord1} and {coord2}")
    return np.linalg.norm(np.array(coord1) - np.array(coord2))

def read_data(df, header_indices, chain_mapping, index_map, coords):
    """
    Read the data from the DataFrame and compute distances.
    """
    output_rows = []
    chain_1 = None
    chain_2 = None

    #Go through all headerblocks
    for start, end in zip(header_indices[:-1], header_indices[1:]):
        header_row = df.iloc[start].dropna().tolist()
        if header_row[0] not in chain_mapping or header_row[1] not in chain_mapping:
            print(f"[WARNING] Chain {header_row[0]} or {header_row[1]} not found in mapping.")
            continue

        #Go through possible combos of chains
        for chain1 in chain_mapping[header_row[0]]:
            for chain2 in chain_mapping[header_row[1]]:

                block = df.iloc[start+1:end].copy()
                cols = list(block.columns)
                if len(cols) >= 2:
                    cols[0] = 'Residue1'
                    cols[1] = 'Residue2'
                block.columns = cols
                
                    
                #Add new identfiers columns
                block['Chain1_Residue1'] = chain1 + '_' + block['Residue1'].astype(int).astype(str)
                block['Chain2_Residue2'] = chain2 + '_' + block['Residue2'].astype(int).astype(str)

                #Get the coordinates
                block['Coord1'] = block.apply(
                    lambda row: get_coords_from_id(index_map, coords, chain1, row['Residue1']), axis=1
                )
                block['Coord2'] = block.apply(
                    lambda row: get_coords_from_id(index_map, coords, chain2, row['Residue2']), axis=1
                )

                #Compute the distances
                block['Distance'] = block.apply(
                    lambda row: compute_distance(row['Coord1'], row['Coord2']), axis=1
                )

                #Reformat
                block['Coord1'] = block['Coord1'].apply(lambda x: f"[{x[0]:.3f}, {x[1]:.3f}, {x[2]:.3f}]" if isinstance(x, (list, np.ndarray)) else np.nan)

                block['Coord2'] = block['Coord2'].apply(lambda x: f"[{x[0]:.3f}, {x[1]:.3f}, {x[2]:.3f}]" if isinstance(x, (list, np.ndarray)) else np.nan)

                output_rows.append(block)
    if output_rows:
        return pd.concat(output_rows, ignore_index=True).dropna(subset=['Distance']).reset_index(drop=True)
    else:
        print(f"[WARNING] No data found")
        
def calc_difference_aligned(structure1, structure2, chain_mapping=None):
    """
    Takes in two aligned structures and outputs a dataframe containing the distances between corresponding residues in the structure
    """
    data = []
    #Go through all atoms in first structure 
    for atom in structure1[0].get_atoms():
        if atom.get_name() == 'CA':
            chain_id = atom.get_parent().get_parent().get_id()
            res_id = atom.get_parent().get_id()
            resnum = res_id[1]
            coord1 = atom.get_coord()
            chain_id2=None
            coord2 = None
            try:
                if chain_mapping:
                    #Need to use correct chain in structure2
                    chain_id2 = chain_mapping[chain_id]
                    coord2 = structure2[0][chain_id2][res_id]['CA'].get_coord()
                    #print(f"[DEBUG] Match found for {chain_id}_{resnum}")
     
                else:
                    coord2 = structure2[0][chain_id][res_id]['CA'].get_coord()
            except KeyError:
                print(f"[WARNING] No matching residue found for {chain_id}_{resnum}")
                continue
            diff_vec = coord2 - coord1
            distance = np.linalg.norm(diff_vec)
            #ADd it all to dataframe
            data.append({
            'ChainID_Resnum1': f'{chain_id}_{resnum}',
            'ChainID_Resnum2': f'{chain_id2}_{resnum}',
            'Coord1': f"[{coord1[0]:.3f}, {coord1[1]:.3f}, {coord1[2]:.3f}]",
            'Coord2': f"[{coord2[0]:.3f}, {coord2[1]:.3f}, {coord2[2]:.3f}]",
            'Diff_Vec': f"[{diff_vec[0]:.3f}, {diff_vec[1]:.3f}, {diff_vec[2]:.3f}]",
            'Distance': distance
            })
    return pd.DataFrame(data).dropna()