import pandas as pd
import numpy as np
from .structure_parsing import get_coords_from_id


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

                output_rows.append(block)
    if output_rows:
        return pd.concat(output_rows, ignore_index=True).dropna(subset=['Distance']).reset_index(drop=True)
    else:
        print(f"[WARNING] No data found")
        
