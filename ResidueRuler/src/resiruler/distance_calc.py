import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist
from .structure_parsing import get_coords_from_id

class DistanceMatrix:
    def __init__(self, coords_list, index_map):
        self.coords = np.array(coords_list)
        self.index_map = dict(sorted(index_map.items())) #(Chain_ID, resnum) -> index
        self.mat = cdist(self.coords, self.coords, metric='euclidean')

    def get_distance(self, start_key, end_key):
        i = self.index_map.get(start_key)
        j = self.index_map.get(end_key)
        if i is None or j is None:
            raise KeyError("Residue key not found.")
        return self.mat[i, j]
    
    def convert_to_df(self):
        keys = list(self.index_map.keys())
        n = len(keys)

        #prevent duplicates in datatable
        triu_i, triu_j = np.triu_indices(n, k=1)

        df = pd.DataFrame({
            'start_key': [f"{keys[i][0]}-{keys[i][1]}" for i in triu_i],
            'end_key': [f"{keys[j][0]}-{keys[j][1]}" for j in triu_j],
            'start_coord': [self.coords[i].tolist() for i in triu_i],  
            'end_coord': [self.coords[j].tolist() for j in triu_j],
            'distance': self.mat[triu_i, triu_j]
        })

        return df
    

class CompareDistanceMatrix:
    def __init__(self, reference_matrix,  target_matrix, res_id_mapping):
        self.shared_keys = set(reference_matrix.index_map.keys()) & set(target_matrix.index_map.keys())
        self.res_id_mapping = res_id_mapping
        self.index_map = {key: i for i, key in enumerate(sorted(self.shared_keys))}
        self.ref_coords = [reference_matrix.coords[reference_matrix.index_map[key]] for key in sorted(self.shared_keys)]
        self.tgt_coords = [target_matrix.coords[target_matrix.index_map[key]] for key in sorted(self.shared_keys)]

        self.ref_mat = cdist(self.ref_coords, self.ref_coords)
        self.tgt_mat = cdist(self.tgt_coords, self.tgt_coords)
        self.mat = self.tgt_mat - self.ref_mat #comparison matrix

    def get_distance_diff(self, start_key, end_key):
        i = self.index_map[start_key]
        j = self.index_map[end_key]
        return self.mat[i, j]
    
    def get_coords(self, key):
        return self.ref_coords[self.index_map[key]], self.tgt_coords[self.index_map[key]]

    def convert_to_df(self):

        keys = list(self.shared_keys)
        n = len(keys)

        triu_i, triu_j = np.triu_indices(n, k=1)

        df = pd.DataFrame({
            'start_key_ref': [f"{keys[i][0]}-{keys[i][1]}" for i in triu_i],
            'end_key_ref': [f"{keys[j][0]}-{keys[j][1]}" for j in triu_j],
            'start_coord_ref': [self.ref_coords[i].tolist() for i in triu_i],  # convert np arrays to lists
            'end_coord_ref': [self.ref_coords[j].tolist() for j in triu_j],
            'start_key_tgt': [f"{self.res_id_mapping[keys[i]][0]}-{self.res_id_mapping[keys[i]][1]}" for i in triu_i],
            'end_key_tgt': [f"{self.res_id_mapping[keys[j]][0]}-{self.res_id_mapping[keys[j]][1]}" for j in triu_j],
            'start_coord_tgt': [self.tgt_coords[i].tolist() for i in triu_i],  # convert np arrays to lists
            'end_coord_tgt': [self.tgt_coords[j].tolist() for j in triu_j],
            '∆ distance': self.mat[triu_i, triu_j]
        })

        return df




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

def calc_difference_from_mapper(structure_mapper,explicit_chain_mapping):
    if explicit_chain_mapping:
        structure_mapper.map_chain_explicit(explicit_chain_mapping)

    else:
        structure_mapper.map_chains(threshold=95)
    data = []

    for chain_mapping in structure_mapper.chain_mappings.values():
        ref_id = chain_mapping.ref_chain.id
        tgt_id = chain_mapping.tgt_chain.id

        chain_mapping.calc_aligned_coords()
        for ref_res_id, tgt_res_id in chain_mapping.res_id_map.items():


            ref_coord = chain_mapping.get_ref_coord(ref_res_id[0], ref_res_id[1])
            tgt_coord = chain_mapping.get_tgt_coord(ref_res_id[0], ref_res_id[1])

            diff_vec = tgt_coord - ref_coord

            distance = np.linalg.norm(diff_vec)

            data.append({
            'ChainID_Resnum1': f'{ref_id}_{ref_res_id[1]}',
            'ChainID_Resnum2': f'{tgt_id}_{tgt_res_id[1]}',
            'Coord1': ref_coord.tolist(),
            'Coord2': tgt_coord.tolist(),
            'Diff_Vec': diff_vec.tolist(),
            'Distance': distance
            })
    return pd.DataFrame(data).dropna()