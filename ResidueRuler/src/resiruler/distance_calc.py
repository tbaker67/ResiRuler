import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist

class DistanceMatrix:
    def __init__(self, coords_list, index_map, res_id_map=None):
        self.res_id_mapping = res_id_map
        self.coords = np.array(coords_list, dtype=np.float32)
        self.index_map = dict(sorted(index_map.items())) #(Chain_ID, resnum) -> index
        self.mat = cdist(self.coords, self.coords, metric='euclidean')

    def get_distance(self, start_key, end_key):
        i = self.index_map.get(start_key)
        j = self.index_map.get(end_key)
        if i is None or j is None:
            raise KeyError("Residue key not found.")
        return self.mat[i, j]
    
    def convert_to_df(self, lower_threshold=None, upper_threshold=None):
        keys = list(self.index_map.keys())
        n = len(keys)

        #prevent duplicates in datatable
        triu_i, triu_j = np.triu_indices(n, k=1)

        distances = self.mat[triu_i, triu_j]

        # allow for thresholding to save memory
        if lower_threshold is not None and upper_threshold is not None:
            mask = (distances > lower_threshold) & (distances < upper_threshold)
            triu_i = triu_i[mask]
            triu_j = triu_j[mask]
            distances = distances[mask]

        df = pd.DataFrame({
            'ChainID_Resnum1': [f"{keys[i][0]}-{keys[i][1][1]}{keys[i][1][2]}" for i in triu_i],
            'ChainID_Resnum2': [f"{keys[j][0]}-{keys[j][1][1]}{keys[j][1][2]}" for j in triu_j],
            #'Coord1': [self.coords[i].tolist() for i in triu_i],  
            #'Coord2': [self.coords[j].tolist() for j in triu_j],
            'Distance': distances
        })

        return df
    def convert_to_df_chunked(self, lower_threshold=None, upper_threshold=None, chunk_size=5000):
        keys = list(self.index_map.keys())
        n = len(keys)

        dfs = []

        for start_i in range(0, n, chunk_size):
            end_i = min(start_i + chunk_size, n)
            
            # Only compute upper triangle for this chunk
            for i in range(start_i, end_i):
                j_start = i + 1
                if j_start >= n:
                    continue
                
                j_indices = np.arange(j_start, n)
                distances = self.mat[i, j_indices]

                # Threshold filter to save memory
                if lower_threshold is not None and upper_threshold is not None:
                    mask = (distances > lower_threshold) & (distances < upper_threshold)
                    if not np.any(mask):
                        continue
                    j_indices = j_indices[mask]
                    distances = distances[mask]

                # Create small DataFrame for this row
                df_chunk = pd.DataFrame({
                    'Res1': [f"{keys[i][0]}-{keys[i][1][1]}{keys[i][1][2]}"] * len(j_indices),
                    'Res2': [f"{keys[j][0]}-{keys[j][1][1]}{keys[j][1][2]}" for j in j_indices],
                    'Distance': distances
                })

                dfs.append(df_chunk)

        # Concatenate only at the end
        return pd.concat(dfs, ignore_index=True)
    

class CompareDistanceMatrix:
    def __init__(self, reference_matrix,  target_matrix, res_id_mapping):
        self.shared_keys = set(reference_matrix.index_map.keys()) & set(target_matrix.index_map.keys())
        self.res_id_mapping = res_id_mapping
        self.index_map = {key: i for i, key in enumerate(sorted(self.shared_keys))}
        self.ref_coords = np.array([reference_matrix.coords[reference_matrix.index_map[key]] for key in sorted(self.shared_keys)], dtype=np.float32)
        self.tgt_coords = np.array([target_matrix.coords[target_matrix.index_map[key]] for key in sorted(self.shared_keys)], dtype=np.float32)

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
            'ChainID_Resnum1_ref': [f"{keys[i][0]}-{keys[i][1][1]}{keys[i][1][2]}" for i in triu_i],
            'ChainID_Resnum2_ref': [f"{keys[j][0]}-{keys[j][1][1]}{keys[j][1][2]}" for j in triu_j],
            'Coord1_ref': [self.ref_coords[i].tolist() for i in triu_i],  # convert np arrays to lists
            'Coord2_ref': [self.ref_coords[j].tolist() for j in triu_j],
            'ChainID_Resnum1_tgt': [f"{self.res_id_mapping[keys[i]][0]}-{self.res_id_mapping[keys[i]][1][1]}{self.res_id_mapping[keys[i]][1][2]}" for i in triu_i],
            'ChainID_Resnum2_tgt': [f"{self.res_id_mapping[keys[j]][0]}-{self.res_id_mapping[keys[j]][1][1]}{self.res_id_mapping[keys[j]][1][2]}" for j in triu_j],
            'Coord1_tgt': [self.tgt_coords[i].tolist() for i in triu_i],  # convert np arrays to lists
            'Coord2_tgt': [self.tgt_coords[j].tolist() for j in triu_j],
            '∆ distance': self.mat[triu_i, triu_j]
        })

        return df
    
def compute_distance(coord1, coord2):
    if coord1 is None or coord2 is None:
        return np.nan
    #print(f"[DEBUG] Computing distance between {coord1} and {coord2}")
    return np.linalg.norm(np.array(coord1) - np.array(coord2))


def calc_difference_from_mapper(structure_mapper,explicit_chain_mapping=None):
    if explicit_chain_mapping:
        structure_mapper.map_chain_explicit(explicit_chain_mapping)

    else:
        structure_mapper.map_chains(threshold=95)
    
    ref_ids = []
    tgt_ids = []
    ref_coords = []
    tgt_coords = []
    diffs = []
    distances = []

    for chain_mapping in structure_mapper.chain_mappings.values():
        ref_id = chain_mapping.ref_chain.id
        tgt_id = chain_mapping.tgt_chain.id

        chain_mapping.calc_aligned_coords()

        for ref_res_id, tgt_res_id in chain_mapping.res_id_map.items():
            try:
                ref_coord = chain_mapping.get_ref_coord(ref_res_id[0], ref_res_id[1])
                tgt_coord = chain_mapping.get_tgt_coord(ref_res_id[0], ref_res_id[1])
            except (KeyError, IndexError):
                continue

            diff_vec = tgt_coord - ref_coord
            dist = np.linalg.norm(diff_vec)

            ref_ids.append(f"{ref_id}_{ref_res_id[1]}{ref_res_id[2]}")
            tgt_ids.append(f"{tgt_id}_{tgt_res_id[1]}{tgt_res_id[2]}")
            ref_coords.append(ref_coord)
            tgt_coords.append(tgt_coord)
            diffs.append(diff_vec)
            distances.append(dist)

    df = pd.DataFrame({
        "ChainID_Resnum1": ref_ids,
        "ChainID_Resnum2": tgt_ids,
        "Coord1": np.array(ref_coords).tolist(),
        "Coord2": np.array(tgt_coords).tolist(),
        "Diff_Vec": np.array(diffs).tolist(),
        "Distance": distances
    })

    return df.dropna()