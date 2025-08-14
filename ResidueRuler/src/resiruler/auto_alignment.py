import Bio
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.PDB import MMCIFParser, MMCIFIO, Structure, Model
from io import StringIO
from structure_parsing import extract_res_from_chain, extract_seq_from_chain
from distance_calc import DistanceMatrix, CompareDistanceMatrix
import numpy as np
from scipy.optimize import linear_sum_assignment
import copy
import pandas as pd 

class ChainMapper:
    """
    ChainMapper class maps two corresponding chains in a reference and target structure
    """
    def __init__(self, ref_chain, ref_seq, tgt_chain, tgt_seq, alignment):
        self.ref_chain = ref_chain
        self.ref_seq = ref_seq
        
        self.tgt_chain = tgt_chain
        self.tgt_seq = tgt_seq

        self.alignment = alignment

        self.aligned_ref_seq = alignment[0]
        self.aligned_tgt_seq = alignment[1]

        self.res_id_mapping = self.calc_residue_mapping()

    def calc_percent_identity(self):
        """
        Calculates percent identity (ignores gaps)
        """
        aligned_ref_str = self.aligned_ref_seq
        aligned_tgt_str = self.aligned_tgt_seq
        assert len(aligned_ref_str) == len(aligned_tgt_str)

        # Count matches (ignore gaps)
        matches = sum(
            a == b and a != '-' and b != '-'
            for a, b in zip(aligned_ref_str, aligned_tgt_str)
        )

        # Identity = matches / length of aligned region (excluding gaps)
        aligned_length = sum(
            a != '-' and b != '-' for a, b in zip(aligned_ref_str, aligned_tgt_str)
        )

        percent_identity = 100 * matches / aligned_length if aligned_length > 0 else 0.0

        return percent_identity

    def calc_residue_mapping(self):
        """
        Gets aligned residues and stores them in a dictionary that is ref(ChainID, ResID) -> tgt(ChainID, ResID)
        ResID are biopython residue ID's which take the form of ('heteroflag','residue number','insertion code')
        """
        ref_res = extract_res_from_chain(self.ref_chain)
        tgt_res = extract_res_from_chain(self.tgt_chain)

        aligned_ref_seq = str(self.alignment[0])
        aligned_tgt_seq = str(self.alignment[1])

        assert len(aligned_ref_seq) == len(aligned_tgt_seq)

        res_id_mapping = {}
        idx_ref = 0
        idx_tgt = 0

        # for every letter (residue) in the aligned sequences
        for i, (ref_char, tgt_char) in enumerate(zip(aligned_ref_seq, aligned_tgt_seq)):
            #No gap in ref, so there is a residue there and we "select" it 
            if ref_char != "-":
                a_res = ref_res[idx_ref]
                idx_ref += 1
            else:
                a_res = None

            if tgt_char != "-":
                #No gap in tgt, so there is a residue there and we "select" it 
                b_res = tgt_res[idx_tgt]
                idx_tgt += 1
            else:
                b_res = None

            #we have a matching set, so we can add it to the aligned coords and update our mapping 
            if a_res is not None and b_res is not None:
                    ref_id = (self.ref_chain.id, a_res.id)
                    tgt_id = (self.tgt_chain.id, b_res.id)
                    res_id_mapping[ref_id] = tgt_id

        return res_id_mapping
    
    def get_aligned_coord_lists(self):
        """
        Calculates the aligned CA coordinate lists for residues between mapped chains
        These can be index via the index map which stores a simple way to use 
        """
        index_map={}
        current_index = 0
        aligned_ref_coords=[]
        aligned_tgt_coords=[]

        for (ref_chain_id, ref_res_id), (tgt_chain_id, tgt_res_id) in self.res_id_mapping.items():
            ref_res = self.ref_chain[ref_res_id]
            tgt_res = self.tgt_chain[tgt_res_id]
            if "CA" in ref_res and "CA" in tgt_res:
                aligned_ref_coords.append(ref_res["CA"].get_coord())
                aligned_tgt_coords.append(tgt_res["CA"].get_coord())
                index_map[(ref_chain_id,ref_res_id)] = current_index
                current_index += 1

        return aligned_ref_coords, aligned_tgt_coords, index_map
        

    def get_ref_coord(self, res_id):
        """
        Gets CA coordinate of residue in reference structure based on Resid) identifier
        """
        return self.ref_chain[res_id]["CA"].get_coord()

    def get_tgt_coord(self, res_id):
        """
        Gets coordinate of residue in reference structure based on reference Resid identifier
        """
        return self.tgt_chain[self.res_id_mapping[res_id]]["CA"].get_coord()

    def calc_rmsd(self, coords1, coords2):
        """
        Calculate rmsd between two sets of aligned coords
        This assumes that the coords are aligned, so make sure you don't put unaligned lists of coordinates here
        """
        coords1 = np.array(coords1)
        coords2 = np.array(coords2)
        assert coords1.shape == coords2.shape
        if coords1.ndim == 1:
            coords1 = coords1.reshape(1, -1)
            coords2 = coords2.reshape(1, -1)

        diff = coords1 - coords2
        return np.sqrt(np.mean(np.sum(diff**2, axis=1)))
        

class StructureMapper:
    """
    Structure Mapper class to map a target structure to a reference
    Mapping is based on PCT identities and tiebreaks based on RMSD. It is most recommended to map only using aligned structures
    Otherwise the map may struggle a bit (This problem only really persists for symmetric assemblies however)
    """
    def __init__(self, ref_structure, tgt_structure, aligner=None):
        #self.map_id = map_id
        self.ref_structure = ref_structure
        self.tgt_structure = tgt_structure
        self.aligner = self.aligner = aligner or self._default_aligner()
        self.chain_mappings = {} # ref_chain_id -> ChainMapping object 

        self.matched_ref_chains = set()
        self.matched_tgt_chains = set()
        self._cached_sequences = {}

    def _default_aligner(self):
        '''
        Pairwise alignment that uses the Blosum62 substition matrix and Needleman-wunsch algorithm
        '''
        aligner = PairwiseAligner()
        aligner.mode = "global"
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        aligner.left_open_gap_score = 1
        
        return aligner
    
    def get_chain_sequence(self, structure, chain):
        key = (id(structure), chain.id)
        if key not in self._cached_sequences:
            seq = extract_seq_from_chain(chain)
            self._cached_sequences[key] = seq
        return self._cached_sequences[key]

    def map_chains(self,threshold):
        """
        Map all chains into ChainMapper Objects based on first sequence, then break ties with RMSD
        Assignmnet of matches occurs via solving the linear sum assignment problem  
        """
        ref_structure = self.ref_structure
        tgt_structure = self.tgt_structure

        ref_chains = list(ref_structure.get_chains())
        tgt_chains = list(tgt_structure.get_chains())

        valid_pairs = []
        scores = {} #(i,j) -> score
        mappings = {} #(i,j) - > ChainMapper

        # Collect valid combinations of chains (above the pct identity threshold)
        # So that we can assign matches via the optimal linear sum assignmnet
        for i, ref_chain in enumerate(ref_chains):
            ref_seq = self.get_chain_sequence(ref_structure, ref_chain)

            for j, tgt_chain in enumerate(tgt_chains):
                tgt_seq = self.get_chain_sequence(tgt_structure, tgt_chain)

                alignment = self.aligner.align(ref_seq, tgt_seq)[0]
                potential_map = ChainMapper(ref_chain, ref_seq, tgt_chain, tgt_seq, alignment)

                percent_identity = potential_map.calc_percent_identity()
                #print(f"Testing chains {ref_chain.id} vs {tgt_chain.id} — %ID: {percent_identity:.2f}")
                

                if percent_identity < threshold:
                    continue

                aligned_ref_coords, aligned_tgt_coords, _ = potential_map.get_aligned_coord_lists()
                rmsd = potential_map.calc_rmsd(
                    aligned_ref_coords,
                    aligned_tgt_coords
                )

                score = percent_identity - 1e-6 * rmsd  # prioritize identity, break ties with RMSD
                valid_pairs.append((i, j))
                scores[(i, j)] = score
                mappings[(i, j)] = potential_map

        if not valid_pairs:
            print(" No valid chain matches above threshold.")
            return

        # Build score matrix using only the valid matchings
        rows = sorted(set(i for i, _ in valid_pairs))
        cols = sorted(set(j for _, j in valid_pairs))
        row_idx_map = {i: idx for idx, i in enumerate(rows)}
        col_idx_map = {j: idx for idx, j in enumerate(cols)}

        reduced_matrix = np.full((len(rows), len(cols)), 1e9)  # high cost by default, this will be the value for an unmatched reference chain


        for (i, j), score in scores.items():
            r = row_idx_map[i]
            c = col_idx_map[j]
            reduced_matrix[r, c] = -score  # negate as linear sum assignment finds the optimum minimum solution

       
        row_ind, col_ind = linear_sum_assignment(reduced_matrix)

        # Recover global i,j which is usable to get the mappings
        for r_idx, c_idx in zip(row_ind, col_ind):
            #Get the indices for the actual mappings
            i = rows[r_idx]
            j = cols[c_idx]

            if (i, j) not in mappings:
                print(f"Skipping assignment ({i}, {j}) — not in valid mappings")
                continue

            mapping = mappings[(i, j)]
            ref_id = mapping.ref_chain.id
            tgt_id = mapping.tgt_chain.id

            self.chain_mappings[ref_id] = mapping
            self.matched_ref_chains.add(ref_id)
            self.matched_tgt_chains.add(tgt_id)
            print(f"Matched {ref_id} → {tgt_id} with score {scores[(i, j)]:.4f}")

    def map_chains_explicit(self, explicit_chain_mapping):
        for ref_chain_id, tgt_chain_id in explicit_chain_mapping.items():
            ref_chain = self.ref_structure[ref_chain_id]
            tgt_chain = self.tgt_structure[tgt_chain_id]

            self.matched_ref_chains.add(ref_chain_id)
            self.matched_tgt_chains.add(tgt_chain_id)

            ref_seq = extract_seq_from_chain(ref_chain)
            tgt_seq = extract_seq_from_chain(tgt_chain)

            alignment = self.aligner.align(ref_seq, tgt_seq)[0]

            self.chain_mappings[ref_chain_id] = ChainMapper(ref_chain, ref_seq, tgt_chain, tgt_seq, alignment)

    def calc_matrices (self, selected_chains=None):
        """
        Create DistanceMatrix and CompareDistanceMatrix for the reference and targets.
        """
        
        coords_ref, coords_tgt, index_map, res_id_map = self.get_selected_mapping(selected_chains)

        if len(coords_ref) == 0:
            raise ValueError("No aligned coordinates were found. Check selected chains and mapping.")

        coords_ref = np.array(coords_ref)
        coords_tgt = np.array(coords_tgt)

        if coords_ref.ndim != 2 or coords_ref.shape[1] != 3:
            raise ValueError(f"ref coords shape: {coords_ref.shape} — expected (N, 3)")

        ref_dm = DistanceMatrix(coords_ref, index_map)
        tgt_dm = DistanceMatrix(coords_tgt, index_map)
        return ref_dm, tgt_dm, CompareDistanceMatrix(ref_dm, tgt_dm, res_id_map)
    
    def get_selected_mapping(self, selected_chains=None):
        """
        Extract a mapping from the ChainMapper Objects based on a selection of chains
        """
        coords_ref_list = []
        coords_tgt_list = []
        index_map = {} 
        res_id_map = {}

        coords_index = 0
        for chain_id, cm in self.chain_mappings.items():
            if selected_chains and chain_id not in selected_chains:
                print(f"[INFO] Skipping chain {chain_id} (not in selected_chains)")
                continue
            
            aligned_ref_coords, aligned_tgt_coords, chain_index_map = cm.get_aligned_coord_lists()
            
            coords_ref_list.append(np.array(aligned_ref_coords))
            coords_tgt_list.append(np.array(aligned_tgt_coords))

            for ref_key, chain_index in chain_index_map.items():
                index_map[ref_key] = coords_index + chain_index

            coords_index += len(chain_index_map)
            res_id_map.update(cm.res_id_mapping)

        if len(coords_ref_list) == 0:
            raise ValueError("No aligned coordinates were found. Check selected chains and mapping.")

        coords_ref = np.vstack(coords_ref_list)
        coords_tgt = np.vstack(coords_tgt_list)
        
        return coords_ref, coords_tgt, index_map, res_id_map


class EnsembleMapper:
    """
    This is an Ensemble Mapper Class designed to handle mutiple structure mapping objects
    As currently implemented, it will store relevant mappings for selected chains calculated via the get_selected_global_coords class method
    All coordinates are aligned to a single set of reference cordinates and thus corresponding coordinates can all be found using the same index
    """

    def __init__(self, ref_structure, aligner):
        self.ref_structure=ref_structure
        self.aligner = aligner
        self.coords_ref = None
        self.structure_mappings = {} # Structure Name -> StructureMapper
        self.global_index_mapping = {} # (ChainID, ResID) -> index (for use with aligned tgt coords)
        self.res_id_mappings = {} # Structure Name -> (ref(ChainID, ResID) -> tgt(ChainID, ResID))
        self.coords_targets_dict= {} # Structure Name -> aligned_tgt_coords
    
    def add_structure(self, tgt_structure_name, tgt_structure, threshold, explicit_mapping=None):
        structure_mapping = StructureMapper(self.ref_structure, tgt_structure, self.aligner)
        if explicit_mapping:
            structure_mapping.map_chains_explicit(explicit_mapping)
        else:
            structure_mapping.map_chains(threshold)
        self.structure_mappings[tgt_structure_name] = structure_mapping
    
    def get_common_ref_residues(self, selected_chains):
        """
        Gets only the residues that are matched in all mapped structures
        """
        if selected_chains is None:
            selected_chains = [chain.id for chain in self.ref_structure.get_chains()]

        if selected_chains is None:
            #Use all
            selected_chains = {chain_id 
                  for struct_map in self.structure_mappings.values() 
                  for chain_id in struct_map.chain_mappings.keys()
                  }
        
        
        #get a set of all matched residues
        mapped_residues_sets = []
        for structure_mapping in self.structure_mappings.values():
            mapped_ref_residues = set()
            for chain_id in selected_chains:
                if chain_id not in structure_mapping.matched_ref_chains:
                    print(f"{chain_id} not mapped")
                    continue
                print(structure_mapping.chain_mappings[chain_id].ref_seq)
                print(structure_mapping.chain_mappings[chain_id].tgt_seq)
                #Add keys of res_id map from chain mapping (matched residues)
                mapped_ref_residues.update(structure_mapping.chain_mappings[chain_id].res_id_mapping.keys())
            
            mapped_residues_sets.append(mapped_ref_residues)
        
        if not mapped_residues_sets:
            return set()
        
        common_ref_residues = set.intersection(*mapped_residues_sets)

        
        return common_ref_residues


    def set_selected_global_coords(self, selected_chains=None):
        if not self.structure_mappings:
            raise ValueError("No target structures added yet.")
            
        if selected_chains is None:
            selected_chains = [chain.id for chain in self.ref_structure.get_chains()]
        
        #now we have a set of residues that got matched in all structure mappings
        common_ref_residues = self.get_common_ref_residues(selected_chains)

        #chains_included = {c for c, _ in common_ref_residues}
       

        sorted_residues = sorted(common_ref_residues) # consistent order (Alphabetical by chain first, then residue number, then insertion code)


        first_mapping = next(iter(self.structure_mappings.values()))
        coords_ref_list = []
        global_index_mapping = {}
        
        coords_ref,_ ,index_map , _ = first_mapping.get_selected_mapping(selected_chains)
        
        #Get the reference coordinates
        for idx, res_id in enumerate(sorted_residues):

            coords_ref_list.append(coords_ref[index_map[res_id]])
            global_index_mapping[res_id] = idx

        self.coords_ref = np.array(coords_ref_list)
        self.global_index_mapping = global_index_mapping


        self.coords_targets_dict.clear()
        self.res_id_mappings.clear()
        #Map the target coordinates for each mapper
        for name, structure_mapping in self.structure_mappings.items():

            _,coords_tgt,index_map,res_id_map = structure_mapping.get_selected_mapping(selected_chains)

            #filter res_id map to only include the relevant residues, now we can assume anything in this is actually mapped in across all structures
            #Only necessary if we iterate through this mapping, otherwise we can get away with not filtering 
            filtered_res_id_map = {chain_res_ref: chain_res_tgt for chain_res_ref, chain_res_tgt in res_id_map.items() if chain_res_ref in common_ref_residues}
            
            remapped_coords_tgt = []

            for ref_id in sorted_residues:
                remapped_coords_tgt.append(coords_tgt[index_map[ref_id]])

            self.coords_targets_dict[name] = np.array(remapped_coords_tgt)
            self.res_id_mappings[name] = filtered_res_id_map
            if filtered_res_id_map is None:
                print(f"RESID MAP IS NONE FOR {name}")
            

    def calc_matrices(self):
        """
        Creates a dictionary of DistanceMatrix and CompareMatrix Objects and returns it 
        This should only be used after set_selected_global_coords has been called/used
        """
        tgt_dms = {}
        compare_dms = {}
        index_map = self.global_index_mapping
        aligned_ref_coords = self.coords_ref
        ref_dm = DistanceMatrix(aligned_ref_coords, index_map)
        for structure_name in self.structure_mappings.keys():
            res_id_map = self.res_id_mappings.get(structure_name, None)
            if not res_id_map:
                print(f"[WARNING] Skipping {structure_name}: no residue ID mapping found.")
                continue  # skip this one
            
            aligned_tgt_coords = self.coords_targets_dict[structure_name]
            tgt_dm = DistanceMatrix(aligned_tgt_coords, index_map, res_id_map)
            compare_dm = CompareDistanceMatrix(ref_dm, tgt_dm, res_id_map)

            tgt_dms[structure_name] = tgt_dm
            compare_dms[structure_name] = compare_dm
                
        return ref_dm, tgt_dms, compare_dms
    
    def calc_movement_dfs_rmsds(self):
        """
        Calculate movement dfs and global rmsd for each structure mapped to the reference, by taking the common residues and finding their difference with corresponding 
        residues in the reference 
        """
        movement_dfs = {}
        rmsds = {}
        ref_coords = self.coords_ref
        global_index_map = self.global_index_mapping
        for structure_name, tgt_coords in self.coords_targets_dict.items():
            
            tgt_coords = self.coords_targets_dict[structure_name]

            #distances/difference vectors from residue in reference structure to corresponding residue in the next
            diff_vecs = ref_coords - tgt_coords
            distances = np.linalg.norm(diff_vecs, axis=1)
            rmsd = np.sqrt(np.mean(np.sum(diff_vecs**2, axis=1)))

            ref_ids, tgt_ids = zip(*self.res_id_mappings[structure_name].items())
            print("Ref ID Length", len(ref_ids))
            print("Tgt ID Length", len(tgt_ids))
            print("Ref Coords Length", len(ref_coords))
            print("TGT Coords Length", len(tgt_coords))
            df = pd.DataFrame({
            "ChainID_Resnum1": ref_ids,
            "ChainID_Resnum2": tgt_ids,
            "Coord1": np.array(ref_coords).tolist(),
            "Coord2": np.array(tgt_coords).tolist(),
            "Diff_Vec": np.array(diff_vecs).tolist(),
            "Distance": distances
            })

            movement_dfs[structure_name] = df
            rmsds[structure_name] = rmsd

        return movement_dfs, rmsds


def write_filtered_structure(structure, matched_chains=None, matched_residues=None):
    """
    Writes a filtered structure which includes either only the chains with an associted match, or only the residues with an associated match between tgt and reference
    """

    #Store new writes
    new_structure = Structure.Structure(f"{structure.id}_filtered")
    new_model = Model.Model(0)
    new_structure.add(new_model)

    for chain in structure.get_chains():
        #No match, or not looking for filtering chain
        if matched_chains is not None and chain.id not in matched_chains:
            continue

        new_chain = chain.__class__(chain.id)
        for res in chain.get_residues():
            res_copy = copy.deepcopy()
            key = (chain.id, res.id)
            #No residue match or not filtering residues
            if matched_residues is not None and key not in matched_residues:
                continue
            new_chain.add(res_copy)

        #If we've got a chain add it the new structure
        if len(new_chain):
            new_model.add(new_chain)

    io_buffer = StringIO()
    io = MMCIFIO()
    io.set_structure(new_structure)
    io.save(io_buffer)
    return io_buffer.getvalue()  # return CIF string

def filter_and_write_aligned_maps(ref_cif, tgt_cif, identity_threshold=95.0):
    """
    Filter aligned models, and write out new models which include only matched residues as well as models which include only matched chains 
    """
    parser = MMCIFParser(QUIET=True)
    ref_structure = parser.get_structure("ref", ref_cif)
    tgt_structure = parser.get_structure("tgt", tgt_cif)

    mapper = StructureMapper(ref_structure, tgt_structure)
    mapper.map_chains(threshold=identity_threshold)

    matched_ref_residues = set()
    matched_tgt_residues = set()
    matched_ref_chains = set()
    matched_tgt_chains = set()

    #Get matches residues and chains 
    for chain_id, cm in mapper.chain_mappings.items():
        for ref_res_id, tgt_res_id in cm.res_id_mapping.items():
            matched_ref_residues.add((ref_res_id)) # (chainID, resID)
            
            matched_tgt_residues.add((tgt_res_id)) # (chainID, resID)

        matched_ref_chains.add(cm.ref_chain.id)
        matched_tgt_chains.add(cm.tgt_chain.id)

    ref_chain_str = write_filtered_structure(ref_structure, matched_chains=matched_ref_chains)
    tgt_chain_str = write_filtered_structure(tgt_structure, matched_chains=matched_tgt_chains)

    ref_residue_str = write_filtered_structure(ref_structure,
                                               matched_chains=matched_ref_chains,
                                               matched_residues=matched_ref_residues)

    tgt_residue_str = write_filtered_structure(tgt_structure,
                                               matched_chains=matched_tgt_chains,
                                               matched_residues=matched_tgt_residues)

    return ref_chain_str, ref_residue_str, tgt_chain_str, tgt_residue_str