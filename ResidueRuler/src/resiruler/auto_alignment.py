import Bio
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.PDB import MMCIFParser, MMCIFIO, Structure, Model
from io import StringIO
from .structure_parsing import extract_res_from_chain, extract_seq_from_chain, get_CA_from_residue, get_CB_from_residue, get_SC_from_residue, get_C1prime_from_residue, ChainCollection
from .distance_calc import DistanceMatrix, CompareDistanceMatrix
import numpy as np
from scipy.optimize import linear_sum_assignment
import copy
import pandas as pd 

class ChainMapper:
    """
    ChainMapper class maps two corresponding chains in a reference and target structure
    """
    def __init__(self, ref_chain, ref_seq, aligned_ref_seq, tgt_chain, tgt_seq, aligned_tgt_seq, alignment, type):
        self.ref_chain = ref_chain
        self.ref_seq = ref_seq
        self.ref_start=None
        
        self.tgt_chain = tgt_chain
        self.tgt_seq = tgt_seq
        self.tgt_start=None

        self.alignment = alignment

        self.aligned_ref_seq = aligned_ref_seq
        self.aligned_tgt_seq = aligned_tgt_seq

        self.type = type
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
        self.ref_start = None
        self.tgt_start = None
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
    
    def get_aligned_coord_lists(self, mode):
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

    

            ref_coord = None
            tgt_coord = None

            if mode == "CA":
                ref_coord = get_CA_from_residue(ref_res)
                tgt_coord = get_CA_from_residue(tgt_res)
            elif mode == "CB":
                ref_coord = get_CB_from_residue(ref_res)
                tgt_coord = get_CB_from_residue(tgt_res)
            
            #measure sidechain distances (heavy carbons not including CA)
            elif mode == "SC":
                ref_coord = get_SC_from_residue(ref_res)
                tgt_coord = get_SC_from_residue(tgt_res)
            
            elif mode =="C1'":
                ref_coord = get_C1prime_from_residue(ref_res)
                tgt_coord = get_C1prime_from_residue(tgt_res)
            
            else:
                print("Improper Mode Selected")

            if ref_coord is None or tgt_coord is None:
                print(f"[WARNING] skipping {ref_res_id} <-> {tgt_res_id} due to missing coordinate")
                continue 

            aligned_ref_coords.append(ref_coord)
            aligned_tgt_coords.append(tgt_coord)
            
            index_map[(ref_chain_id, ref_res_id)] = current_index
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
    def __init__(self, ref_structure, tgt_structure, protein_aligner, nucleotide_aligner):
        #self.map_id = map_id
        self.ref_structure = ref_structure
        self.tgt_structure = tgt_structure
        self.chain_mappings = {} # ref_chain_id -> ChainMapping object 

        self.protein_aligner = protein_aligner
        self.nucleotide_aligner = nucleotide_aligner

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

        ref_collection = ChainCollection(ref_structure)
        tgt_collection = ChainCollection(tgt_structure)

        valid_pairs = []
        scores = {} #(i,j) -> score
        mappings = {} #(i,j) - > ChainMapper

        # Collect valid combinations of chains (above the pct identity threshold)
        # So that we can assign matches via the optimal linear sum assignmnet

        dna_aligner = PairwiseAligner()
        dna_aligner.mode='global'
        
        dna_aligner.open_gap_score = -10
        dna_aligner.extend_gap_score = -1
        dna_aligner.match_score=5
        dna_aligner.mismatch_score=-4

        for ref_id, tgt_id, ref_info, tgt_info in ref_collection.valid_pairs(tgt_collection):

        
            
            if ref_info.type == "protein":
                alignment = self.protein_aligner.align(ref_info.seq, tgt_info.seq)[0]
                aligned_ref_seq = alignment[0]
                aligned_tgt_seq = alignment[1]
            elif ref_info.type == "dna":
                alignment = self.nucleotide_aligner.align(ref_info.seq, tgt_info.seq)[0]
                aligned_ref_seq = alignment[0]
                aligned_tgt_seq = alignment[1]

            elif ref_info.type == 'rna':
                #Substitution Matrices for the pairwise aligners are all set up to use T rather than U
                #So we must replace U with T before aligning and then shift it back
                converted_ref_seq = ref_info.seq.replace("U","T")
                converted_tgt_seq = tgt_info.seq.replace("U","T")

                alignment = self.nucleotide_aligner.align(converted_ref_seq, converted_tgt_seq)[0]
                #Change T's back to U's
                aligned_ref_seq = alignment[0].replace("T","U")
                aligned_tgt_seq = alignment[1].replace("T","U")
            
            #unknown chain type so just skip it
            else:
                print(f"[WARNING] skipping chain {ref_id} due to an unknown chain type")
                continue

            potential_map = ChainMapper(ref_info.chain, ref_info.seq,aligned_ref_seq, tgt_info.chain, tgt_info.seq,aligned_tgt_seq, alignment, ref_info.type)

            percent_identity = potential_map.calc_percent_identity()
            
                

            if percent_identity < threshold:
                continue
            
            # rmsd calculated using CA for protein and C1' for nucleic acids 
            if potential_map.type == "protein":
                aligned_ref_coords, aligned_tgt_coords, _ = potential_map.get_aligned_coord_lists("CA")
            
            elif potential_map.type == "dna" or potential_map.type == "rna":
                aligned_ref_coords, aligned_tgt_coords, _ = potential_map.get_aligned_coord_lists("C1'")

            rmsd = potential_map.calc_rmsd(
                aligned_ref_coords,
                aligned_tgt_coords
            )


            #no division by zero unless pct id threshold <= 0

            score = (alignment.score / len(aligned_ref_seq.replace("-", ""))) - 0.1 * rmsd # prioritize alignment score, break ties with RMSD
            valid_pairs.append((ref_id, tgt_id))
            scores[(ref_id, tgt_id)] = score
            mappings[(ref_id, tgt_id)] = potential_map

        if not valid_pairs:
            print(" No valid chain matches above threshold.")
            return
        
        # Build score matrix using only the valid matchings
        rows = sorted(set(ref_id for ref_id, _ in valid_pairs))
        cols = sorted(set(tgt_id for _, tgt_id in valid_pairs))

        row_idx_map = {ref_id: i for i, ref_id in enumerate(rows)}

        col_idx_map = {tgt_id: j for j, tgt_id in enumerate(cols)}

        cost_matrix = np.full((len(rows), len(cols)), 1e9)  # high cost by default, this will be the value for an unmatched reference chain


        for (ref_id, tgt_id), score in scores.items():
            i = row_idx_map[ref_id]
            j = col_idx_map[tgt_id]
            cost_matrix[i, j] = -score  # negate as linear sum assignment finds the optimum minimum solution

       
        row_ind, col_ind = linear_sum_assignment(cost_matrix)

        # Recover global i,j which is usable to get the mappings
        for i, j in zip(row_ind, col_ind):

            if cost_matrix[i, j] != 1e9: #valid match
                ref_id = rows[i]
                tgt_id = cols[j]

                if (ref_id, tgt_id) not in mappings:
                    print(f"Skipping assignment ({i}, {j}) — not in valid mappings")
                    continue

                mapping = mappings[(ref_id, tgt_id)]
                self.chain_mappings[ref_id] = mapping

                self.matched_ref_chains.add(ref_id)
                self.matched_tgt_chains.add(tgt_id)
                print(f"Matched {ref_id} → {tgt_id} with score {scores[(ref_id, tgt_id)]:.4f}")

    def map_chains_explicit(self, explicit_chain_mapping):
        for ref_chain_id, (tgt_chain_id, type)in explicit_chain_mapping.items():
            ref_chain = self.ref_structure[0][ref_chain_id]
            tgt_chain = self.tgt_structure[0][tgt_chain_id]

            self.matched_ref_chains.add(ref_chain_id)
            self.matched_tgt_chains.add(tgt_chain_id)

            ref_seq = extract_seq_from_chain(ref_chain)
            tgt_seq = extract_seq_from_chain(tgt_chain)

            alignment = self.aligner.align(ref_seq, tgt_seq)[0]

            self.chain_mappings[ref_chain_id] = ChainMapper(ref_chain, ref_seq, tgt_chain, tgt_seq, alignment, type)

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
    
    def get_selected_mapping(self, selected_chains=None, protein_mode="CA", nucleic_mode="C1'"):
        """
        Extract a mapping from the ChainMapper Objects based on a selection of chains
        """
        #dictionary to automatically choose correct mode
        type_to_mode = {
            "protein":protein_mode,
            "dna":nucleic_mode,
            "rna":nucleic_mode
        }

        # current # of residues that have been mapped so far
        total_residues = 0
        chain_info = []

        for chain_id, cm in self.chain_mappings.items():
            if selected_chains and chain_id not in selected_chains:
                print(f"[INFO] Skipping chain {chain_id} (not in selected_chains)")
                continue
            
            #get residue lists for an individual chain pairing in the structure
            aligned_ref_coords, aligned_tgt_coords, chain_index_map = cm.get_aligned_coord_lists(type_to_mode[cm.type])

            n_residues = len(aligned_ref_coords)
            if n_residues == 0:
                continue
            
            chain_info.append((cm, aligned_ref_coords, aligned_tgt_coords, chain_index_map, total_residues))
            total_residues += n_residues

        if total_residues == 0:
            raise ValueError("No aligned coordinates found for selected chains.")

        # pre-allocate arrays to store whole structure mapping
        coords_ref = np.zeros((total_residues, 3), dtype=float)
        coords_tgt = np.zeros((total_residues, 3), dtype=float)
        index_map_global = {}
        res_id_map_global = {}

        for cm, ref_coords, tgt_coords, index_map, offset in chain_info:
            n_res = len(ref_coords)
            coords_ref[offset:offset+n_res] = ref_coords
            coords_tgt[offset:offset+n_res] = tgt_coords

            
            for ref_key, local_idx in index_map.items():
                index_map_global[ref_key] = offset + local_idx

            
            res_id_map_global.update(cm.res_id_mapping)

        return coords_ref, coords_tgt, index_map_global, res_id_map_global


class EnsembleMapper:
    """
    This is an Ensemble Mapper Class designed to handle mutiple structure mapping objects
    As currently implemented, it will store relevant mappings for selected chains calculated via the get_selected_global_coords class method
    All coordinates are aligned to a single set of reference cordinates and thus corresponding coordinates can all be found using the same index
    """

    def __init__(self, ref_structure, protein_aligner, nucleotide_aligner):
        self.ref_structure=ref_structure
        self.protein_aligner = protein_aligner
        self.nucleotide_aligner = nucleotide_aligner
        self.coords_ref = None
        self.structure_mappings = {} # Structure Name -> StructureMapper
        self.global_index_mapping = {} # (ChainID, ResID) -> index (for use with aligned tgt coords)
        self.res_id_mappings = {} # Structure Name -> (ref(ChainID, ResID) -> tgt(ChainID, ResID))
        self.coords_targets_dict= {} # Structure Name -> aligned_tgt_coords
    
    def add_structure(self, tgt_structure_name, tgt_structure, threshold, explicit_mapping=None):
        structure_mapping = StructureMapper(self.ref_structure, tgt_structure, self.protein_aligner, self.nucleotide_aligner)
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
                #Add keys of res_id map from chain mapping (matched residues)
                mapped_ref_residues.update(structure_mapping.chain_mappings[chain_id].res_id_mapping.keys())
            
            mapped_residues_sets.append(mapped_ref_residues)
        
        if not mapped_residues_sets:
            return set()
        
        common_ref_residues = set.intersection(*mapped_residues_sets)

        
        return common_ref_residues
    
    def set_selected_global_coords(self, selected_chains=None, protein_mode="CA", nucleic_mode="C1'"):
        """
        Precompute global coordinates for reference and targets, ensuring consistent ordering,
        and handling missing coordinates safely (no NaNs).
        """
        if not self.structure_mappings:
            raise ValueError("No target structures added yet.")

        if selected_chains is None:
            selected_chains = [chain.id for chain in self.ref_structure.get_chains()]

        #get residues mapped in all structures
        common_ref_residues = self.get_common_ref_residues(selected_chains)
        if not common_ref_residues:
            raise ValueError("No residues mapped in all structures.")

        coords_ref_per_structure = {}
        coords_tgt_per_structure = {}
        index_maps_per_structure = {}
        res_id_maps_per_structure = {}

        for name, structure_mapping in self.structure_mappings.items():
            coords_ref, coords_tgt, index_map, res_id_map = structure_mapping.get_selected_mapping(
                selected_chains, protein_mode, nucleic_mode
            )

            # only keep residues that actually have coordinates
            valid_residues = set(index_map.keys()) & common_ref_residues
            if not valid_residues:
                continue

            # filter coordinate arrays to include only valid residues
            coords_ref_filtered = np.array([coords_ref[index_map[r]] for r in valid_residues])
            coords_tgt_filtered = np.array([coords_tgt[index_map[r]] for r in valid_residues])

            coords_ref_per_structure[name] = coords_ref_filtered
            coords_tgt_per_structure[name] = coords_tgt_filtered

            # map only valid residues
            index_maps_per_structure[name] = {r: i for i, r in enumerate(valid_residues)}
            res_id_maps_per_structure[name] = {r: res_id_map[r] for r in valid_residues}

        if not coords_ref_per_structure:
            raise ValueError("No valid residues with coordinates in any structure.")

        # determine the set of residues present in all structures
        all_valid_residue_sets = [set(idx_map.keys()) for idx_map in index_maps_per_structure.values()]
        global_residues = set.intersection(*all_valid_residue_sets)
        if not global_residues:
            raise ValueError("No residues with coordinates in all structures.")

        sorted_residues = sorted(global_residues)

        # get reference structure info since it will be reused
        first_name = next(iter(coords_ref_per_structure.keys()))
        first_index_map = index_maps_per_structure[first_name]
        first_coords_ref = coords_ref_per_structure[first_name]

        self.coords_ref = np.array([first_coords_ref[first_index_map[r]] for r in sorted_residues])
        self.global_index_mapping = {res_id: idx for idx, res_id in enumerate(sorted_residues)}

        #build target coordinates per structure
        self.coords_targets_dict.clear()
        self.res_id_mappings.clear()

        for name in self.structure_mappings.keys():
            index_map_tgt = index_maps_per_structure[name]
            coords_tgt_all = coords_tgt_per_structure[name]
            res_id_map = res_id_maps_per_structure[name]

            # put into global index scheme
            self.coords_targets_dict[name] = np.array([coords_tgt_all[index_map_tgt[r]] for r in sorted_residues])
            self.res_id_mappings[name] = {r: res_id_map[r] for r in sorted_residues}

                                

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
    
    def _compute_diffs(self, ref_coords, tgt_coords):
        diff_vecs = ref_coords - tgt_coords
        distances = np.linalg.norm(diff_vecs, axis=1)
        rmsd = np.sqrt(np.mean(np.sum(diff_vecs**2, axis=1)))
        return diff_vecs, distances, rmsd
    
    def calc_rmsds(self):
        """
        Calculate global RMSD for each structure.
        """
        rmsds = {}
        ref_coords = self.coords_ref

        for structure_name, tgt_coords in self.coords_targets_dict.items():
           _,_,rmsd = self._compute_diffs(ref_coords, tgt_coords)
           rmsds[structure_name] = rmsd
        return rmsds

    def calc_movement_dfs(self):
        """
        Calculate movement dfs for each structure mapped to the reference, by taking the common residues and finding their difference with corresponding 
        residues in the reference 
        """
        movement_dfs = {}
        ref_coords = self.coords_ref

        for structure_name, tgt_coords in self.coords_targets_dict.items():
            
            tgt_coords = self.coords_targets_dict[structure_name]

            diff_vecs, distances,_ = self._compute_diffs(ref_coords, tgt_coords)

            ref_ids, tgt_ids = zip(*self.res_id_mappings[structure_name].items())

            ref_strs = [f"{chain}-{res[1]}{res[2]}" for chain, res in ref_ids]
            tgt_strs = [f"{chain}-{res[1]}{res[2]}" for chain, res in tgt_ids]
            
            df = pd.DataFrame({
            "ChainID_Resnum1": ref_strs,
            "ChainID_Resnum2": tgt_strs,
            "Coord1": np.array(ref_coords).tolist(),
            "Coord2": np.array(tgt_coords).tolist(),
            "Diff_Vec": np.array(diff_vecs).tolist(),
            "Distance": distances
            })

            movement_dfs[structure_name] = df
            
        return movement_dfs


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
            res_copy = res
            key = (chain.id, res.id)
            #No residue match or not filtering residues
            if matched_residues is None or key not in matched_residues:
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