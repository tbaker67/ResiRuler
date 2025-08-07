import Bio
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.PDB import MMCIFParser, MMCIFIO, Structure, Model
from io import StringIO
from .structure_parsing import extract_res_from_chain, extract_seq_from_chain
from .distance_calc import DistanceMatrix, CompareDistanceMatrix
import numpy as np
from scipy.optimize import linear_sum_assignment

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

        self.index_map = {} # (Chain_ID, Resnum) -> index
        self.res_id_map = {} # (Chain_ID, Resnum) -> (Chain_ID, Resnum)

        self.aligned_ref_coords = []
        self.aligned_tgt_coords = []

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

    def calc_aligned_coords(self):
        """
        Gets the aligned coord lists, stores back in the mapper, the aligned coord list will only contain coords for residues in both the reference and target structures
        """
        ref_res = extract_res_from_chain(self.ref_chain)
        tgt_res = extract_res_from_chain(self.tgt_chain)

        aligned_ref_seq = str(self.alignment[0])
        aligned_tgt_seq = str(self.alignment[1])

        assert len(aligned_ref_seq) == len(aligned_tgt_seq)

        aligned_ref_coords = []
        aligned_tgt_coords = []
        res_id_mapping = {}
        index_map = {}

        idx_ref = 0
        idx_tgt = 0
        coord_index = 0

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
                if "CA" in a_res and "CA" in b_res:
                    aligned_ref_coords.append(a_res["CA"].get_coord())
                    aligned_tgt_coords.append(b_res["CA"].get_coord())
                    ref_id = (self.ref_chain.id, a_res.id[1])
                    tgt_id = (self.tgt_chain.id, b_res.id[1])
                    res_id_mapping[ref_id] = tgt_id
                    index_map[ref_id] = coord_index
                    coord_index += 1

        self.aligned_ref_coords = aligned_ref_coords
        self.aligned_tgt_coords = aligned_tgt_coords
        self.index_map = index_map
        self.res_id_map = res_id_mapping

    def get_ref_coord(self, chain, resnum):
        """
        Gets coordinate of residue in reference structure based on (Chain, Resnum) identifier
        """
        return self.aligned_ref_coords[self.index_map[(chain, resnum)]]

    def get_tgt_coord(self, chain,resnum):
        """
        Gets coordinate of residue in reference structure based on (Chain, Resnum) identifier
        We should use the reference (Chain, Resnum) identifier as that is what is in the index map, and the aligned coordinates ensure we get a nice 1-to-1 mapping
        """
        return self.aligned_tgt_coords[self.index_map[(chain,resnum)]]



        

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
        self.ref_structure = ref_structure
        self.tgt_structure = tgt_structure
        self.aligner = self.aligner = aligner or self._default_aligner()
        self.chain_mappings = {} # ref_chain_id -> ChainMapping object 

        self.matched_ref_chains = set()
        self.matched_tgt_chains = set()

    def _default_aligner(self):
        '''
        Pairwise alignment that uses the Blosum62 substition matrix and Needleman-wunsch algorithm
        '''
        aligner = PairwiseAligner()
        aligner.mode = "global"
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        aligner.left_open_gap_score = 1.0
        return aligner

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
            ref_seq = extract_seq_from_chain(ref_chain)

            for j, tgt_chain in enumerate(tgt_chains):
                tgt_seq = extract_seq_from_chain(tgt_chain)

                alignment = self.aligner.align(ref_seq, tgt_seq)[0]
                potential_map = ChainMapper(ref_chain, ref_seq, tgt_chain, tgt_seq, alignment)

                percent_identity = potential_map.calc_percent_identity()
                #print(f"Testing chains {ref_chain.id} vs {tgt_chain.id} — %ID: {percent_identity:.2f}")
                

                if percent_identity < threshold:
                    continue

                potential_map.calc_aligned_coords()
                rmsd = potential_map.calc_rmsd(
                    potential_map.aligned_ref_coords,
                    potential_map.aligned_tgt_coords
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
        coords_ref = []
        coords_tgt = []
        index_map = {} 
        res_id_map = {}

        coord_index = 0
        for chain_id, cm in self.chain_mappings.items():
            if selected_chains and chain_id not in selected_chains:
                print(f"[INFO] Skipping chain {chain_id} (not in selected_chains)")
                continue

            for key, idx in cm.index_map.items():
                try:
                    ref_coord = cm.aligned_ref_coords[idx]
                    tgt_coord = cm.aligned_tgt_coords[idx]
                except IndexError:
                    print(f"[WARN] Index {idx} out of bounds in chain {chain_id}")
                    continue

                coords_ref.append(ref_coord)
                coords_tgt.append(tgt_coord)
                index_map[key] = coord_index
                coord_index += 1
            
            res_id_map = res_id_map | cm.res_id_map

        if len(coords_ref) == 0:
            raise ValueError("No aligned coordinates were found. Check selected chains and mapping.")

        coords_ref = np.array(coords_ref)
        coords_tgt = np.array(coords_tgt)

        if coords_ref.ndim != 2 or coords_ref.shape[1] != 3:
            raise ValueError(f"ref coords shape: {coords_ref.shape} — expected (N, 3)")

        ref_dm = DistanceMatrix(coords_ref, index_map)
        tgt_dm = DistanceMatrix(coords_tgt, index_map)
        return ref_dm, tgt_dm, CompareDistanceMatrix(ref_dm, tgt_dm, res_id_map)


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
            key = (chain.id, res.id[1])
            #No residue match or not filtering residues
            if matched_residues is not None and key not in matched_residues:
                continue
            new_chain.add(res)

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
        for ref_res_id, tgt_res_id in cm.res_id_map.items():
            matched_ref_residues.add((ref_res_id)) # (chain, resnum)
            
            matched_tgt_residues.add((tgt_res_id)) # (chain, resnum)

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