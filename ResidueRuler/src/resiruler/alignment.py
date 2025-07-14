import Bio
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.PDB import MMCIFParser
from Bio.SeqUtils import seq1
from structure_parsing import load_structure
from distance_calc import DistanceMatrix, CompareDistanceMatrix

class ResidueMapper:
    def __init__(self, reference_chains, target_chains, aligner=None):
        self.aligner = aligner or self.default_aligner()
        self.reference_chains = reference_chains
        self.target_chains = target_chains

        self.mapping = {}      # (ref_chain_id, resnum) -> (tgt_chain_id, resnum)
        self.ref_coords = []   # List of CA coords from reference residues (joint)
        self.tgt_coords = []   # List of CA coords from target residues (joint)
        self.index_map = {}    # Map residue keys to global index in coords
        self.aligned_seqA = "" # concatenated aligned seq for all chains ref
        self.aligned_seqB = "" # concatenated aligned seq for all chains target

        self._build_joint_mapping()

    def default_aligner(self):
        from Bio.Align import PairwiseAligner, substitution_matrices
        aligner = PairwiseAligner()
        aligner.mode = "global"
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        return aligner

    def _extract_seq_and_res(self, chain):
        from Bio.SeqUtils import seq1
        seq = ""
        res_list = [res for res in chain.get_residues() if res.id[0] == ' ']  # standard AA residues
        for res in res_list:
            try:
                seq += seq1(res.get_resname())
            except KeyError:
                seq += "X"
        return seq, res_list

    def _align_sequences(self, seqA, seqB):
        alignment = self.aligner.align(seqA, seqB)[0]
        return alignment

    def _align_seq_and_res(self, seqA, seqB, res_listA, res_listB, alignment):
        aligned_seqA = ""
        aligned_seqB = ""
        mapping = {}

        aligned_blockA, aligned_blockB = alignment.aligned[0], alignment.aligned[1]
        idxA, idxB = 0, 0

        for (startA, endA), (startB, endB) in zip(aligned_blockA, aligned_blockB):
            while idxA < startA:
                aligned_seqA += seqA[idxA]
                aligned_seqB += "-"
                idxA += 1
            while idxB < startB:
                aligned_seqA += "-"
                aligned_seqB += seqB[idxB]
                idxB += 1

            for i in range(endA - startA):
                aligned_seqA += seqA[startA + i]
                aligned_seqB += seqB[startB + i]
                mapping[res_listA[startA + i]] = res_listB[startB + i]

            idxA, idxB = endA, endB

        while idxA < len(seqA):
            aligned_seqA += seqA[idxA]
            aligned_seqB += "-"
            idxA += 1
        while idxB < len(seqB):
            aligned_seqA += "-"
            aligned_seqB += seqB[idxB]
            idxB += 1

        return aligned_seqA, aligned_seqB, mapping

    def _build_joint_mapping(self):
        index = 0

        for chainA in self.reference_chains:
            seqA, res_listA = self._extract_seq_and_res(chainA)

            # Find the best matching chain in target (simple max score)
            best_score = -float('inf')
            best_chainB = None
            best_seqB = None
            best_res_listB = None
            best_alignment = None

            for chainB in self.target_chains:
                seqB, res_listB = self._extract_seq_and_res(chainB)
                alignment = self._align_sequences(seqA, seqB)
                if alignment.score > best_score:
                    best_score = alignment.score
                    best_chainB = chainB
                    best_seqB = seqB
                    best_res_listB = res_listB
                    best_alignment = alignment

            aligned_seqA, aligned_seqB, mapping = self._align_seq_and_res(seqA, best_seqB, res_listA, best_res_listB, best_alignment)

            self.aligned_seqA += aligned_seqA
            self.aligned_seqB += aligned_seqB

            for ref_res, tgt_res in mapping.items():
                key_ref = (chainA.id, ref_res.id[1])
                key_tgt = (best_chainB.id, tgt_res.id[1])

                # Only consider residues with CA atoms present
                if "CA" in ref_res and "CA" in tgt_res:
                    self.mapping[key_ref] = key_tgt
                    self.ref_coords.append(ref_res["CA"].get_coord())
                    self.tgt_coords.append(tgt_res["CA"].get_coord())
                    self.index_map[key_ref] = index
                    index += 1


        ##TODO ADD IN FUNCTIONALITY TO BUILD A DISTANCE MATRIX CLASS WHICH WILL STORE THE DISTANCE MATRIX AS WELL AS HOW TO INDEX IT
        #The eventual goal might be 

        ##WE can make it so to create comparison matrix, just make a (reference_coords,reference_coors) and a (target_coords, target_coords) matrix and then subtract

        ##For a Movement matrix we would just have (reference_coords, target_coords matrix) It actually won't even need to be a matrix because we're just using pairwise


class StructureMapper:

    def __init__(self,reference_structure, chains=None):
        self.reference_structure=reference_structure
        if chains is None:
            self.ref_chains = self.get_chains(self.reference_structure)
        else:
            self.ref_chains = [chain for chain in self.get_chains(self.reference_structure) if chain.id in chains]
        
       
        
        self.chain_mappers = {}

    def get_chains(self, structure):
        chains = []
        for chain in structure.get_chains():
            chains.append(chain)
        return chains
    
    def map_target_structure(self, target_structure):
    
    # Extract all chains from target structure
        target_chains = [chain for chain in target_structure.get_chains()]
    
    # Create a single joint ResidueMapper for the whole structure
        self.structure_mapper = ResidueMapper(self.ref_chains, target_chains)
    
        return self.structure_mapper

    def map_ensemble(self,target_ensemble):
        for target_id, target_structure in target_ensemble:
            mapped_target = self.map_target_structure(target_structure)



def calculate_percent_identity(seq1, seq2):
    matches = sum(a == b for a, b in zip(seq1, seq2) if a != '-' and b != '-')
    aligned = sum(1 for a, b in zip(seq1, seq2) if a != '-' and b != '-')
    return 100 * matches / aligned if aligned > 0 else 0.0

reference = load_structure('/Users/tbaker/Desktop/ResiRuler/aligned1.cif')
target = load_structure('/Users/tbaker/Desktop/ResiRuler/aligned2.cif')
target2 = load_structure('/Users/tbaker/Desktop/ResiRuler/fold_cagx_monomer_model_0.cif')

mapper = StructureMapper(reference, ["AX","AY","Am"])

# Create joint residue mapping for whole target structure
residue_mapper = mapper.map_target_structure(target)

# Build coordinate lists and mappings
ref_coords, target_coords, index_map, residue_pair_map = residue_mapper.ref_coords, residue_mapper.tgt_coords, residue_mapper.index_map, residue_mapper.mapping

print(f"Number of matched residues (ref coords): {len(ref_coords)}")
print(f"Number of matched residues (target coords): {len(target_coords)}")
print(f"Number of entries in index map: {len(index_map)}")
print(f"Number of residue pairs in mapping: {len(residue_pair_map)}")

reference_mat = DistanceMatrix(ref_coords, index_map)
print(reference_mat.get_distance(("AX", 514), ("AY", 1898)))
target_mat = DistanceMatrix(target_coords, index_map)
print(target_mat.get_distance(("AX", 514), ("AY", 1898)))
compare_mat = CompareDistanceMatrix(reference_mat, target_mat)
print(compare_mat.get_distance_diff(("AX", 514), ("AY", 1898)))


####THESE WILL BE OFF, CURRENTLY THE WAY SEQUENCE MATCHING IS WORKING MEAN WE SHOULD ONLY DO THIS WITH MONOMERS
print(reference_mat.get_distance(("AX", 365), ("Am", 291)))
print(target_mat.get_distance(("AX", 365), ("Am", 291)))
print(compare_mat.get_distance_diff(("AX", 365), ("Am", 291)))