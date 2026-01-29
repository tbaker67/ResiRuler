import streamlit as st
from ui_components.utils import chain_selector_ui, load_structure_if_new, get_threshold, load_structures_if_new, get_chain_mappings_for_targets, create_ensemble_mapper, aligner_ui, show_alignments,get_measurement_mode
from src.resiruler.plotting import plot_distance_difference, plot_interactive_contact_map, plot_all_matrices_ensemble
import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices


def show_compare_tab():
    st.header("Compare Distances within two structures")

    
    ref_cif = st.file_uploader("Upload Aligned Reference CIF", type=["cif"], key="reference1")
    tgt_cifs = st.file_uploader("Upload Aligned Target CIFs", type=["cif"], key="tgts", accept_multiple_files=True)

    
    ref_structure = load_structure_if_new(ref_cif, "compare_name1", "compare_structure1")
    tgt_structures = load_structures_if_new(tgt_cifs, "compare_name2", "compare_structure2")

   
    if ref_structure and tgt_structures:
        ref_chains = [ref_chain.id for ref_chain in ref_structure[0].get_chains()]
        chain_mappings = get_chain_mappings_for_targets(tgt_structures,ref_chains, key="compare_mappings")

       
    #get threshold and do alignments
    
    
    st.session_state.setdefault("mapper", None)
    with st.container():
        col1, col2 = st.columns(2)
        with col1:
            st.subheader("Protein Pairwise Aligner Settings")
            protein_aligner = aligner_ui(protein=True, key_prefix="protein compare aligner")
        with col2:
            st.subheader("Nucleotide Pairwise Aligner Settings")
            nucleotide_aligner = aligner_ui(protein=False, key_prefix="nucleotide compare aligner")

    pct_id_threshold = get_threshold("Set a Minimum Percent Identity Threshold for Matching Chains Together", "95.0","compare_pct_id")

    if st.button("Map Chains", key = "map compare chains"):
        st.session_state.mapper = create_ensemble_mapper(ref_structure, tgt_structures, chain_mappings, pct_id_threshold, protein_aligner, nucleotide_aligner)
    
    if st.session_state.mapper is not None:
        show_alignments(st.session_state.mapper, key="compare_alignment")
    ##TODO: Allow for filtering to only matched chains across all
    selected_chains = chain_selector_ui(ref_structure, "Select Chains in reference to compare")

    protein_mode, nucleic_mode = get_measurement_mode(key="compare_measurement_mode")

    if st.button("Compare") and st.session_state.mapper:
        st.session_state.mapper.set_selected_global_coords(selected_chains, protein_mode=protein_mode, nucleic_mode=nucleic_mode)
        ref_dm, tgt_dms_dict, compare_dms_dict = st.session_state.mapper.calc_matrices()
        
        ref_fig, tgt_figs_dict, compare_figs_dict = plot_all_matrices_ensemble(
            ref_dm, tgt_dms_dict, compare_dms_dict, lower_threshold=-10000, upper_threshold=10000
        )

        
        st.session_state.ref_dm = ref_dm
        st.session_state.tgt_dms_dict = tgt_dms_dict
        st.session_state.compare_dms_dict = compare_dms_dict
        st.session_state.ref_fig = ref_fig
        st.session_state.tgt_figs_dict = tgt_figs_dict
        st.session_state.compare_figs_dict = compare_figs_dict
        st.session_state.compare_dfs = {
            tgt_name: compare_dms_dict[tgt_name].convert_to_df()
            for tgt_name in compare_dms_dict
        }

   
    if "ref_fig" in st.session_state and st.session_state.ref_fig:
        target_names = list(st.session_state.tgt_dms_dict.keys())
        selected_target = st.selectbox(
            "Select target structure to display",
            target_names,
            key="selected_target_display"
        )

        #show plots
        st.plotly_chart(st.session_state.ref_fig, use_container_width=False)
        st.plotly_chart(st.session_state.tgt_figs_dict[selected_target], use_container_width=False)
        st.plotly_chart(st.session_state.compare_figs_dict[selected_target], use_container_width=False)

        compare_df = st.session_state.compare_dfs[selected_target]
        st.dataframe(
            compare_df.drop(columns=['Coord1_ref', 'Coord2_ref', 'Coord1_tgt', 'Coord2_tgt']),
            use_container_width=True
        )