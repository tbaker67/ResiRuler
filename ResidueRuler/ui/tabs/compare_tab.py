"""Compare tab for analyzing distance differences between structures."""
import os

import streamlit as st

from src.resiruler.viz.plotting import (
    plot_all_matrices_ensemble,
    plot_comparison_with_contact_filter,
    plot_contacts_gained,
    plot_contacts_lost,
)
from ui.widgets.utils import (
    aligner_ui,
    chain_selector_ui,
    create_ensemble_mapper,
    display_chain_pair_selector,
    distance_threshold_ui,
    get_chain_mappings_for_targets,
    get_measurement_mode,
    get_threshold,
    load_structure_if_new,
    load_structures_if_new,
    show_alignments,
)

def show_compare_tab():
    st.header("Compare Distances within two structures")

    
    ref_cif = st.file_uploader("Upload Aligned Reference CIF", type=["cif"], key="reference1")
    tgt_cifs = st.file_uploader("Upload Aligned Target CIFs", type=["cif"], key="tgts", accept_multiple_files=True)

    
    ref_structure = load_structure_if_new(ref_cif, "compare_name1", "compare_structure1")
    tgt_structures = load_structures_if_new(tgt_cifs, "compare_name2", "compare_structure2")

   
    if ref_structure and tgt_structures:
        ref_chains = [ref_chain.id for ref_chain in ref_structure[0].get_chains()]
        chain_mappings = get_chain_mappings_for_targets(tgt_structures,ref_chains, key="compare_mappings")

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

    selected_chains = chain_selector_ui(ref_structure, "Select Chains in reference to compare")

    protein_mode, nucleic_mode = get_measurement_mode(key="compare_measurement_mode")

    st.subheader("Display Settings")
    
    enable_threshold, lower_threshold, upper_threshold = distance_threshold_ui(
        key_prefix="compare_contact"
    )
    
    contact_threshold = upper_threshold if enable_threshold else None

    if st.button("Compare") and st.session_state.mapper:
        st.session_state.mapper.set_selected_global_coords(selected_chains, protein_mode=protein_mode, nucleic_mode=nucleic_mode)
        ref_dm, tgt_dms_dict, compare_dms_dict = st.session_state.mapper.calc_matrices()
        st.session_state.ref_dm = ref_dm
        st.session_state.tgt_dms_dict = tgt_dms_dict
        st.session_state.compare_dms_dict = compare_dms_dict
    if "ref_dm" in st.session_state and st.session_state.ref_dm is not None:

        target_names = list(st.session_state.tgt_dms_dict.keys())
        selected_target = st.selectbox(
            "Select target structure to display",
            target_names,
            key="selected_target_display"
        )

        residue_count = len(st.session_state.ref_dm.index_map)
        if residue_count > 2000:
            st.info(f"Large structure detected ({residue_count:,} residues). Chain-pair view recommended for better performance.")
        
        display_mode = st.radio(
            "Display Mode",
            ["Full Contact Maps", "Chain Pair View"],
            index=1 if residue_count > 2000 else 0,
            horizontal=True,
            help="Full view shows entire contact map. Chain pair view lets you select specific chain interactions for large structures.",
            key="display_mode"
        )
        
        if display_mode == "Full Contact Maps":
            if "ref_fig" not in st.session_state or st.session_state.ref_fig is None:
                with st.spinner("Generating full contact maps (this may take a moment for large structures)..."):
                    ref_fig, tgt_figs_dict, compare_figs_dict = plot_all_matrices_ensemble(
                        st.session_state.ref_dm, st.session_state.tgt_dms_dict, 
                        st.session_state.compare_dms_dict, 
                        lower_threshold=lower_threshold, upper_threshold=upper_threshold
                    )
                    st.session_state.ref_fig = ref_fig
                    st.session_state.tgt_figs_dict = tgt_figs_dict
                    st.session_state.compare_figs_dict = compare_figs_dict
            
            st.plotly_chart(st.session_state.ref_fig, use_container_width=False)
            st.plotly_chart(st.session_state.tgt_figs_dict[selected_target], use_container_width=False)
            
            if enable_threshold and contact_threshold is not None:
                compare_dm = st.session_state.compare_dms_dict[selected_target]
                
                st.markdown("#### Distance Changes in Shared Contacts")
                shared_fig = plot_comparison_with_contact_filter(
                    compare_dm, contact_threshold=contact_threshold,
                    title=f"ΔDistance (shared contacts < {contact_threshold}Å)"
                )
                st.plotly_chart(shared_fig, use_container_width=False)
                
                col1, col2 = st.columns(2)
                with col1:
                    st.markdown("#### Contacts Gained in Target")
                    gained_fig = plot_contacts_gained(
                        compare_dm, contact_threshold=contact_threshold
                    )
                    gained_fig.update_layout(width=600, height=600)
                    st.plotly_chart(gained_fig, use_container_width=True)
                    
                with col2:
                    st.markdown("#### Contacts Lost in Target")
                    lost_fig = plot_contacts_lost(
                        compare_dm, contact_threshold=contact_threshold
                    )
                    lost_fig.update_layout(width=600, height=600)
                    st.plotly_chart(lost_fig, use_container_width=True)
            else:
                st.plotly_chart(st.session_state.compare_figs_dict[selected_target], use_container_width=False)
            
            st.caption("To refresh plots with new threshold settings, click 'Compare' again.")
        else:
            display_chain_pair_selector(
                st.session_state.ref_dm,
                st.session_state.tgt_dms_dict[selected_target],
                st.session_state.compare_dms_dict[selected_target],
                selected_target,
                lower_threshold, upper_threshold,
                contact_threshold=contact_threshold if enable_threshold else None,
            )

        st.subheader("Export Data")
        
        compare_dm = st.session_state.compare_dms_dict[selected_target]
        total_pairs = len(compare_dm.index_map) * (len(compare_dm.index_map) - 1) // 2
        
        st.caption(f"Full dataset: {total_pairs:,} residue pairs")
        
        output_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "output")
        os.makedirs(output_dir, exist_ok=True)
        
        csv_filename = f"comparison_{selected_target}.csv"
        csv_filepath = os.path.join(output_dir, csv_filename)
        
        if st.button("Export Full Comparison to CSV", key="export_csv"):
            with st.spinner(f"Exporting {total_pairs:,} pairs to CSV (streaming to avoid memory issues)..."):
                rows_written = compare_dm.export_to_csv_streaming(csv_filepath)
                st.success(f"Exported {rows_written:,} rows to: `{csv_filepath}`")
                
                # For smaller files direct download
                if total_pairs < 500000:
                    with open(csv_filepath, 'r') as f:
                        csv_data = f.read()
                    st.download_button(
                        label="Download CSV",
                        data=csv_data,
                        file_name=csv_filename,
                        mime="text/csv",
                        key="download_csv"
                    )
                else:
                    st.info(f"File is large ({total_pairs:,} pairs). Download from: `{csv_filepath}`")
