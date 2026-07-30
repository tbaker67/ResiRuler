"""Compare tab for analyzing distance differences between structures."""

import os
from io import BytesIO

import streamlit as st
from src.resiruler.viz.plotting import (
    plot_all_matrices_ensemble,
    plot_comparison_with_contact_filter,
    plot_contacts_gained,
    plot_contacts_lost,
    plot_chain_average_map
)

from ui.widgets.utils import (
    chain_selector_ui,
    color_scale_ui,
    create_ensemble_mapper,
    distance_threshold_ui,
    full_aligner_ui,
    get_chain_mappings_for_targets,
    get_measurement_mode,
    get_threshold,
    load_structure_if_new,
    load_structures_if_new,
    reset_downstream,
    show_alignments,
)


def add_svg_download(fig, filename):
    svg_buffer = BytesIO()
    fig.write_image(svg_buffer, format="svg")

    st.download_button(
        label="Download SVG",
        data=svg_buffer.getvalue(),
        file_name=filename,
        mime="image/svg+xml",
        key=f"download_{filename}",
    )


def show_compare_tab():
    st.header("Compare Distances Within Two Structures")

    ref_cif = st.file_uploader(
        "Upload Aligned Reference CIF", type=["cif"], key="reference1"
    )
    tgt_cifs = st.file_uploader(
        "Upload Aligned Target CIFs",
        type=["cif"],
        key="tgts",
        accept_multiple_files=True,
    )

    ref_structure = load_structure_if_new(
        ref_cif, "compare_name1", "compare_structure1"
    )
    tgt_structures = load_structures_if_new(
        tgt_cifs, "compare_name2", "compare_structure2"
    )

    if ref_structure and tgt_structures:
        ref_chains = [ref_chain.id for ref_chain in ref_structure[0].get_chains()]
        chain_mappings = get_chain_mappings_for_targets(
            tgt_structures, ref_chains, key="compare_mappings"
        )
        print(chain_mappings)

    st.session_state.setdefault("mapper", None)
    protein_aligner, nucleotide_aligner = full_aligner_ui(key="compare")

    pct_id_threshold = get_threshold(
        "Set a Minimum Percent Identity Threshold for Matching Chains Together",
        "95.0",
        "compare_pct_id",
    )

    if st.button("Map Chains", key="map compare chains"):
        reset_downstream(
            "mapper",
            "ref_dm",
            "tgt_dms_dict",
            "compare_dms_dict",
            "small_ref_dm",
            "small_tgt_dms_dict",
            "small_compare_dms_dict",
            "ref_fig",
            "tgt_figs_dict",
            "compare_figs_dict",
            "_last_selected_chains",
        )
        st.session_state.mapper = create_ensemble_mapper(
            ref_structure,
            tgt_structures,
            chain_mappings,
            pct_id_threshold,
            protein_aligner,
            nucleotide_aligner,
        )

    if st.session_state.get("mapper") is not None:
        show_alignments(st.session_state.mapper, key="compare_alignment")

    protein_mode, nucleic_mode = get_measurement_mode(key="compare_measurement_mode")

    if st.button("Compare") and st.session_state.get("mapper"):
        # the previous full-resolution ensemble matrices (and any chain-filtered
        # "small" copies/figures derived from them) are about to be superseded
        reset_downstream(
            "ref_dm",
            "tgt_dms_dict",
            "compare_dms_dict",
            "small_ref_dm",
            "small_tgt_dms_dict",
            "small_compare_dms_dict",
            "ref_fig",
            "tgt_figs_dict",
            "compare_figs_dict",
            "_last_selected_chains",
        )
        st.session_state.mapper.set_selected_global_coords(
            selected_chains=None,
            protein_mode=protein_mode,
            nucleic_mode=nucleic_mode,
        )
        ref_dm, tgt_dms_dict, compare_dms_dict = st.session_state.mapper.calc_matrices()
        st.session_state.ref_dm = ref_dm
        st.session_state.tgt_dms_dict = tgt_dms_dict
        st.session_state.compare_dms_dict = compare_dms_dict

    if "ref_dm" in st.session_state and st.session_state.ref_dm is not None:
        target_names = list(st.session_state.tgt_dms_dict.keys())
        selected_target = st.selectbox(
            "Select Target Structure to Display",
            target_names,
            key="selected_target_display",
        )

        compare_dm = st.session_state.compare_dms_dict[selected_target]

        # Reserve the plot's position first, render its color-scale control later
        overview_slot = st.empty()
        chain_avg_scale = color_scale_ui(
            purpose="diverging",
            key_prefix="compare_chainavg_scale",
            default_colorscale="RdBu_r",
            min= float(-max(abs(v.mat).max() for v in st.session_state.compare_dms_dict.values())),
            max= float(max(abs(v.mat).max() for v in st.session_state.compare_dms_dict.values())),
        )
        overview_plot = plot_chain_average_map(
            compare_dm,
            colorscale=chain_avg_scale["colorscale"],
            vmin=chain_avg_scale["vmin"],
            vmax=chain_avg_scale["vmax"],
        )
        overview_slot.plotly_chart(overview_plot)

        total_pairs = len(compare_dm.index_map) * (len(compare_dm.index_map) - 1) // 2
        st.caption(f"Full dataset: {total_pairs:,} residue pairs")

        output_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "output")
        os.makedirs(output_dir, exist_ok=True)

        csv_filename = f"comparison_{selected_target}.csv"
        csv_filepath = os.path.join(output_dir, csv_filename)

        if st.button("Export Full Comparison to CSV", key="export_csv"):
            with st.spinner(f"Exporting {total_pairs:,} pairs to CSV"):
                rows_written = compare_dm.export_to_csv_streaming(csv_filepath)
                st.success(f"Exported {rows_written:,} rows to: `{csv_filepath}`")

        selected_chains = chain_selector_ui(
            ref_structure,
            "Select Chains in Reference to View in Full Resolution",
            key_prefix="compare_small_chain_selector",
        )
        st.subheader("Display Settings")
        enable_threshold, lower_threshold, upper_threshold = distance_threshold_ui(
            key_prefix="compare_contact"
        )
        contact_threshold = upper_threshold if enable_threshold else None
        if selected_chains and selected_chains != st.session_state.get("_last_selected_chains"):
            st.session_state._last_selected_chains = selected_chains
            st.session_state.small_ref_dm = st.session_state.ref_dm.get_submatrix(selected_chains)
            
            st.session_state.small_tgt_dms_dict = {}
            for tgt, tgt_dm in st.session_state.tgt_dms_dict.items():
                st.session_state.small_tgt_dms_dict[tgt] = tgt_dm.get_submatrix(selected_chains)

            st.session_state.small_compare_dms_dict = {}
            for tgt, sub_compare_dm in st.session_state.compare_dms_dict.items():
                st.session_state.small_compare_dms_dict[tgt] = sub_compare_dm.get_submatrix(selected_chains)

        if selected_chains and "small_ref_dm" in st.session_state:
            high_res_selected_target = st.selectbox(
                "Select Target Structure to Display",
                target_names,
                key="selected_target_display_small",
            )

            full_res_slot_distance = st.empty()
            dist_scale = color_scale_ui(
                purpose="sequential",
                key_prefix="compare_dist_scale",
                default_colorscale="Viridis",
                min= float(0),
                max= upper_threshold if enable_threshold else max(
                    st.session_state.small_ref_dm.mat.max(),
                    float(max(v.mat.max() for v in st.session_state.small_tgt_dms_dict.values()))
                    )
            )
            full_res_slot_distance_diff = st.empty()
            diff_scale = color_scale_ui(
                purpose="diverging",
                key_prefix="compare_diff_scale",
                default_colorscale="RdBu_r",
                min = float(-max(abs(v.mat).max() for v in st.session_state.small_compare_dms_dict.values())),
                max = float(max(abs(v.mat).max() for v in st.session_state.small_compare_dms_dict.values()))
            )

            ref_fig, tgt_figs_dict, compare_figs_dict = plot_all_matrices_ensemble(
                st.session_state.small_ref_dm,
                st.session_state.small_tgt_dms_dict,
                st.session_state.small_compare_dms_dict,
                lower_threshold=lower_threshold,
                upper_threshold=upper_threshold,
                dist_colorscale=dist_scale["colorscale"],
                diff_colorscale=diff_scale["colorscale"],
                dist_vmin=dist_scale["vmin"],
                dist_vmax=dist_scale["vmax"],
                diff_vmax=diff_scale["vmax"],
            )
            st.session_state.ref_fig = ref_fig
            st.session_state.tgt_figs_dict = tgt_figs_dict
            st.session_state.compare_figs_dict = compare_figs_dict

            shared_fig = plot_comparison_with_contact_filter(
                                st.session_state.small_compare_dms_dict[selected_target],
                                contact_threshold=contact_threshold if contact_threshold else 0,
                                title=f"ΔDistance (shared contacts < {contact_threshold}Å)",
                                colorscale=diff_scale["colorscale"],
                                min_val=diff_scale["vmin"],
                                max_val=diff_scale["vmax"],
                            )

            with full_res_slot_distance.container():
                col_ref, col_tgt = st.columns(2)

                with col_ref:
                    st.plotly_chart(ref_fig, use_container_width=False)
                    add_svg_download(ref_fig, "full_ref_fig_map.svg")

                with col_tgt:
                    st.plotly_chart(
                        tgt_figs_dict[high_res_selected_target],
                        use_container_width=False,
                    )
                    add_svg_download(
                        tgt_figs_dict[high_res_selected_target],
                        "full_tgt_fig_map.svg",
                    )
            
            with full_res_slot_distance_diff.container():
                st.plotly_chart(
                    compare_figs_dict[high_res_selected_target] if not enable_threshold else shared_fig,
                    use_container_width=True,
                )
                add_svg_download(
                    compare_figs_dict[high_res_selected_target],
                    "compare_fig_map.svg",
                )

            if enable_threshold and contact_threshold is not None:
                small_compare_dm = st.session_state.small_compare_dms_dict[high_res_selected_target]
            
                gained_lost_slot = st.empty()
            

                gained_fig = plot_contacts_gained(
                    small_compare_dm,
                    contact_threshold=contact_threshold,
                )
                gained_fig.update_layout(width=600, height=600)

                lost_fig = plot_contacts_lost(
                    small_compare_dm,
                    contact_threshold=contact_threshold,
                )
                lost_fig.update_layout(width=600, height=600)

                with gained_lost_slot.container():
                    col1, col2 = st.columns(2)
                    with col1:
                        st.markdown("#### Contacts Gained in Target")
                        st.plotly_chart(gained_fig, use_container_width=True)
                        add_svg_download(gained_fig, "full_gained_fig_map.svg")

                    with col2:
                        st.markdown("#### Contacts Lost in Target")
                        st.plotly_chart(lost_fig, use_container_width=True)
                        add_svg_download(lost_fig, "full_lost_fig_map.svg")
