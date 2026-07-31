"""displacement analysis tab for comparing structural shifts between aligned structures."""

import json
import os
from pathlib import Path

import streamlit as st
from src.resiruler.viz.export_visualizations import (
    generate_arrow_dicts,
    generate_multiple_displacement_scripts,
    save_bild_files_and_generate_chimerax_script,
    save_pml_files_and_generate_pymol_script,
)

from ui.viewers.molstar_viewers import (
    create_distance_shift_builder,
    write_displacement_annotations,
)
from ui.viewers.plotly_viewer import plot_vectors_plotly
from ui.widgets.color_mapping_utils import (
    build_gradient_cmap,
    gradient_palette_picker,
    show_gradient_bar,
)
from ui.widgets.utils import (
    chain_selector_ui,
    create_downloadable_zip_grouped,
    create_ensemble_mapper,
    filter_df_by_chains,
    full_aligner_ui,
    get_chain_mappings_for_targets,
    get_measurement_mode,
    get_threshold,
    load_structure_if_new,
    load_structures_if_new,
    reset_downstream,
    show_alignments,
    struct_to_temp_cif,
)


def show_displacement_tab():
    st.header("Displacement Analysis Between Aligned Structures")

    st.session_state.setdefault("displacement_dfs", None)

    ref_cif = st.file_uploader(
        "Upload Aligned Reference CIF",
        type=["cif"],
        key="displacement_reference",
    )
    tgt_cifs = st.file_uploader(
        "Upload Aligned Target CIFs",
        type=["cif"],
        key="displacement_tgts",
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
            tgt_structures, ref_chains, key="displacement_mappings"
        )

    st.session_state.setdefault("disp_mapper", None)
    protein_aligner, nucleotide_aligner, rmsd_scaling_factor = full_aligner_ui(key="displacement")

    pct_id_threshold = get_threshold(
        "Set a Minimum Percent Identity Threshold for Matching Chains Together",
        "95.0",
        "displacement_pct_id",
    )

    if st.button("Map Chains", key="map displacement chains", help="Find chain-chain and residue-residue correspondences"):
        reset_downstream(
            "disp_mapper",
            "displacement_dfs",
            "vector_view",
            "molstar_builder",
            "defatt",
            "chimera_script",
            "pml_script",
            "bild_scripts",
            "pml_arrows",
            "zip_buffer",
        )
        st.session_state.disp_mapper = create_ensemble_mapper(
            ref_structure,
            tgt_structures,
            chain_mappings,
            pct_id_threshold,
            protein_aligner,
            nucleotide_aligner,
            rmsd_scaling_factor
        )

    if st.session_state.get("disp_mapper") is not None:
        show_alignments(st.session_state.disp_mapper, key="displacement_alignment")

    protein_mode, nucleic_mode = get_measurement_mode(
        key="displacement_measurement_mode"
    )

    if st.button("Analyze Displacement"):
        reset_downstream(
            "displacement_dfs",
            "vector_view",
            "molstar_builder",
            "defatt",
            "chimera_script",
            "pml_script",
            "bild_scripts",
            "pml_arrows",
            "zip_buffer",
        )
        st.session_state.disp_mapper.set_selected_global_coords(
            protein_mode=protein_mode, nucleic_mode=nucleic_mode
        )
        st.session_state.displacement_dfs = (
            st.session_state.disp_mapper.calc_displacement_dfs()
        )
        st.success("Displacement Analysis Complete")

    if st.session_state.get("displacement_dfs") is not None:
        selected_chains = chain_selector_ui(
            ref_structure,
            "Select Chains in Reference to Compare",
            key_prefix="",
        )
        filtered_displacement_dfs = {
            struct_name: filter_df_by_chains(df, selected_chains)
            for struct_name, df in st.session_state.displacement_dfs.items()
        }
        mins = [df["Distance"].min() for df in filtered_displacement_dfs.values()]
        min(mins)

        maxes = [df["Distance"].max() for df in filtered_displacement_dfs.values()]
        vmax = max(maxes)

        default_colors = [
            "#00008B",
            "#20073a",
            "#6d1950",
            "#bd4545",
            "#d48849",
            "#f0d171",
        ]
        # Show gradient color picker and preview
        palette, positions = gradient_palette_picker(
            0,
            vmax + 1,
            default_colors=default_colors,
            key="displacement_palette_picker",
        )
        min_val = min(positions)
        max_val = max(positions)
        show_gradient_bar(palette, positions, min_val=min_val, max_val=max_val)
        cmap_obj = build_gradient_cmap(palette, positions, min_val, max_val)

        structure_choices = {f.name[:-4]: f for f in tgt_cifs}

        selected_structure = st.selectbox(
            "Select Structure For Visualization",
            options=list(structure_choices.keys()),
        )

        fidelity = st.slider(
            "Vector Fidelity (show every Nth vector)",
            min_value=1,
            max_value=20,
            value=5,
            step=1,
            key="displacement_fidelity",
        )

        st.session_state.vector_view = plot_vectors_plotly(
            filtered_displacement_dfs[selected_structure],
            cmap_obj,
            min_val,
            max_val,
            fidelity=fidelity,
        )
        st.subheader("Displacement Vectors Preview Visualization")
        st.plotly_chart(st.session_state.vector_view, use_container_width=True)

        # create a molstar view builder
        if "molstar_builder" not in st.session_state:
            st.session_state.molstar_builder = create_distance_shift_builder()

        builder = st.session_state.molstar_builder

        # create temp files for the builder to read in
        ref_cif_path = struct_to_temp_cif(ref_structure)
        tgt_cif_path = struct_to_temp_cif(tgt_structures[selected_structure])

        annotations = write_displacement_annotations(
            filtered_displacement_dfs[selected_structure],
            cmap_obj,
            min_val,
            max_val,
        )
        annotations_json = json.dumps(annotations)

        with (
            struct_to_temp_cif(ref_structure) as ref_cif_path,
            struct_to_temp_cif(tgt_structures[selected_structure]) as tgt_cif_path,
        ):
            ref_cif_data = open(ref_cif_path).read()
            tgt_cif_data = open(tgt_cif_path).read()
            builder.molstar_streamlit(
                data={
                    "local.cif": ref_cif_data,
                    "annotations.json": annotations_json,
                }
            )

            annotations2 = write_displacement_annotations(
                filtered_displacement_dfs[selected_structure],
                cmap_obj,
                min_val,
                max_val,
                ref=False,
            )
            annotations_json2 = json.dumps(annotations2)
            builder.molstar_streamlit(
                data={
                    "local.cif": tgt_cif_data,
                    "annotations.json": annotations_json2,
                },
            )

        st.subheader("Displacement Data")
        st.dataframe(filtered_displacement_dfs[selected_structure])

        st.subheader("Download Full Data and Scripts Folder")

        if st.button("Generate Download Package", key="generate_displacement_zip"):
            st.session_state.ref_name = os.path.splitext(ref_cif.name)[0]
            defatt, chimera_cxc, full_pml_script = generate_multiple_displacement_scripts(
                filtered_displacement_dfs,
                st.session_state.ref_name,
                palette,
                positions,
            )

            st.session_state.defatt = defatt
            st.session_state.chimera_script = chimera_cxc
            st.session_state.pml_script = full_pml_script
            st.session_state.bild_scripts, st.session_state.pml_arrows = (
                generate_arrow_dicts(
                    filtered_displacement_dfs,
                    cmap_obj,
                    min_val,
                    max_val,
                    fidelity=fidelity,
                )
            )

            root_files = {
                "full_defattr.defattr": st.session_state.defatt,
                "chimera_coloring_script.cxc": st.session_state.chimera_script,
                "coloring_script.pml": st.session_state.pml_script,
            }

            # --- CSV files per structure ---
            csv_dict = {}
            for struct_name, df in filtered_displacement_dfs.items():
                csv_filename = f"{struct_name}_displacement.csv"
                csv_dict[csv_filename] = df.to_csv(index=False)

            # --- Arrow files ---
            bild_dict = st.session_state.bild_scripts  # {filename: content}
            pml_arrows_dict = st.session_state.pml_arrows

            chimerax_script_content = save_bild_files_and_generate_chimerax_script(
                bild_dict
            )
            bild_dict["open_all_bilds.cxc"] = chimerax_script_content

            pymol_arrows_script_content = save_pml_files_and_generate_pymol_script(
                pml_arrows_dict
            )
            pml_arrows_dict["run_all_arrows.pml"] = pymol_arrows_script_content

            # --- models (CIFs) ---
            models_dict = {}
            if ref_cif:
                ref_name = Path(ref_cif.name).name
                models_dict[ref_name] = ref_cif.getvalue().decode("utf-8")
            for cif in tgt_cifs:
                cif_name = Path(cif.name).name
                models_dict[cif_name] = cif.getvalue().decode("utf-8")

            # --- combine into grouped ZIP ---
            file_groups = {
                None: root_files,  # root folder
                "csv": csv_dict,
                "bild": bild_dict,
                "pml_arrows": pml_arrows_dict,
                "models": models_dict,
            }

            st.session_state.zip_buffer = create_downloadable_zip_grouped(file_groups)

        if st.session_state.get("zip_buffer") is not None:
            st.download_button(
                "Download All as ZIP",
                data=st.session_state.zip_buffer,
                file_name="resi_ruler_displacement_output.zip",
                mime="application/zip",
            )