import streamlit as st
from ui_components.pymol_viewers import draw_movement_shift_pymol, start_pymol_viewer, draw_movement_vectors_py3dmol, plot_vectors_plotly
from ui_components.utils import create_downloadable_zip, create_ensemble_mapper, load_structure_if_new, get_threshold, load_structures_if_new, get_chain_mappings_for_targets, create_downloadable_zip_grouped, aligner_ui, show_alignments, get_measurement_mode,  struct_to_temp_cif
from ui_components.color_mapping_utils import gradient_palette_picker, build_gradient_cmap, show_gradient_bar, get_coloring_values
from ui_components.molstar_viewers import create_distance_shift_builder, write_movement_annotations
from src.resiruler.distance_calc import calc_difference_from_mapper
from src.resiruler.plotting import plot_colorbar
from src.resiruler.chimera_export import generate_arrow_dicts, generate_multiple_movement_scripts, generate_pml_arrows
import os
from pathlib import Path
from Bio.Align import PairwiseAligner, substitution_matrices
import json

def show_movement_tab():
    st.header("Movement Analysis Between Aligned Structures")

    st.session_state.setdefault("movement_dfs", None)
    
    ref_cif = st.file_uploader("Upload Aligned Reference CIF", type=["cif"], key="movement_reference")
    tgt_cifs = st.file_uploader("Upload Aligned Target CIFs", type=["cif"], key="movement_tgts", accept_multiple_files=True)

    
    ref_structure = load_structure_if_new(ref_cif, "compare_name1", "compare_structure1")
    tgt_structures = load_structures_if_new(tgt_cifs, "compare_name2", "compare_structure2")


    if ref_structure and tgt_structures:
        ref_chains = [ref_chain.id for ref_chain in ref_structure[0].get_chains()]
        chain_mappings = get_chain_mappings_for_targets(tgt_structures,ref_chains, key="movement_mappings")


    #get threshold and do alignments
    pct_id_threshold = get_threshold("Set a minimum pct Identity Threshold for matching chains together", "95.0")
    
    st.session_state.setdefault("mapper", None)

    aligner = aligner_ui("movement aligner")

    if st.button("Map Chains"):
        st.session_state.mapper = create_ensemble_mapper(ref_structure, tgt_structures, chain_mappings, pct_id_threshold, aligner)
        

        
        
    if st.session_state.mapper is not None:
        show_alignments(st.session_state.mapper, key="movement_alignment")

    mode = get_measurement_mode(key="movement_measurement_mode")
    
    if st.button("Analyze Movement"):

        st.session_state.mapper.set_selected_global_coords(mode=mode)

        st.session_state.movement_dfs = st.session_state.mapper.calc_movement_dfs()



        st.success("Movement analysis complete!")

    
    if st.session_state.movement_dfs is not None:

        
        mins = [df['Distance'].min() for df in st.session_state.movement_dfs.values()]
        vmin = min(mins)
        
        maxes = [df['Distance'].max() for df in st.session_state.movement_dfs.values()]
        vmax = max(maxes)

        default_colors = ["#00008B","#20073a", "#6d1950", "#bd4545", "#d48849", "#f0d171"]
        # Show gradient color picker and preview

        palette = gradient_palette_picker(default_colors=default_colors, key="movement_palette_picker")
        min_val,max_val = get_coloring_values(vmin, vmax)

        show_gradient_bar(palette, min_val=min_val, max_val=max_val)
        cmap_obj = build_gradient_cmap(palette, max_val, min_val)

        structure_choices = {
        f.name[:-4]:f for f in tgt_cifs
        }
        

        selected_structure = st.selectbox(
        "Select structure for visualization",
        options=list(structure_choices.keys())
        )
        
        
        if "vector_view" in st.session_state:
            del st.session_state.vector_view

        
        viewer1 = start_pymol_viewer(ref_cif)
        st.session_state.vector_view = plot_vectors_plotly(st.session_state.movement_dfs[selected_structure], cmap_obj,min_val,max_val)
        st.subheader("Movement Vectors Preview Visualization")
        st.plotly_chart(st.session_state.vector_view, use_container_width=True)

        # create a molstar view builder
        if "molstar_builder" not in st.session_state:
            st.session_state.molstar_builder = create_distance_shift_builder()
            
        builder = st.session_state.molstar_builder

        # create temp files for the builder to read in
        ref_cif_path = struct_to_temp_cif(ref_structure)
        tgt_cif_path = struct_to_temp_cif(tgt_structures[selected_structure])

        annotations = write_movement_annotations(st.session_state.movement_dfs[selected_structure], cmap_obj, min_val, max_val )
        annotations_json = json.dumps(annotations)

        with struct_to_temp_cif(ref_structure) as ref_cif_path, \
             struct_to_temp_cif(tgt_structures[selected_structure]) as tgt_cif_path:

            ref_cif_data = open(ref_cif_path).read()
            tgt_cif_data = open(tgt_cif_path).read()
            builder.molstar_streamlit(
                data={"local.cif": ref_cif_data, "annotations.json": annotations_json},
                width=1600,
                height=600,
            )

            annotations2 = write_movement_annotations(
                st.session_state.movement_dfs[selected_structure], cmap_obj, min_val, max_val, ref=False, 
            )
            annotations_json2 = json.dumps(annotations2)
            builder.molstar_streamlit(
                data={"local.cif": tgt_cif_data, "annotations.json": annotations_json2},
            )

        
        st.subheader("Movement Data")
        st.dataframe(st.session_state.movement_dfs[selected_structure])


        st.session_state.ref_name = os.path.splitext(ref_cif.name)[0]
        defatt, chimera_cxc, full_pml_script = generate_multiple_movement_scripts(st.session_state.movement_dfs, st.session_state.ref_name, palette, min_val, max_val)
        
        st.session_state.defatt = defatt
        st.session_state.chimera_script = chimera_cxc
        st.session_state.pml_script = full_pml_script
        st.session_state.bild_scripts, st.session_state.pml_arrows = generate_arrow_dicts(st.session_state.movement_dfs, cmap_obj, vmin, vmax)

        root_files = {
        f"full_defattr.defattr": st.session_state.defatt,
        "chimera_coloring_script.cxc": st.session_state.chimera_script,
        "coloring_script.pml":st.session_state.pml_script
    }

        # --- CSV files per structure ---
        csv_dict = {}
        for struct_name, df in st.session_state.movement_dfs.items():
            csv_filename = f"{struct_name}_movement.csv"
            csv_dict[csv_filename] = df.to_csv(index=False)

        # --- Arrow files ---
        bild_dict = st.session_state.bild_scripts  # assuming {filename: content}
        pml_arrows_dict = st.session_state.pml_arrows
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
            None: root_files,   # root folder
            "csv": csv_dict,
            "bild": bild_dict,
            "pml_arrows":pml_arrows_dict,
            "models": models_dict
        }

        zip_buffer = create_downloadable_zip_grouped(file_groups)

        st.subheader("Download Full Data and Scripts Folder")
        st.download_button(
            "Download All as ZIP",
            data=zip_buffer,
            file_name="resi_ruler_movement_output.zip",
            mime="application/zip"
        )


