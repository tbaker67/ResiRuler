import streamlit as st
from ui_components.pymol_viewers import draw_movement_shift_pymol, start_pymol_viewer, draw_movement_vectors_py3dmol
from ui_components.utils import json_mapping_input, create_downloadable_zip, create_ensemble_mapper, load_structure_if_new, get_threshold, load_structures_if_new, get_chain_mappings_for_targets, create_downloadable_zip_grouped
from src.resiruler.distance_calc import calc_difference_from_mapper
from src.resiruler.plotting import plot_colorbar
from src.resiruler.chimera_export import generate_multiple_bilds, generate_multiple_movement_scripts
import os
from pathlib import Path
from Bio.Align import PairwiseAligner, substitution_matrices

def show_movement_tab():
    st.header("Movement Analysis Between Aligned Structures")

    st.session_state.setdefault("movement_dfs", None)
    
    ref_cif = st.file_uploader("Upload Aligned Reference CIF", type=["cif"], key="movement_reference")
    tgt_cifs = st.file_uploader("Upload Aligned Target CIFs", type=["cif"], key="movement_tgts", accept_multiple_files=True)

    
    ref_structure = load_structure_if_new(ref_cif, "compare_name1", "compare_structure1")
    tgt_structures = load_structures_if_new(tgt_cifs, "compare_name2", "compare_structure2")

   
    if ref_structure and tgt_structures:
        default_mapping = '''{
            "AA": "ZZ",
            "BB": "YY",
            "CC": "XX",
            "DD": "WW"
        }'''
        chain_mappings = get_chain_mappings_for_targets(tgt_structures, default_mapping)


    #get threshold and do alignments
    pct_id_threshold = get_threshold("Set a minimum pct Identity Threshold for matching chains together", "95.0")
    
    st.session_state.setdefault("mapper", None)

    ##TODO: Make functino to allow user to set aligner with parameters like gap penalties etc
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.left_open_gap_score = 1
    aligner.open_gap_score = -10

    if st.button("Map Chains"):
        st.session_state.mapper = create_ensemble_mapper(ref_structure, tgt_structures, chain_mappings, pct_id_threshold, aligner)
        st.session_state.mapper.set_selected_global_coords()


    if st.button("Analyze Movement"):

        st.session_state.movement_dfs = st.session_state.mapper.calc_movement_dfs()

        st.session_state.ref_name = os.path.splitext(ref_cif.name)[0]
    
        
        defatt, chimera_cxc = generate_multiple_movement_scripts(st.session_state.movement_dfs, st.session_state.ref_name)
        st.session_state.defatt = defatt
        st.session_state.chimera_script = chimera_cxc
        st.session_state.bild_scripts = generate_multiple_bilds(st.session_state.movement_dfs)


        st.success("Movement analysis complete!")

    
    if st.session_state.movement_dfs is not None:

        root_files = {
        f"{st.session_state.ref_name}_colors.defattr": st.session_state.defatt,
        "chimera_coloring_script.cxc": st.session_state.chimera_script
    }

        # --- CSV files per structure ---
        csv_dict = {}
        for struct_name, df in st.session_state.movement_dfs.items():
            csv_filename = f"{struct_name}_movement.csv"
            csv_dict[csv_filename] = df.to_csv(index=False)

        # --- BILD files ---
        bild_dict = st.session_state.bild_scripts  # assuming {filename: content}

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

        st.subheader("Color Legend (Distance Shift)")
        #min_shift = st.session_state.movement_df['Distance'].min()
        #max_shift = st.session_state.movement_df['Distance'].max()
        #fig = plot_colorbar(min_shift, max_shift, cmap_name="plasma")  
        #st.pyplot(fig)


        structure_choices = {
        f.name[:-4]:f for f in tgt_cifs
        }
        

        selected_structure = st.selectbox(
        "Select structure for visualization",
        options=list(structure_choices.keys())
        )

        # Initialize viewer from chosen file
        viewer = start_pymol_viewer(ref_cif)

        # Draw both movement visualizations
        st.session_state.movement_view = draw_movement_shift_pymol(
        st.session_state.movement_dfs[selected_structure], viewer

        )
        st.session_state.vector_view = draw_movement_vectors_py3dmol(
        st.session_state.movement_dfs[selected_structure], viewer
        )
        st.subheader("Distance Shift Visualization")
        html = st.session_state.movement_view._make_html()
        st.components.v1.html(html, height=600, width=1000)
        

        st.subheader("Movement Vectors Pymol Visualization")
        vector_html = st.session_state.vector_view._make_html()
        st.components.v1.html(vector_html, height=600, width=1000)

        st.subheader("Movement Data")
        st.dataframe(st.session_state.movement_dfs[selected_structure])


