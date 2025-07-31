import streamlit as st
import pandas as pd
import yaml
import json
from ui_components.pymol_viewers import draw_movement_shift_pymol, start_pymol_viewer, draw_movement_vectors_py3dmol
from ui_components.utils import save_temp_file, json_mapping_input, create_downloadable_zip
from src.resiruler import load_structure
from src.resiruler.distance_calc import calc_difference_aligned
from src.resiruler.plotting import plot_colorbar
from src.resiruler.chimera_export import generate_cxc_scripts, generate_bild_string
import os
from pathlib import Path

def show_movement_tab():
    st.header("Movement Analysis Between Aligned Structures")

    
    cif1 = st.file_uploader("Upload Aligned CIF #1", type=["cif"], key="aligned1")
    cif2 = st.file_uploader("Upload Aligned CIF #2", type=["cif"], key="aligned2")
    #yaml_file = st.file_uploader("Upload Chain Mapping YAML (Optional)", type=["yaml", "yml"])

    
    st.session_state.setdefault("movement_df", None)
    st.session_state.setdefault("chimera_script", None)
    st.session_state.setdefault("bild_script", None)
    st.session_state.setdefault("defatt1", None)
    st.session_state.setdefault("defatt2", None)
    st.session_state.setdefault("movement_view", None)
    st.session_state.setdefault("vector_view", None)
    st.session_state.setdefault("structure1_name", None)
    st.session_state.setdefault("structure2_name", None)

    chain_mapping = None

    label = "Enter a chain mapping in the following format such that to the left of the ':' is a chain id in structure 1, and to the right is the corresponding chain id in structure 2" 

    default = '''{ 
            "AA":"ZZ",
            "BB":"YY",
            "CC":"XX",
            "DD":"WW"
    }'''

    key = 'movement'
    chain_mapping = json_mapping_input(label,default,key)
    
    if st.button("Analyze Movement"):
        if not cif1 or not cif2:
            st.error("Please upload both structure files.")
            return

        
        cif1_path = save_temp_file(cif1)
        cif2_path = save_temp_file(cif2)

       
        structure1 = load_structure(cif1_path)
        structure2 = load_structure(cif2_path)

       
        
        df = calc_difference_aligned(structure1, structure2, chain_mapping)
        st.session_state.movement_df = df

        st.session_state.structure1_name = os.path.splitext(cif1.name)[0]
        st.session_state.structure2_name = os.path.splitext(cif2.name)[0]
        st.session_state.movement_view =  draw_movement_shift_pymol(df, start_pymol_viewer(cif1_path))
        st.session_state.vector_view = draw_movement_vectors_py3dmol(df, start_pymol_viewer(cif1_path))
        defatt1, defatt2, chimera_cxc = generate_cxc_scripts(df, cif1.name, cif2.name,
                                                                           cif1_path, cif2_path, 
                                                                           st.session_state.structure1_name, 
                                                                           st.session_state.structure2_name, 
                                                                           chain_mapping)
        st.session_state.defatt1 = defatt1
        st.session_state.defatt2 = defatt2
        st.session_state.chimera_script = chimera_cxc
        st.session_state.bild_script = generate_bild_string(df)

        st.success("Movement analysis complete!")

    
    if st.session_state.movement_df is not None:
        st.subheader("Movement Data")
        st.dataframe(st.session_state.movement_df)

        st.subheader("Color Legend (Distance Shift)")
        min_shift = st.session_state.movement_df['Distance'].min()
        max_shift = st.session_state.movement_df['Distance'].max()
        fig = plot_colorbar(min_shift, max_shift, cmap_name="plasma")  
        st.pyplot(fig)

        st.subheader("Distance Shift Visualization")
        html = st.session_state.movement_view._make_html()
        st.components.v1.html(html, height=600, width=1000)
        


        st.subheader("Movement Vectors Pymol Visualization")
        vector_html = st.session_state.vector_view._make_html()
        st.components.v1.html(vector_html, height=600, width=1000)


        defatt1_name = f"{st.session_state.structure1_name}_colors.defattr"
        defatt2_name = f"{st.session_state.structure2_name}_colors.defattr"
        chimera_filename = "chimera_coloring_script.cxc"
        bild_filename = "colored_vectors.bild"
        csv_filename = "residue_movement.csv"

        # read back uploaded CIF content
        cif1_filename = Path(cif1.name).name
        cif2_filename = Path(cif2.name).name
        cif1_content = cif1.getvalue().decode("utf-8")
        cif2_content = cif2.getvalue().decode("utf-8")

        # create downloadable ZIP
        files_to_zip = {
            chimera_filename: st.session_state.chimera_script,
            defatt1_name:  st.session_state.defatt1,
            defatt2_name:  st.session_state.defatt2,
            cif1_filename: cif1_content,
            cif2_filename: cif2_content,
            csv_filename: st.session_state.movement_df.to_csv(index=False),
            bild_filename: st.session_state.bild_script 
        }
        zip_buffer = create_downloadable_zip(files_to_zip)

        st.session_state.zip_buffer = zip_buffer

       

        st.subheader("Download Full Data and Scripts Folder")

        st.download_button("Download All as ZIP",
                           data=st.session_state.zip_buffer,
                           file_name="movement_analysis_package.zip",
                           mime="application/zip")