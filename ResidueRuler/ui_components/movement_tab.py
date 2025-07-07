import streamlit as st
import pandas as pd
import yaml
import json
from ui_components.pymol_viewers import draw_movement_shift_pymol, start_pymol_viewer, draw_movement_vectors_py3dmol
from ui_components.utils import save_temp_file, json_mapping_input
from src.resiruler import load_structure
from src.resiruler.distance_calc import calc_difference_aligned
from src.resiruler.plotting import plot_movement_shift, plot_movement_vectors
from src.resiruler.chimera_export import generate_chimera_color_script, generate_bild_string


def show_movement_tab():
    st.header("Movement Analysis Between Aligned Structures")

    
    cif1 = st.file_uploader("Upload Aligned CIF #1", type=["cif"], key="aligned1")
    cif2 = st.file_uploader("Upload Aligned CIF #2", type=["cif"], key="aligned2")
    #yaml_file = st.file_uploader("Upload Chain Mapping YAML (Optional)", type=["yaml", "yml"])

    
    st.session_state.setdefault("movement_df", None)
    st.session_state.setdefault("chimera_script", None)
    st.session_state.setdefault("bild_script", None)
    st.session_state.setdefault("movement_view", None)
    st.session_state.setdefault("vectors_view", None)

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

       
        
        
        #if yaml_file:
            #yaml_path = save_temp_file(yaml_file)
            #with open(yaml_path) as f:
                #chain_mapping = yaml.safe_load(f).get("chain_mapping")

        
        df = calc_difference_aligned(structure1, structure2, chain_mapping)
        st.session_state.movement_df = df

        
        st.session_state.movement_view =  draw_movement_shift_pymol(df, start_pymol_viewer(cif1_path))
        st.session_state.vector_view = draw_movement_vectors_py3dmol(df, start_pymol_viewer(cif1_path))
        st.session_state.chimera_script = generate_chimera_color_script(df)
        st.session_state.bild_script = generate_bild_string(df)

        st.success("Movement analysis complete!")

    
    if st.session_state.movement_df is not None:
        st.subheader("Movement Data")
        st.dataframe(st.session_state.movement_df)

        #Plotly Visualization
        #st.subheader("Δ Distance Plot")
        #st.plotly_chart(plot_movement_shift(st.session_state.movement_df, plotly=True))

        st.subheader("Δ Distance Visualization")
        html = st.session_state.movement_view._make_html()
        st.components.v1.html(html, height=600, width=1000)
        
        #Plotly Visualization
        #st.subheader("Movement Vectors")
        #st.plotly_chart(plot_movement_vectors(st.session_state.movement_df, plotly=True))

        st.subheader("Movement Vectors Pymol Visualization")
        vector_html = st.session_state.vector_view._make_html()
        st.components.v1.html(vector_html, height=600, width=1000)

        st.subheader("Download Options")

        csv_name = st.text_input("CSV Filename", value="residue_movement.csv")
        chimera_name = st.text_input("Chimera Coloring Script Filename", value="chimera_script.cxc")
        bild_name = st.text_input("Chimera Vector Script Filename", value="chimera_vectors.bild")
    
        st.download_button("Download CSV", data=st.session_state.movement_df.to_csv(index=False),
                           file_name=csv_name, mime="text/csv")

        st.download_button("Download Chimera Coloring Script",
                           data=st.session_state.chimera_script,
                           file_name=chimera_name, mime="text/plain")

        st.download_button("Download Chimera Vector (BILD) Script",
                           data=st.session_state.bild_script,
                           file_name=bild_name, mime="text/plain")
