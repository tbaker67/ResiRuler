import streamlit as st
from ui_components.utils import save_temp_file
from src.resiruler.structure_parsing import load_structure, extract_CA_coords
from src.resiruler.distance_calc import get_header_indices, read_data
from src.resiruler.chimera_export import generate_chimera_link_script
import yaml
import pandas as pd
from pathlib import Path

def show_run_tab():
    st.header("Run Distance Extraction")

    excel_file = st.file_uploader("Upload Excel File", type=["xlsx"])
    cif_file = st.file_uploader("Upload Structure File (.cif)", type=["cif"])
    yaml_file = st.file_uploader("Upload Config YAML", type=["yaml", "yml"])
    #visualize = st.checkbox("Generate ChimeraX script")

    st.session_state.setdefault("run_df", None)
    st.session_state.setdefault("run_script", None)
  

    if st.button("Run"):
        if not all([excel_file, cif_file, yaml_file]):
            st.error("Please upload all required files.")
            return

        excel_path = save_temp_file(excel_file)
        cif_path = save_temp_file(cif_file)
        yaml_path = save_temp_file(yaml_file)

        structure = load_structure(cif_path)
        index_map, coords = extract_CA_coords(structure)

        with open(yaml_path) as f:
            config = yaml.safe_load(f)

        df = pd.read_excel(excel_path).replace('', pd.NA).dropna(how='all').reset_index(drop=True)
        header_indices = get_header_indices(df)
        result_df = read_data(df, header_indices, config['chain_mapping'], index_map, coords)

        st.success("Distance table created!")
        st.dataframe(result_df)
        st.session_state.run_df = result_df

        st.session_state.run_script = generate_chimera_link_script(result_df)

        if st.session_state.run_df is not None:

            st.subheader("Distance Table")
            st.dataframe(st.session_state.run_df)

            st.subheader("Download Options")

            csv_name = st.text_input("CSV Filename", value="residue_movement.csv")
            chimera_name = st.text_input("Chimera Coloring Script Filename", value="chimera_script.cxc")

            st.download_button("Download CSV", data=st.session_state.movement_df.to_csv(index=False),
                           file_name=csv_name, mime="text/csv")

            st.download_button("Download Chimera Link Script",
                           data=st.session_state.chimera_script,
                           file_name=chimera_name, mime="text/plain")

        