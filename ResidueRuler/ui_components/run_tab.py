import streamlit as st
from ui_components.utils import save_temp_file,json_mapping_input
from src.resiruler.structure_parsing import load_structure, extract_CA_coords
from src.resiruler.distance_calc import get_header_indices, read_data
from src.resiruler.chimera_export import generate_chimera_link_script
from ui_components.pymol_viewers import start_pymol_viewer, draw_links_pymol
import yaml
import pandas as pd
from pathlib import Path

def show_run_tab():
    st.header("Run Distance Extraction")

    excel_file = st.file_uploader("Upload Excel File", type=["xlsx"])
    cif_file = st.file_uploader("Upload Structure File (.cif)", type=["cif"])
    #yaml_file = st.file_uploader("Upload Config YAML", type=["yaml", "yml"])
    #visualize = st.checkbox("Generate ChimeraX script")

    st.session_state.setdefault("run_df", None)
    st.session_state.setdefault("run_script", None)

    chain_mapping = None

    label = "Enter a chain mapping in the following format such that to the left of the ':' is the chain label specific in the input excel file, and to the right is the chain id which shows up in chimera/the cif. If you have multiple chains to match up to an identifer, simply separate them inside the '[]' with commas, but always make sure to wrap the right side in '[]'" 

    default = '''
    { 
            "Protein A":["AA","DA"],
            "Protein B":["AB","CB"],
            "Protein C":["CC"],
            "Protein D":["DD"]
    }
    '''

    key = 'run_mapping'

    chain_mapping = json_mapping_input(label,default,key)
  
    thresholds_label = "Enter thresholding for yellow and green coloring by simply editing what is on the right of the ':'"

    thresholds_default = '''
    {
            "green": 23,
            "yellow": 33
    }
    '''

    key =  'run_thresholds'
    thresholds = json_mapping_input(thresholds_label, thresholds_default,key)
    if st.button("Run"):
        if not all([excel_file, cif_file]):
            st.error("Please upload all required files.")
            return

        excel_path = save_temp_file(excel_file)
        cif_path = save_temp_file(cif_file)
        

        structure = load_structure(cif_path)
        index_map, coords = extract_CA_coords(structure)


        df = pd.read_excel(excel_path).replace('', pd.NA).dropna(how='all').reset_index(drop=True)
        header_indices = get_header_indices(df)

        result_df = read_data(df, header_indices, chain_mapping, index_map, coords)

        st.success("Distance table created!")
        st.session_state.run_df = result_df

        view = draw_links_pymol(result_df, start_pymol_viewer(cif_path), thresholds)
        


        st.session_state.run_script = generate_chimera_link_script(result_df, chain_mapping, thresholds)



        if st.session_state.run_df is not None:

            st.subheader("Distance Table")
            st.dataframe(st.session_state.run_df)

            st.subheader("PyMol Visualization")
            html = view._make_html()
            st.components.v1.html(html, height=600, width=1000)

            st.subheader("Download Options")

            csv_name = st.text_input("CSV Filename", value="residue_movement.csv")
            chimera_name = st.text_input("Chimera Link Script Filename", value="chimera_script.cxc")

            st.download_button("Download CSV", data=st.session_state.run_df.to_csv(index=False),
                           file_name=csv_name, mime="text/csv")

            st.download_button("Download Chimera Link Script",
                           data=st.session_state.run_script,
                           file_name=chimera_name, mime="text/plain")

        