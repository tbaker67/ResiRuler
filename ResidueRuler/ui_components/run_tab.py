import streamlit as st
from ui_components.utils import save_temp_file,json_mapping_input
from src.resiruler.structure_parsing import load_structure, extract_residues_from_structure, convert_to_CA_coords_list
from src.resiruler.distance_calc import DistanceMatrix
from src.resiruler.plotting import plot_interactive_contact_map
from src.resiruler.chimera_export import generate_chimera_link_script
from ui_components.pymol_viewers import start_pymol_viewer, draw_links_pymol
import pandas as pd

def show_run_tab():
    cif_file = st.file_uploader("Upload Structure File (.cif)", type=["cif"])
    st.session_state.setdefault("run_clicked", False)
    st.session_state.setdefault("run_df", None)
    st.session_state.setdefault("run_contact_map", None)
    st.session_state.setdefault("run_matrix", None)
    
    default_data = pd.DataFrame({
        "Chain": ["A"], 
    })

    st.markdown("### Select Chains to Visualize")
    edited_df = st.data_editor(default_data, num_rows="dynamic", use_container_width=True)

    #Allow for selection based on input into table, and handle blank table 
    selected_chains = edited_df['Chain'].dropna().astype(str).str.strip()
    selected_chains = selected_chains[selected_chains != ""].tolist()
    if not selected_chains:
        selected_chains = None


    
    if st.button("Run"):
        if cif_file is not None:
            cif_path = save_temp_file(cif_file)
            structure = load_structure(cif_path)
            res_list = extract_residues_from_structure(structure, selected_chains)
            coords_list, index_map = convert_to_CA_coords_list(res_list)
            matrix = DistanceMatrix(coords_list, index_map)
            df = matrix.convert_to_df()
            contact_map = plot_interactive_contact_map(matrix, title="Contact Map")

            # save important data for plotting/tables
            st.session_state.run_clicked = True
            st.session_state.run_matrix = matrix
            st.session_state.run_df = df
            st.session_state.run_contact_map = contact_map
        else:
            st.error("No File Uploaded")


    # Show results if run was successful
    if st.session_state.run_clicked and st.session_state.run_df is not None:
        st.markdown("### Contact Map")
        st.plotly_chart(st.session_state.run_contact_map)

        st.markdown("### Distance Data")
        st.dataframe(st.session_state.run_df, use_container_width=True)

    st.markdown("### Edit or Paste Desired Residue Pairs")
    default_pairs = pd.DataFrame({
            'ChainID_Resnum1': ['A-58'],
            'ChainID_Resnum2': ['B-102']
        })
    


    link_selections = st.data_editor(default_pairs, num_rows="dynamic", use_container_width=True, key="residue_pairs_editor")
    thresholds_label = "Enter thresholding for yellow and green coloring by simply editing what is on the right of the ':'"
    thresholds_default = '''
    {
            "green": 23,
            "yellow": 33
    }
    '''
    
    key =  'run_thresholds'
    thresholds = json_mapping_input(thresholds_label, thresholds_default,key)


    st.session_state.setdefault("draw_links_clicked", False)
    st.session_state.setdefault("links_df", None)
    st.session_state.setdefault("link_script", None)
    st.session_state.setdefault("link_pymol", None)

    if st.button("Draw Links"):
        
        links_df = st.session_state.run_df.merge(link_selections, on=['ChainID_Resnum1', 'ChainID_Resnum2'])
        link_script = generate_chimera_link_script(links_df, thresholds=thresholds)

        st.session_state.draw_links_clicked = True
        st.session_state.links_df = links_df
        st.session_state.link_script = link_script
        st.session_state.link_pymol = draw_links_pymol(links_df, start_pymol_viewer(cif_file), thresholds=thresholds)

    if st.session_state.draw_links_clicked:

        st.markdown("### PyMol Visualization Preview")
        link_html = st.session_state.link_pymol._make_html()
        st.components.v1.html(link_html, height=600, width=1000)
        

        st.markdown("### Links Subset Data")
        st.dataframe(st.session_state.links_df)

        st.markdown("### Download ChimeraX Visualization Script")
        st.download_button("Download CXC Script",
                           data = st.session_state.link_script,
                           file_name = "link_script.cxc",
                           mime="text/plain")



