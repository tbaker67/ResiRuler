import streamlit as st
from ui_components.utils import json_mapping_input, create_mapper, chain_selector_ui, load_structure_if_new, get_threshold
from src.resiruler.plotting import plot_distance_difference, plot_interactive_contact_map
import numpy as np


def show_compare_tab():
    st.header("Compare Distances within two structures")

   
    cif1 = st.file_uploader("Upload Aligned CIF #1", type=["cif"], key="aligned1")
    cif2 = st.file_uploader("Upload Aligned CIF #2", type=["cif"], key="aligned2")
 

    chain_mapping = None

    label = "Enter a chain mapping in the following format such that to the left of the ':' is a chain id in structure 1, and to the right is the corresponding chain id in structure 2" 

    default = '''{ 
            "AA":"ZZ",
            "BB":"YY",
            "CC":"XX",
            "DD":"WW"
    }'''

    key = 'compare'

    chain_mapping = json_mapping_input(label,default,key)

    structure1 = load_structure_if_new(cif1, name_key="compare_name1", struct_key="compare_structure1")
    structure2 = load_structure_if_new(cif2, name_key="compare_name2", struct_key="compare_structure2")

    selected_chains = chain_selector_ui(structure1, "Select Chains in reference to compare")

    pct_id_threshold = get_threshold("Set a minimum Pct Identity Threshold for matching chains together", "95.0")

    if st.button("Compare"):
       
        try:
            mapper = create_mapper(structure1, structure2, chain_mapping, pct_id_threshold)
            
            
        except ValueError as e:
            print(f"Error creating mapper: {e}")

        
        ref_dm, tgt_dm, compare_dm = mapper.calc_matrices(selected_chains)
        
        min_val = min(np.nanmin(ref_dm.mat), np.nanmin(tgt_dm.mat))
        max_val = max(np.nanmax(ref_dm.mat), np.nanmax(tgt_dm.mat))

        st.session_state.ref_contact_map = plot_interactive_contact_map(ref_dm, title="Reference Contact Map", min=min_val, max=max_val)

        st.session_state.tgt_contact_map = plot_interactive_contact_map(tgt_dm, title = "Target Contact Map", min=min_val,max=max_val)

        st.session_state.compare_contact_map = plot_interactive_contact_map(compare_dm, title = "Distance Difference Contact Map (Target - Reference)")

        st.session_state.compare_df = compare_dm.convert_to_df()

        #if merged is None or merged.empty:
           # st.warning("No valid pairs found to compare.")
            #return

        # Store results in session
       # st.session_state.compare_merged = merged
        #st.session_state.compare_fig = fig

        st.success("Comparison complete!")

    
        
        st.dataframe(st.session_state.compare_df.drop(columns=['Coord1_ref', 'Coord2_ref',
                                                               'Coord1_tgt', 'Coord2_tgt']), use_container_width=True)

        st.plotly_chart(st.session_state.ref_contact_map, use_container_width=False)

        st.plotly_chart(st.session_state.tgt_contact_map,use_container_width=False)

        st.plotly_chart(st.session_state.compare_contact_map, use_container_width=False)
        


    #if "compare_merged" in st.session_state:
       # st.subheader("Distance Differences")
       # st.plotly_chart(st.session_state.compare_fig)
       # st.dataframe(st.session_state.compare_merged)

        
      #  st.subheader("Download Options")
       # csv_name = st.text_input("CSV Output Filename", value="distance_difference.csv")

       
       # st.download_button("Download CSV",
       #                    data=st.session_state.compare_merged.to_csv(index=False),
       #                    file_name=csv_name,
       #                    mime="text/csv")