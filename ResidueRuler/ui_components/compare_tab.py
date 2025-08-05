import streamlit as st
from ui_components.utils import json_mapping_input, create_mapper
from src.resiruler.plotting import plot_distance_difference, plot_interactive_contact_map
import numpy as np


def show_compare_tab():
    st.header("Compare Two CSV Outputs")

   
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

    if st.button("Compare"):
       
        try:
            mapper = create_mapper(cif1, cif2, chain_mapping)
            
        except ValueError as e:
            print(f"Error creating mapper: {e}")

        
        ref_dm, tgt_dm, compare_dm = mapper.calc_matrices()
        
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

    
        
        st.dataframe(st.session_state.compare_df)

        st.plotly_chart(st.session_state.ref_contact_map)

        st.plotly_chart(st.session_state.tgt_contact_map)

        st.plotly_chart(st.session_state.compare_contact_map)
        


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