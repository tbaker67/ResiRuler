import streamlit as st
import yaml
from ui_components.utils import save_temp_file, json_mapping_input
from src.resiruler.plotting import plot_distance_difference
from src.resiruler.chimera_export import generate_chimera_color_script

def show_compare_tab():
    st.header("Compare Two CSV Outputs")

   
    csv1 = st.file_uploader("Upload First CSV", type=["csv"], key="compare_csv1")
    csv2 = st.file_uploader("Upload Second CSV", type=["csv"], key="compare_csv2")
    yaml_file = st.file_uploader("Upload Config YAML (Optional)", type=["yaml", "yml"])

    chain_mapping = None

    label = "Enter a chain mapping in the following format such that to the left of the ':' is a chain id in structure 1, and to the right is the corresponding chain id in structure 2" 

    default = '''{ 
            "AA":"ZZ",
            "BB":"YY",
            "CC":"XX",
            "DD":"WW"
    }'''

    chain_mapping = json_mapping_input(label,default)

    if st.button("Compare"):
        if not csv1 or not csv2:
            st.error("Please upload both CSV files.")
            return

       
        csv1_path = save_temp_file(csv1)
        csv2_path = save_temp_file(csv2)

        chain_mapping = None
        if yaml_file:
            yaml_path = save_temp_file(yaml_file)
            with open(yaml_path) as f:
                chain_mapping = yaml.safe_load(f).get("chain_mapping")

        # Perform comparison
        fig, merged = plot_distance_difference(csv1_path, csv2_path, chain_map=chain_mapping, plotly=True)

        if merged is None or merged.empty:
            st.warning("No valid pairs found to compare.")
            return

        # Store results in session
        st.session_state.compare_merged = merged
        st.session_state.compare_fig = fig

        st.success("Comparison complete!")

    
    if "compare_merged" in st.session_state:
        st.subheader("Distance Differences")
        st.plotly_chart(st.session_state.compare_fig)
        st.dataframe(st.session_state.compare_merged)

        
        st.subheader("Download Options")
        csv_name = st.text_input("CSV Output Filename", value="distance_difference.csv")

       
        st.download_button("Download CSV",
                           data=st.session_state.compare_merged.to_csv(index=False),
                           file_name=csv_name,
                           mime="text/csv")