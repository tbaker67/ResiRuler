import streamlit as st
from ui_components.utils import save_temp_file
from src.resiruler.plotting import plot_distance_difference
import os
import yaml

def show_compare_tab():
    st.header("Compare Two CSV Outputs")

    csv1 = st.file_uploader("Upload First CSV", key="compare_csv1")
    csv2 = st.file_uploader("Upload Second CSV", key="compare_csv2")
    yaml_file = st.file_uploader("Upload Config YAML", type=["yaml", "yml"])

    if st.button("Compare"):
        if not csv1 or not csv2:
            st.error("Please upload both CSV files.")
            return


        csv1_path = save_temp_file(csv1)
        csv2_path = save_temp_file(csv2)
        plot_path = csv1_path.replace(".csv", ".png")
        if yaml_file:
            yaml_path = save_temp_file(yaml_file)
            with open(yaml_path) as f:

                chain_mapping = yaml.safe_load(f)['chain_mapping']
                st.plotly_chart(plot_distance_difference(csv1_path,csv2_path,chain_map=chain_mapping,plotly=True))

        st.success("Comparison complete.")
        