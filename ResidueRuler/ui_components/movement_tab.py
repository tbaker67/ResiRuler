import streamlit as st
from ui_components.utils import save_temp_file, json_mapping_input
from src.resiruler import load_structure
from src.resiruler.distance_calc import calc_difference_aligned
from src.resiruler.plotting import plot_movement_shift, plot_movement_vectors
from src.resiruler.chimera_export import chimera_color_shift_from_csv, chimera_movement_vectors_from_csv
import pandas as pd
import yaml

def show_movement_tab():
    st.header("Movement Analysis Between Aligned Structures")

    cif1 = st.file_uploader("Upload Aligned CIF #1", type=["cif"], key="aligned1")
    cif2 = st.file_uploader("Upload Aligned CIF #2", type=["cif"], key="aligned2")
    yaml_file = st.file_uploader("Upload Chain Mapping YAML (optional)", type=["yaml", "yml"])

    if st.button("Analyze Movement"):
        if not cif1 or not cif2:
            st.error("Please upload both structure files.")
            return

        cif1_path = save_temp_file(cif1)
        cif2_path = save_temp_file(cif2)

        structure1 = load_structure(cif1_path)
        structure2 = load_structure(cif2_path)

        chain_mapping = None
        if yaml_file:
            yaml_path = save_temp_file(yaml_file)
            with open(yaml_path) as f:
                chain_mapping = yaml.safe_load(f)['chain_mapping']

        df = calc_difference_aligned(structure1, structure2, chain_mapping)

        st.success("Movement calculated.")
        st.dataframe(df)

        st.plotly_chart(plot_movement_shift(df,plotly=True))

        st.plotly_chart(plot_movement_vectors(df,plotly=True))
        

        csv_output_path = cif1.name.replace(".cif", "_movement.csv")
        
        st.download_button("Download CSV", data=df.to_csv(index=False), file_name=csv_output_path)

        chimera_color_shift_from_csv(csv_output_path, csv_output_path.replace(".csv", ".cxc"), chain_mapping)
        chimera_movement_vectors_from_csv(csv_output_path, csv_output_path.replace(".csv", ".bild"))