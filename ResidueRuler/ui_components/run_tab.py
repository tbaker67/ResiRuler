import streamlit as st
import pandas as pd
import matplotlib.colors as mcolors
from pathlib import Path
from ui_components.utils import chain_selector_ui, create_downloadable_zip, load_structure_if_new, get_threshold
from src.resiruler.structure_parsing import extract_residues_from_structure, convert_to_CA_coords_list
from src.resiruler.distance_calc import DistanceMatrix
from src.resiruler.plotting import plot_interactive_contact_map
from src.resiruler.chimera_export import generate_chimera_link_script
from ui_components.pymol_viewers import start_pymol_viewer, draw_links_pymol
from ui_components.color_mapping_utils import (
    get_coloring_values,
    discrete_palette_picker,
    show_gradient_bar,
    show_discrete_bar,
    build_gradient_cmap,
    sort_discrete_mapping,
    gradient_palette_picker
)
import pandas as pd
from pathlib import Path

def show_run_tab():
    cif_file = st.file_uploader("Upload Structure File (.cif)", type=["cif"])
    st.session_state.setdefault("run_clicked", False)
    st.session_state.setdefault("run_df", None)
    st.session_state.setdefault("run_contact_map", None)
    st.session_state.setdefault("run_matrix", None)
    
    st.markdown("### Select Chains to Visualize")

    structure = load_structure_if_new(cif_file, name_key="name1", struct_key="structure1")

    selected_chains = chain_selector_ui(structure)
    
    lower_threshold = get_threshold("Minimum Distance Threshold", "10.0")

    upper_threshold = get_threshold("Maximum Distance Threshold", "100.0")
    
    if st.button("Run"):
            res_list = extract_residues_from_structure(structure, selected_chains)
            coords_list, index_map = convert_to_CA_coords_list(res_list)
            matrix = DistanceMatrix(coords_list, index_map)
            df = matrix.convert_to_df(lower_threshold, upper_threshold)
            contact_map = plot_interactive_contact_map(matrix, title="Contact Map", lower_threshold=lower_threshold,
                                                       upper_threshold=upper_threshold)

            # save important data for plotting/tables
            st.session_state.run_clicked = True
            st.session_state.run_matrix = matrix
            st.session_state.run_df = df
            st.session_state.run_contact_map = contact_map
        


    # Show results if run was successful
    if st.session_state.run_clicked:
        st.markdown("### Contact Map")
        st.plotly_chart(st.session_state.run_contact_map, use_container_width=False)

        st.markdown("### Distance Data")
        #get rid of coord columns to save memory for display
        st.dataframe(st.session_state.run_df.drop(columns=['Coord1', 'Coord2']), use_container_width=True)

    st.markdown("### Edit or Paste Desired Residue Pairs")
    default_pairs = pd.DataFrame({
            'ChainID_Resnum1': ['A-58'],
            'ChainID_Resnum2': ['B-102']
        })
    


    link_selections = st.data_editor(default_pairs, num_rows="dynamic", use_container_width=True, key="residue_pairs_editor")

    st.session_state.setdefault("draw_links_clicked", False)
    st.session_state.setdefault("links_df", None)
    st.session_state.setdefault("link_script", None)
    st.session_state.setdefault("link_pymol", None)

    # Use utility to get coloring mode and value range from user
    st.title("Color Mapping")


    mode = st.radio("Coloring Mode", ["Gradient", "Discrete Mapping"], horizontal=True)

    min_val, max_val = get_coloring_values()

    if min_val >= max_val:
        st.error("Min value must be less than max value")

    if mode == "Gradient":
        default_colors = ["#44ce1b","#e51f1f"]
        # Show gradient color picker and preview
        palette = gradient_palette_picker()
        show_gradient_bar(palette, min_val, max_val)
        cmap_obj = build_gradient_cmap(palette, min_val, max_val)
    else:
        # Show discrete mapping editors and preview
        discrete_mapping = discrete_palette_picker(
        default_thresholds=[20, 30, 50],
        default_colors=["#44ce1b", "#f7e379", "#e51f1f"])
        discrete_mapping = sort_discrete_mapping(discrete_mapping)
        show_discrete_bar(discrete_mapping, min_val, max_val)

    if st.button("Draw Links"):
        # Make sure run data exists
        if "run_df" not in st.session_state or st.session_state.run_df is None or st.session_state.run_df.empty:
            st.error("Please run the analysis first!")
            return

        # Prepare augmented links accounting for reversed pairs
        reversed_links = link_selections.rename(columns={
            'ChainID_Resnum1': 'ChainID_Resnum2',
            'ChainID_Resnum2': 'ChainID_Resnum1'
        })

        run_df_clean = st.session_state.run_df.copy()
        run_df_clean['ChainID_Resnum1'] = run_df_clean['ChainID_Resnum1'].astype(str).str.strip()
        run_df_clean['ChainID_Resnum2'] = run_df_clean['ChainID_Resnum2'].astype(str).str.strip()

        link_selections = link_selections.applymap(lambda x: str(x).strip())
        augmented_links = pd.concat([link_selections, reversed_links], ignore_index=True).drop_duplicates()

        # Merge user selections with run dataframe on chain-residue pairs
        links_df = run_df_clean.merge(augmented_links, on=['ChainID_Resnum1', 'ChainID_Resnum2'])

        if mode == "Gradient":
            cmap_obj = mcolors.LinearSegmentedColormap.from_list("custom_gradient", [hex for hex, _ in palette])
            link_script = generate_chimera_link_script(
                links_df,
                chains=selected_chains,
                color_mode="gradient",
                cmap=cmap_obj,
                min_val=min_val,
                max_val=max_val
            )
            pymol_view = draw_links_pymol(
                links_df,
                start_pymol_viewer(cif_file),
                color_mode="gradient",
                cmap=cmap_obj,
                min_val=min_val,
                max_val=max_val
            )
        else:
            sorted_thresholds = sorted(discrete_mapping, key=lambda x: x[0])
            link_script = generate_chimera_link_script(
                links_df,
                chains=selected_chains,
                color_mode="discrete",
                thresholds=sorted_thresholds
            )
            pymol_view = draw_links_pymol(
                links_df,
                start_pymol_viewer(cif_file),
                color_mode="discrete",
                thresholds=sorted_thresholds
            )

        # Store and display results, offer downloads, etc.
        st.session_state.link_script = link_script
        st.session_state.draw_links_clicked = True
        st.session_state.links_df = links_df
        st.session_state.link_pymol = pymol_view

    if st.session_state.draw_links_clicked:
        st.markdown("### PyMol Visualization Preview")
        link_html = st.session_state.link_pymol._make_html()
        st.components.v1.html(link_html, height=600, width=1000)

        st.markdown("### Links Subset Data")
        st.dataframe(st.session_state.links_df)

        st.markdown("### Download ChimeraX Visualization Script")
        st.download_button("Download CXC Script",
                           data=st.session_state.link_script,
                           file_name="link_script.cxc",
                           mime="text/plain")

        chimera_filename = "chimera_link_script.cxc"
        cif_filename = Path(cif_file.name).name
        cif_content = cif_file.getvalue().decode("utf-8")
        csv_filename = "run.csv"

        files_to_zip = {
            chimera_filename: st.session_state.link_script,
            cif_filename: cif_content,
            csv_filename: st.session_state.run_df.to_csv(index=False),
        }

        zip_buffer = create_downloadable_zip(files_to_zip)
        st.session_state.run_zip_buffer = zip_buffer

        st.subheader("Download Full Data and Scripts Folder")

        st.download_button("Download All as ZIP",
                           data=st.session_state.run_zip_buffer,
                           file_name="run_analysis_package.zip",
                           mime="application/zip")


