"""Utility functions and UI widgets for the ResiRuler Streamlit interface."""
import io
import os
import tempfile
import zipfile
from contextlib import contextmanager
import numpy as np
import pandas as pd
import streamlit as st
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.PDB.mmcifio import MMCIFIO

from src.resiruler.core.auto_alignment import EnsembleMapper, StructureMapper
from src.resiruler.core.structure_parsing import extract_res_from_chain, load_structure
from src.resiruler.viz.plotting import (
    plot_comparison_with_contact_filter,
    plot_contacts_gained,
    plot_contacts_lost,
    plot_interactive_contact_map,
)

@contextmanager 
def struct_to_temp_cif(structure):
    """
    Create a temporary CIF file for a Bio.PDB structure.
    Automatically deletes the file after use.
    """
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".cif")
    try:
        io = MMCIFIO()
        io.set_structure(structure)
        io.save(tmp.name)
        tmp.close()
        yield tmp.name
    finally:
        if os.path.exists(tmp.name):
            os.unlink(tmp.name)

def save_temp_file(uploaded_file):
    with tempfile.NamedTemporaryFile(delete=False, suffix=os.path.splitext(uploaded_file.name)[-1]) as tmp_file:
        tmp_file.write(uploaded_file.getbuffer())
        return tmp_file.name


def create_downloadable_zip(files_dict):
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, 'w') as z:
        for filename, content in files_dict.items():
            z.writestr(filename, content)
    buffer.seek(0)
    return buffer

def create_ensemble_mapper(ref_structure, tgt_structures, chain_mappings, threshold, protein_aligner, nucleotide_aligner):
    ensemble_mapper = EnsembleMapper(ref_structure, protein_aligner=protein_aligner, nucleotide_aligner=nucleotide_aligner)

    for structure_name, tgt_structure in tgt_structures.items():
        ensemble_mapper.add_structure(structure_name, tgt_structure,threshold, chain_mappings[structure_name])
    
    return ensemble_mapper

def chain_selector_ui(structure, label="Select Chains to Visualize (Chains with All will be deselected)", default_all=True, key_prefix="chain_selector"):
    """
    Display a Streamlit multiselect box for choosing chains from a structure or structure mapper.
    If "All" is selected along with individual chains, returns all chains EXCEPT those individual chains.
    As such, structure should only be a BioPython structure object, or a StructureMapper object from autoalignments.py
    """
    if structure is None:
        st.warning("No structure loaded. Please upload a file first.")
        return None
    
    chain_options = None
    if isinstance(structure, StructureMapper):
        chain_options = sorted(structure.matched_ref_chains)
    else:
        chain_options = sorted({chain.id for model in structure for chain in model})
    
    multiselect_options = ["All"] + chain_options
    default_selection = ["All"] if default_all else []
    selected = st.multiselect(label, multiselect_options, default=default_selection, key=key_prefix)
    
    if "All" in selected and len(selected) > 1:
        excluded_chains = [c for c in selected if c != "All"]
        selected_chains = [c for c in chain_options if c not in excluded_chains]
        return selected_chains if selected_chains else None
    elif "All" in selected:
        return chain_options
    elif selected:
        return selected
    else:
        return None
    
def load_structure_if_new(cif_file, name_key, struct_key):
    if cif_file is None:
        return None

    if (name_key not in st.session_state 
        or st.session_state[name_key] != cif_file.name):
        
        cif_path = save_temp_file(cif_file)
        structure = load_structure(cif_path)


        st.session_state[struct_key] = structure
        st.session_state[name_key] = cif_file.name
    else:
        # same file
        structure = st.session_state.get(struct_key)

    return structure

def load_structures_if_new(cif_files, name_key_prefix, struct_key_prefix):
    """
    Load one or more structures, caching them in session_state.
    Returns a dict mapping file name -> structure object.
    """
    if not cif_files:
        return {}

    structures = {}
    for i, cif_file in enumerate(cif_files):
        name_key = f"{name_key_prefix}_{i}"
        struct_key = f"{struct_key_prefix}_{i}"

        if (name_key not in st.session_state 
            or st.session_state[name_key] != cif_file.name):

            cif_path = save_temp_file(cif_file)
            structure = load_structure(cif_path)

            st.session_state[struct_key] = structure
            st.session_state[name_key] = cif_file.name
        else:
            structure = st.session_state.get(struct_key)

        #remove .cif for names
        structures[cif_file.name[:-4]] = structure

    return structures


def chain_mapping_input(ref_chains, tgt_chains, default=None, key="chain_mapping"):
    """
    Editable table for mapping reference → target chains.
    - Blank entries are ignored.
    - if all entries blank returns None (auto mapping mode).
    """
    st.caption("Map reference chains → target chains (leave all blank for auto mapping)")

    # Build default DataFrame
    with st.expander("Explicit Chain Mapping Options"):
        df = pd.DataFrame(
            [(ref, "","") for ref in ref_chains],
            columns=["Reference Chain", "Target Chain", "Chain Type"]
            )

        # Editable table
        edited_df = st.data_editor(
            df,
            num_rows="fixed",
            key=key,
            column_config={
                "Reference Chain": st.column_config.TextColumn("Reference Chain", disabled=True),
                
                "Target Chain": st.column_config.SelectboxColumn(
                    "Target Chain",
                    options=[""] + list(tgt_chains),  # allow blank
                    required=False,
                ),

                "Chain Type": st.column_config.SelectboxColumn(
                    "Chain Type",
                    options=[""] + ["protein", "dna", "rna"]
                ),
            },
            use_container_width=True,
            hide_index=True,
        )

        # Convert to dict, skipping blanks
        mapping = {
            ref: (tgt, type) 
            for ref, tgt, type in zip(edited_df["Reference Chain"], edited_df["Target Chain"], edited_df['Chain Type'])
            if tgt  # only keep non-blank
        }

        # If user left everything blank → return None (auto mode)
        if not mapping:
            st.info("No explicit mapping specified → auto mapping will be used.")
            return None

        return mapping


def get_chain_mappings_for_targets(structures_dict, ref_chains, default_mapping=None, key="mapping"):
    """
    Build chain mapping table for each target structure.
    """
    mappings = {}
    for filename, tgt_structure in structures_dict.items():
        st.subheader(f"Chain mapping for target: {filename}")

        tgt_chains = [chain.id for chain in tgt_structure.get_chains()]  
        mapping = chain_mapping_input(
            ref_chains,
            tgt_chains,
            default=default_mapping,
            key=f"{key}_{filename}"
        )
        mappings[filename] = mapping

    return mappings

def get_threshold(label, default, key):
    threshold_input = st.text_input(label, value=default, key=key)
    try:
        return float(threshold_input)
    except ValueError:
        st.error("Please enter a valid number.")
        return None


def create_downloadable_zip_grouped(file_groups):
    """
    Create a ZIP buffer from multiple file groups.
    file_groups = dict where
        key   = folder name or None for root
        value = dict mapping {filename: content}
    
    Special case:
        key=None → files go in the root of the ZIP
        key="csv"/"bild"/"models" → files go in that subfolder
    Supports str, bytes, StringIO, BytesIO as content.
    """
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, 'w', zipfile.ZIP_DEFLATED) as z:
        for folder, files in file_groups.items():
            for filename, content in files.items():
                # Convert buffers to bytes
                if hasattr(content, "getvalue"):
                    content = content.getvalue()
                if isinstance(content, str):
                    content = content.encode("utf-8")
                # Determine path inside ZIP
                arcname = f"{folder}/{filename}" if folder else filename
                z.writestr(arcname, content)
    buffer.seek(0)
    return buffer


PROTEIN_SUBSTITUTION_MATRICES = [
    "BLOSUM62","BLOSUM45", "BLOSUM50", "BLOSUM89", "BLOSUM90", "BLASTP", "DAYHOFF", "FENG", "GENETIC", "GONNET1992", "JOHNSON", "JONES","LEVIN","MCLACHLAN", "MDM78", 
    "BENNER22","BENNER6", "BENNER74", "PAM250", "PAM30", "PAM70", "RAO", "RISLER"
]

NUCLEOTIDE_SUBSTITUTION_MATRICES = [
    "MEGABLAST", "BLASTN", "NUC.4.4", "HOXD70"
]


def aligner_ui(protein, key_prefix="aligner"):
    """
    Function to allow user to specify a biopython Pairwise Aligner Object
    """
    with st.container():
        # Alignment mode and matrix 
        col1, col2 = st.columns(2)
        with col1:
            mode = st.selectbox("Alignment mode", ["global", "local"], key=f"{key_prefix}_mode")
        with col2:
            if protein:
                matrix_name = st.selectbox("Substitution matrix", PROTEIN_SUBSTITUTION_MATRICES , key=f"{key_prefix}_matrix")
            else:
                matrix_name = st.selectbox("Substitution matrix", NUCLEOTIDE_SUBSTITUTION_MATRICES , key=f"{key_prefix}_matrix")

        # Gap scores in one row
        col3, col4 = st.columns(2)
        with col3:
            open_gap_score = st.number_input("Gap opening penalty", value=-10, step=1, key=f"{key_prefix}_open_gap")
        with col4:
            extend_gap_score = st.number_input("Gap extension penalty", value=-1, step=1, key=f"{key_prefix}_extend_gap")

        # optional end-gap penalties 
        with st.expander("Advanced: End-gap penalties"):
            use_end_gaps = st.checkbox("Enable end-gap penalties?", key=f"{key_prefix}_use_end_gaps")
            left_open_gap_score = None
            right_open_gap_score = None
            if use_end_gaps:
                col5, col6 = st.columns(2)
                with col5:
                    left_open_gap_score = st.number_input("Left end gap penalty", value=0, step=1, key=f"{key_prefix}_left_end_gap")
                with col6:
                    right_open_gap_score = st.number_input("Right end gap penalty", value=0, step=1, key=f"{key_prefix}_right_end_gap")

    
    aligner = PairwiseAligner()
    aligner.mode = mode
    aligner.substitution_matrix = substitution_matrices.load(matrix_name)
    aligner.open_gap_score = open_gap_score
    aligner.extend_gap_score = extend_gap_score
    if use_end_gaps:
        aligner.left_open_gap_score = left_open_gap_score
        aligner.right_open_gap_score = right_open_gap_score

    return aligner

def show_alignments(ensemble_mapper, key="alignment_select"):
    """
    shows the alignment present in all structure mappers of a selected chiain in the reference
    """
    selected_chain = st.selectbox(
        "Select Chain to view alignment of",
        options=[chain.id for chain in ensemble_mapper.ref_structure.get_chains()],
        key=key
    )

    for structure_name, mapper in ensemble_mapper.structure_mappings.items():
        chain_mapper = mapper.chain_mappings.get(selected_chain, None)
        if chain_mapper is None or chain_mapper.alignment is None:
            continue

        
        st.subheader(f"Alignment for structure: {structure_name}, chain: {selected_chain}")
        
        show_chain_alignment(chain_mapper)


def show_chain_alignment(chain_mapper):
    """
    Creates html for an alignment viewer, with highlighted colors indicating mismatches(red), gaps(grey), and matches(green)
    additionally, by hovering over a residue in the sequence the residue number will be shown
    """
    if chain_mapper is None or chain_mapper.alignment is None:
        st.warning("No alignment available")
        return
    


    seqs = [seq for seq in chain_mapper.alignment]
    if len(seqs) < 2:
        st.warning("Alignment has less than 2 sequences")
        return

    ref_seq = chain_mapper.aligned_ref_seq
    tgt_seq = chain_mapper.aligned_tgt_seq
    ref_res = extract_res_from_chain(chain_mapper.ref_chain)
    tgt_res = extract_res_from_chain(chain_mapper.tgt_chain)

    idx_ref, idx_tgt = 0, 0

    ref_line, tgt_line = "", ""

    for a, b in zip(ref_seq, tgt_seq):
        #color based on matches
        color = "lightgray" if a == "-" or b == "-" else "lightgreen" if a == b else "lightcoral"

        # reference residue
        if a != "-":
            res = ref_res[idx_ref]
            ref_num = f"{res.id[1]}{res.id[2].strip()}"  # includes insertion code if present
            idx_ref += 1
            ref_line += f"<span class='residue' style='background-color:{color}'>" \
                        f"{a}<span class='tooltip'>{ref_num}</span></span>"
        else:
            ref_line += f"<span class='residue' style='background-color:{color}'>{a}</span>"

        # target residue
        if b != "-":
            res = tgt_res[idx_tgt]
            tgt_num = f"{res.id[1]}{res.id[2].strip()}"
            idx_tgt += 1
            tgt_line += f"<span class='residue' style='background-color:{color}'>" \
                        f"{b}<span class='tooltip'>{tgt_num}</span></span>"
        else:
            tgt_line += f"<span class='residue' style='background-color:{color}'>{b}</span>"

    #Creates an html object that can be viewed in streamlit
    html_content = f"""
    <html>
    <head>
    <style>
    .residue {{
        position: relative;
        display: inline-block;
        width: 1.3ch;  /* fixed width per residue */
        font-family: monospace;
        text-align: center;
        cursor: pointer;
    }}
    .tooltip {{
        display: none;
        position: absolute;
        top: -0.8em;
        left: 50%;
        transform: translateX(-50%);
        background: #eee;
        border: 1px solid #999;
        font-size: 10px;
        padding: 1px 3px;
        border-radius: 3px;
        white-space: nowrap;
        z-index: 10;
    }}
    .residue:hover .tooltip {{
        display: block;
    }}
    .sequence {{
        white-space: pre;
    }}
    .label {{
        display: inline-block;
        width: 80px;   /* fixed width for both labels */
        text-align: right;
        font-weight: bold;
    }}
    </style>
    </head>
    <body>
        <div class='sequence'><span class='label'>Ref ({chain_mapper.ref_chain.id}):</span> {ref_line}</div>
        <div class='sequence'><span class='label'>Tgt ({chain_mapper.tgt_chain.id}):</span> {tgt_line}</div>
    </body>
    </html>
    """

    st.components.v1.html(html_content, height=50, scrolling=True)

def get_measurement_mode(key="measurement_mode"):
    st.markdown("Choose Where to Measure Distances From")
    protein_option_to_mode = {"C-α":"CA", "C-β":"CB", "Side-Chain Centroid (Heavy Atoms Excluding CA)":"SC"}
    nucleic_option_to_mode = {"C1'":"C1'", "Centroid ()":"NC"}

    col1, col2 = st.columns(2)
    with col1:
        selected_protein_mode = st.selectbox("Choose Which Atoms To Measure Proteins From", protein_option_to_mode.keys(), key=f"{key}_protein")
    with col2:
        selected_nucleic_mode = st.selectbox("Choose Which Atoms To Measure Nucleic Acids From", nucleic_option_to_mode.keys(), key=f"{key}_nucleic")

    return protein_option_to_mode[selected_protein_mode], nucleic_option_to_mode[selected_nucleic_mode]

def filter_df_by_chains(df, selected_chains):
    """Filter dataframe to only include residues from selected chains in the reference structure."""
    if df.empty or not selected_chains:
        return df
    chain1 = df['ChainID_Resnum1'].str.split('-').str[0]
    mask = chain1.isin(selected_chains)
    return df[mask].reset_index(drop=True)

def distance_threshold_ui(key_prefix="distance_threshold"):
    """
    Function to allow user input for distance thresholding with a double-ended slider.
    """
    enable_threshold = st.checkbox(
        "Enable contact filtering",
        value=False,
        help="Filter displayed contacts by distance range. Shows shared contacts, gained, and lost.",
        key=f"{key_prefix}_enable"
    )
    
    if not enable_threshold:
        return False, None, None
    
    st.markdown("**Distance Range (Å)**")
    
    slider_values = st.slider(
        "Distance range",
        min_value=0.0, max_value=100.0,
        value=(0.0, 15.0),
        step=0.5,
        key=f"{key_prefix}_range_slider",
        label_visibility="collapsed"
    )
    
    col_low, col_high = st.columns(2)
    with col_low:
        lower_num = st.number_input(
            "Lower (Å)",
            min_value=0.0, max_value=100.0,
            value=slider_values[0],
            step=0.5,
            key=f"{key_prefix}_lower_num"
        )
    with col_high:
        upper_num = st.number_input(
            "Upper (Å)",
            min_value=0.0, max_value=200.0,
            value=slider_values[1],
            step=0.5,
            key=f"{key_prefix}_upper_num"
        )
    
    # Prefer numeric input
    lower_threshold = lower_num if lower_num != slider_values[0] else slider_values[0]
    upper_threshold = upper_num if upper_num != slider_values[1] else slider_values[1]
    
    if lower_threshold >= upper_threshold:
        st.warning("Lower threshold must be less than upper threshold.")
    
    st.caption(f"Showing contacts with distances between {lower_threshold:.1f}Å and {upper_threshold:.1f}Å")
    
    return enable_threshold, lower_threshold, upper_threshold

def get_chain_pair_options(index_map):
    """Generate list of chain pair options from index_map."""
    chains = sorted(set(chain for chain, _ in index_map.keys()))
    pairs = []
    for i, c1 in enumerate(chains):
        for c2 in chains[i:]:  # Include self-pairs and upper triangle
            if c1 == c2:
                pairs.append((c1, c2, f"{c1} (intra-chain)"))
            else:
                pairs.append((c1, c2, f"{c1} × {c2}"))
    return chains, pairs


def display_chain_pair_selector(ref_dm, tgt_dm, compare_dm, selected_target, 
                                 lower_threshold, upper_threshold, contact_threshold=None):
    """Display chain-pair selector for large structures."""
    chains, pairs = get_chain_pair_options(ref_dm.index_map)
    
    st.markdown("### Select Chain Interactions to Display")
    st.caption("For large structures, select specific chain pairs to view interactively.")
    pair_labels = [p[2] for p in pairs]
    default_selection = pair_labels[:min(3, len(pair_labels))]
    
    selected_pairs = st.multiselect(
        "Chain pairs",
        pair_labels,
        default=default_selection,
        help="Select which chain interactions to display. Each pair will show as a separate heatmap.",
        key="chain_pair_selector"
    )
    
    if not selected_pairs:
        st.warning("Select at least one chain pair to display.")
        return
    
    if st.button("Show Selected Chain Pairs", key="show_chain_pairs"):
        # Map labels back to chain IDs
        label_to_chains = {p[2]: (p[0], p[1]) for p in pairs}
        
        for pair_label in selected_pairs:
            c1, c2 = label_to_chains[pair_label]
            chain_ids = [c1] if c1 == c2 else [c1, c2]
            st.markdown(f"#### {pair_label}")
            
            try:
                ref_sub = ref_dm.get_submatrix(chain_ids=chain_ids)
                tgt_sub = tgt_dm.get_submatrix(chain_ids=chain_ids)
                compare_sub = compare_dm.get_submatrix(chain_ids=chain_ids)
                
                sub_size = ref_sub.get_size()
                st.caption(f"{sub_size} residues")
                
                ref_mat_filtered = ref_sub.mat
                tgt_mat_filtered = tgt_sub.mat
                if lower_threshold is not None and upper_threshold is not None:
                    ref_mat_filtered = np.where((ref_sub.mat > lower_threshold) & (ref_sub.mat < upper_threshold), ref_sub.mat, np.nan)
                    tgt_mat_filtered = np.where((tgt_sub.mat > lower_threshold) & (tgt_sub.mat < upper_threshold), tgt_sub.mat, np.nan)
                
                shared_min = min(np.nanmin(ref_mat_filtered), np.nanmin(tgt_mat_filtered))
                shared_max = max(np.nanmax(ref_mat_filtered), np.nanmax(tgt_mat_filtered))
                
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    ref_fig = plot_interactive_contact_map(
                        ref_sub, lower_threshold=lower_threshold, upper_threshold=upper_threshold,
                        title=f"Reference",
                        min=shared_min, max=shared_max
                    )
                    ref_fig.update_layout(width=400, height=400)
                    st.plotly_chart(ref_fig, use_container_width=True)
                
                with col2:
                    tgt_fig = plot_interactive_contact_map(
                        tgt_sub, lower_threshold=lower_threshold, upper_threshold=upper_threshold,
                        title=f"Target",
                        min=shared_min, max=shared_max
                    )
                    tgt_fig.update_layout(width=400, height=400)
                    st.plotly_chart(tgt_fig, use_container_width=True)
                
                with col3:
                    if contact_threshold is not None:
                        compare_fig = plot_comparison_with_contact_filter(
                            compare_sub, contact_threshold=contact_threshold,
                            title=f"Δ Distance (shared)"
                        )
                    else:
                        compare_fig = plot_interactive_contact_map(
                            compare_sub, lower_threshold=lower_threshold, upper_threshold=upper_threshold,
                            title=f"Δ Distance"
                        )
                    compare_fig.update_layout(width=400, height=400)
                    st.plotly_chart(compare_fig, use_container_width=True)
                
                if contact_threshold is not None:
                    col_g, col_l = st.columns(2)
                    with col_g:
                        gained_fig = plot_contacts_gained(compare_sub, contact_threshold=contact_threshold)
                        gained_fig.update_layout(width=400, height=400)
                        st.plotly_chart(gained_fig, use_container_width=True)
                    with col_l:
                        lost_fig = plot_contacts_lost(compare_sub, contact_threshold=contact_threshold)
                        lost_fig.update_layout(width=400, height=400)
                        st.plotly_chart(lost_fig, use_container_width=True)
                
                with st.expander(f"View {pair_label} data table"):
                    pair_df = compare_sub.convert_to_df()
                    st.dataframe(pair_df, use_container_width=True)
                    
            except ValueError as e:
                st.error(f"Could not display {pair_label}: {e}")
