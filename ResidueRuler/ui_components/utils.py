import tempfile
import os
import streamlit as st
import json
import re 
import base64
import io
from io import BytesIO
import zipfile
from src.resiruler.auto_alignment import StructureMapper, EnsembleMapper
from src.resiruler.structure_parsing import load_structure, extract_res_from_chain

def save_temp_file(uploaded_file):
    with tempfile.NamedTemporaryFile(delete=False, suffix=os.path.splitext(uploaded_file.name)[-1]) as tmp_file:
        tmp_file.write(uploaded_file.getbuffer())
        return tmp_file.name
    
def json_mapping_input(label, default,key):
    st.write("Enter JSON mapping")
    if isinstance(default, dict):
        default_json = json.dumps(default, indent=2)
    elif isinstance(default, str):
        default_json = default
    else:
        default_json = '{}'

    user_input = st.text_area(label, height=200, value=default_json, key=key)

    try:
        mapping = json.loads(user_input)
        if isinstance(mapping, dict):
            st.success("JSON parsed successfully!")
            st.json(mapping)
            return mapping
        else:
            st.error("Input must be a JSON object (dictionary).")
    except json.JSONDecodeError as e:
        st.error(f"Invalid JSON: {e.msg}")

def embed_local_images(markdown_text, base_path):
    """
    Finds all ![alt](relative_path) and embeds them as base64-encoded images.
    """
    def replacer(match):
        alt_text, rel_path = match.groups()
        image_path = (base_path / rel_path).resolve()
        if image_path.exists():
            # Guess content type based on extension
            ext = image_path.suffix.lower().lstrip(".")
            mime = f"image/{'jpeg' if ext == 'jpg' else ext}"
            img_bytes = image_path.read_bytes()
            b64 = base64.b64encode(img_bytes).decode()
            return f'![{alt_text}](data:{mime};base64,{b64})'
        else:
            return f"*Image not found: {rel_path}*"

    # Replace all image links in markdown
    return re.sub(r'!\[(.*?)\]\((.*?)\)', replacer, markdown_text)

def create_downloadable_zip(files_dict):
    buffer = BytesIO()
    with zipfile.ZipFile(buffer, 'w') as z:
        for filename, content in files_dict.items():
            z.writestr(filename, content)
    buffer.seek(0)
    return buffer

def create_mapper(structure1, structure2, chain_mapping, threshold):
    if not structure1 or not structure2:
        raise ValueError("Missing a structure.")
    
    mapper = StructureMapper(structure1, structure2)

    if chain_mapping:
        mapper.map_chains_explicit(chain_mapping)
    else:
        mapper.map_chains(threshold=threshold)
    
    return mapper

def create_ensemble_mapper(ref_structure, tgt_structures, chain_mappings, threshold, aligner):
    ensemble_mapper = EnsembleMapper(ref_structure, aligner)

    for structure_name, tgt_structure in tgt_structures.items():
        ensemble_mapper.add_structure(structure_name, tgt_structure,threshold, chain_mappings[structure_name])
    
    return ensemble_mapper

def chain_selector_ui(structure, label="Select Chains to Visualize", default_all=True):
    """
    Display a Streamlit multiselect box for choosing chains from a structure or structure mapper
    As such, structure should only be a BioPython structure object, or a StructureMapper object from autoalignments.py
    """
    if structure is None:
        st.warning("No structure loaded. Please upload a file first.")
        return None
    
    chain_options = None
    if isinstance(structure, StructureMapper):
        chain_options = sorted(structure.matched_ref_chains)
        #
    else:
        chain_options = sorted({chain.id for model in structure for chain in model})
    multiselect_options = ["All"] + chain_options

    default_selection = ["All"] if default_all else []
    selected = st.multiselect(label, multiselect_options, default=default_selection)

    if "All" in selected:
        return chain_options
    elif selected:
        return selected
    else:
        return None
    
def load_structure_if_new(cif_file, name_key, struct_key):
    if cif_file is None:
        return None  # Nothing uploaded


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

import pandas as pd

def chain_mapping_input(ref_chains, tgt_chains, default=None, key="chain_mapping"):
    """
    Editable table for mapping reference → target chains.
    - Blank entries are ignored.
    - if all entries blank returns None (auto mapping mode).
    """
    st.caption("Map reference chains → target chains (leave all blank for auto mapping)")

    # Build default DataFrame
    with st.expander("Explicit Chain Mapping Options"):
        if isinstance(default, dict):
            df = pd.DataFrame(
                [(ref, default.get(ref, "")) for ref in ref_chains],
                columns=["Reference Chain", "Target Chain"]
            )
        else:
            df = pd.DataFrame(
                [(ref, "") for ref in ref_chains],
                columns=["Reference Chain", "Target Chain"]
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
            },
            use_container_width=True,
            hide_index=True,
        )

        # Convert to dict, skipping blanks
        mapping = {
            ref: tgt
            for ref, tgt in zip(edited_df["Reference Chain"], edited_df["Target Chain"])
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

def get_threshold(label, default):
    threshold_input = st.text_input(label, value=default)
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

from Bio.Align import PairwiseAligner, substitution_matrices

def aligner_ui(key_prefix="aligner"):
    """
    Function to allow user to specify a biopython Pairwise Aligner Object
    """
    st.subheader("Pairwise Aligner Settings")

    
    with st.container():
        # Alignment mode and matrix 
        col1, col2 = st.columns(2)
        with col1:
            mode = st.selectbox("Alignment mode", ["global", "local"], key=f"{key_prefix}_mode")
        with col2:
            matrix_name = st.selectbox("Substitution matrix", ["BLOSUM62", "BLOSUM80", "PAM30", "PAM250"], key=f"{key_prefix}_matrix")

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
        
        show_chain_alignment(chain_mapper, ref_start=chain_mapper.ref_start, tgt_start=chain_mapper.tgt_start)


def show_chain_alignment(chain_mapper, ref_start=1, tgt_start=1):
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
    option_to_mode = {"C-α":"CA", "C-β":"CB", "Side-Chain Centroid (Heavy Atoms Excluding CA)":"SC"}
    selected_mode = st.selectbox(
        "Choose Where to Measure Distances From",
        options=["C-α", "C-β", "Side-Chain Centroid (Heavy Atoms Excluding CA)"],
        key=key
    )
    return option_to_mode[selected_mode]