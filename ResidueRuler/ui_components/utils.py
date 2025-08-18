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
from src.resiruler.structure_parsing import load_structure

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


def get_chain_mappings_for_targets(structures_dict, default_mapping):
    """
    Ask user for chain mapping JSON for each target structure.
    Returns dict {filename: mapping_dict}.
    """
    mappings = {}
    for filename in structures_dict.keys():
        label = f"Chain mapping for target: {filename}"
        mapping = json_mapping_input(label, default_mapping, key=f"mapping_{filename}")
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
