import tempfile
import os
import streamlit as st
import json
import re 
import base64
from io import BytesIO
import zipfile
from src.resiruler.auto_alignment import StructureMapper
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

def create_mapper(cif1, cif2, chain_mapping):
    if not cif1 or not cif2:
        raise ValueError("Both cif1 and cif2 must be provided.")
    cif1_path = save_temp_file(cif1)
    cif2_path = save_temp_file(cif2)

    structure1 = load_structure(cif1_path)
    structure2 = load_structure(cif2_path)

    mapper = StructureMapper(structure1, structure2)

    if chain_mapping:
        mapper.map_chains_explicit(chain_mapping)
    else:
        mapper.map_chains(threshold=95)
    
    return mapper