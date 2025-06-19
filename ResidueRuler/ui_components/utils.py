import tempfile
import os
import streamlit as st
import json

def save_temp_file(uploaded_file):
    with tempfile.NamedTemporaryFile(delete=False, suffix=os.path.splitext(uploaded_file.name)[-1]) as tmp_file:
        tmp_file.write(uploaded_file.getbuffer())
        return tmp_file.name
    
def json_mapping_input(label, default):
    st.write("Enter JSON mapping")
    if isinstance(default, dict):
        default_json = json.dumps(default, indent=2)
    elif isinstance(default, str):
        default_json = default
    else:
        default_json = '{}'

    user_input = st.text_area(label, height=200, value=default_json)

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