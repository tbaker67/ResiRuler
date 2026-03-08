import streamlit as st
from pathlib import Path
import re

def show_readme_tab():
    readme_path = Path(__file__).parent / "UIReadME.md"
    base_path = readme_path.parent  # for relative images

    if readme_path.exists():
        with open(readme_path, 'r') as f:
            lines = f.readlines()

        #Store markdown image syntax
        image_pattern = re.compile(r'!\[(.*?)\]\((.*?)\)')

        for line in lines:
            #Look for any markdown image syntax
            match = image_pattern.search(line)
            if match:
                #Get the stuff inside [] and actual image path
                alt_text, img_path = match.groups()
                
                #Make full image path and display in streamlit
                img_file = (base_path / img_path).resolve()
                if img_file.exists():
                    st.image(str(img_file), caption=alt_text, use_container_width=True)
                else:
                    st.warning(f"Image not found: {img_path}")
            else:
                # Display markdown for normal text lines
                if line.strip(): 
                    st.markdown(line)
    else:
        st.error("README.md not found.")