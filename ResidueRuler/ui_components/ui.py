import sys
import os
from pathlib import Path
sys.path.insert(0,os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import streamlit as st
from align_tab import show_align_tab
from run_tab import show_run_tab
from compare_tab import show_compare_tab
from movement_tab import show_movement_tab
from readme_tab import show_readme_tab
st.set_page_config(layout="wide")
st.title("ResiRuler Prototype UI")


tab0,tab1,tab2,tab3,tab4 = st.tabs(["README","Run","Align", "Compare","Movement"])
    
with tab0:
    show_readme_tab()
with tab1:
    show_run_tab()
with tab2:
    show_align_tab()
with tab3:
    show_compare_tab()
with tab4: 
    show_movement_tab()