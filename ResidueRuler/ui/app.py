"""Main Streamlit application entry point for ResiRuler UI."""
import os
import sys

# Add ResidueRuler directory to path for imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import streamlit as st

from ui.tabs.align_tab import show_align_tab
from ui.tabs.compare_tab import show_compare_tab
from ui.tabs.movement_tab import show_movement_tab
from ui.tabs.readme_tab import show_readme_tab
from ui.tabs.run_tab import show_run_tab
st.set_page_config(layout="wide")
st.title("ResiRuler Prototype UI")


tab0,tab1,tab2,tab3,tab4 = st.tabs(["README","Run","Align", "Compare","Movement"])
    
with tab0:
    show_readme_tab()
#with tab1:
    #show_run_tab()
with tab2:
    show_align_tab()
with tab3:
    show_compare_tab()
with tab4: 
    show_movement_tab()