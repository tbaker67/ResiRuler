"""Main Streamlit application entry point for ResiRuler UI."""
import os
import sys

# Add ResidueRuler directory to path for imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import streamlit as st

from ui.tabs.align_tab import show_align_tab
from ui.tabs.compare_tab import show_compare_tab
from ui.tabs.displacement_tab import show_displacement_tab
st.set_page_config(layout="wide")
st.title("ResiRuler")

tab0,tab1,tab2 = st.tabs(["Align", "Compare","Movement"])
with tab0:
    show_align_tab()
with tab1:
    show_compare_tab()
with tab2: 
    show_displacement_tab()
