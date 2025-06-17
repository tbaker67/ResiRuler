import sys
import os
sys.path.insert(0,os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import streamlit as st
from run_tab import show_run_tab
from compare_tab import show_compare_tab
from movement_tab import show_movement_tab
st.set_page_config(layout="wide")
st.title("ResiRuler Prototype UI")

mode = st.sidebar.radio("Choose Mode", ["Run", "Compare", "Movement"])

if mode == "Run":
    show_run_tab()
elif mode == "Compare":
    show_compare_tab()
elif mode == "Movement":
    show_movement_tab()