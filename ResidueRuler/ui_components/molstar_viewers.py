import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from src.resiruler.chimera_export import get_color_discrete, get_color_gradient
import molviewspec as mvs

def create_distance_shift_builder(
    
) :
    """
    Create a Mol* builder for visualizing distance shifts between two structures.
    
    Args:
        ref_cif_data: CIF content of the reference structure (string).
        tgt_cif_data: CIF content of the target structure (string).
        ref_annotations_json: JSON string with reference residue colors.
        tgt_annotations_json: JSON string with target residue colors.
    
    Returns:
        mvs.Builder object ready for .molstar_streamlit() rendering.
    """
    builder = mvs.create_builder()
    
    # --- Reference structure ---
    structure = (
        builder.download(url='local.cif')
        .parse(format='mmcif')
        .model_structure()
    )
    cartoon = structure.component(selector='polymer').representation(type='cartoon')
    cartoon.color_from_uri(uri='annotations.json', format='json', schema='all_atomic')
    
    
    return builder

def write_movement_annotations(df, cmap=plt.cm.get_cmap("plasma"), ref=True):
    """Generate residue color annotations for Mol* viewer."""
    annotations = []
    df = df.dropna(subset=["Coord1", "Coord2", "Distance"])
    vmin, vmax = df["Distance"].min(), df["Distance"].max()
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    for _, row in df.iterrows():
        chain, resi = (
            row["ChainID_Resnum1"].split("-")
            if ref
            else row["ChainID_Resnum2"].split("-")
        )
        dist = row["Distance"]
        r, g, b = [int(255 * c) for c in cmap(norm(dist))[:3]]
        color_hex = f"#{r:02x}{g:02x}{b:02x}"
        annotations.append(
            {"auth_asym_id": chain, "auth_seq_id": resi, "color": color_hex}
        )
    return annotations