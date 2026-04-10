"""Export utilities for ChimeraX scripts, BILD files, and other visualizations."""
import ast
import io
import os
from os.path import basename
import math
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd


def safe_eval(val):
    """
    Extract the values from a dataframe safely and handle possibility of them
    being stored as strings or already-parsed objects.
    """
    # handle numpy arrays / lists before isna check
    if isinstance(val, (np.ndarray, list, tuple)):
        return val
    if pd.isna(val):
        return np.nan
    if isinstance(val, (int, float, dict)):
        return val
    try:
        return ast.literal_eval(val)
    except (ValueError, SyntaxError):
        return np.nan

def get_color_discrete(distance, thresholds):
    if thresholds is None:
            raise ValueError("Discrete mode requires thresholds")
    for threshold, color_hex in thresholds:
        if distance <= threshold:
            return color_hex
    # If distance above all thresholds, return last color
    return thresholds[-1][1]

def get_color_gradient(distance, cmap, min_val, max_val):
    if None in (cmap, min_val, max_val):
        raise ValueError("Gradient mode requires cmap, min_val, and max_val")
    norm = mcolors.Normalize(vmin=min_val, vmax=max_val)
    rgba = cmap(norm(distance))  # tuple with floats (r,g,b,a)
    r, g, b = (int(255 * c) for c in rgba[:3])
    return f"#{r:02x}{g:02x}{b:02x}"


def generate_chimera_link_script(df, chains = None, color_mode="discrete", **kwargs):
    output = io.StringIO()
    for _, row in df.iterrows():
        coord1, coord2, dist = row['Coord1'], row['Coord2'], row['Distance']
        if any(pd.isna([coord1, coord2, dist])):
            continue
        color = None
        if color_mode == "discrete":
            color = get_color_discrete(dist, **kwargs)
        else:
            color = get_color_gradient(dist, **kwargs)
        output.write(
            f"shape cylinder radius 1 fromPoint {coord1[0]},{coord1[1]},{coord1[2]} "
            f"toPoint {coord2[0]},{coord2[1]},{coord2[2]} color {color}\n"
        )

    if chains:
        output.write("hide\n")
        for chain in chains:
            output.write(f"cartoon /{chain}\n")

    return output.getvalue()

def generate_chimera_coloring_palette_string(palette, positions):

    palette_string = ""
    for color, position in zip(palette, positions):
        palette_string += f"{position:.2f},{color}:"
    
    return palette_string[:-1]


def generate_pml_palette_string(palette):
    palette_string = ""
    for hex_color in palette:
        #regular hex codes have #00008B, for pymol it needs to be in the form 0x00008B
        pymol_color = hex_color.replace("#", "0x")
        palette_string += f"{pymol_color} "
    return palette_string


def generate_multiple_movement_scripts(movement_dfs, ref_name, palette, positions):

    full_def_attr = io.StringIO()
    full_cxc_script = io.StringIO()
    full_pml_script = io.StringIO()

    ids = 1
    for tgt_structure_name, movement_df in movement_dfs.items():
        #Assumes that the names are simply the cif files without the extension
        ref_cif_name = ref_name + ".cif"
        tgt_cif_name = tgt_structure_name + ".cif"
        def_attr, cxc_script, pml_script = generate_shift_scripts(movement_df, ref_cif_name, ref_name, tgt_cif_name, tgt_structure_name, ids)

        full_def_attr.write(def_attr + '\n')
        full_cxc_script.write(cxc_script + '\n')
        full_pml_script.write(pml_script + '\n')

        ids += 2

    chimera_coloring_palette_string = generate_chimera_coloring_palette_string(palette, positions)
    full_cxc_script.write("open full_defattr.defattr \n")
    full_cxc_script.write(f"color #{1}-{ids - 1} grey \n")
    full_cxc_script.write(f"color byattribute r:distance #{1}-{ids - 1} target scab palette {chimera_coloring_palette_string}\n")
    chimera_coloring_palette_string = chimera_coloring_palette_string.split(",", 1)[1]
    chimera_key_string = chimera_coloring_palette_string.replace(",", " ") + ":"
    full_cxc_script.write(f"key {chimera_key_string }\n")

    pml_palette_string = generate_pml_palette_string(palette)
    full_pml_script.write(f'spectrum properties["distance"], {pml_palette_string}')

    return full_def_attr.getvalue(), full_cxc_script.getvalue(), full_pml_script.getvalue()

def generate_shift_scripts(df, cif1_name, structure_name1, cif2_name, structure_name2, first_structure_id):
    """
    Generate defattr files, a bild file, a cxc chimera script and a pml script to color models corresponding to distance between corresponding residues in the reference and target structures
    """
    distances = df['Distance'].apply(safe_eval)
    ids_ref = df['ChainID_Resnum1']
    ids_tgt = df['ChainID_Resnum2']


    name1 = structure_name1
    name2 = structure_name2

    pml = io.StringIO()
    pml.write(f"load models/{cif1_name}, {name1}-{name2} \n")
    pml.write(f"load models/{cif2_name}, {name2}-{name1} \n")

    defattr = io.StringIO()
    if first_structure_id == 1:
        defattr.write("attribute: distance\nrecipient: residues\n")
    

    #write out the defattr files, basically just assign the disance as an attribute to each residue
    for id_ref, id_tgt, dist in zip(ids_ref, ids_tgt, distances):
        chain1, resnum1 = id_ref.split("-")
        chain2, resnum2 = id_tgt.split("-")

        pml.write(f'set_atom_property distance, {dist}, ({name1}-{name2} and chain {chain1} and resi {resnum1}), proptype=3 \n')
        pml.write(f'set_atom_property distance, {dist}, ({name2}-{name1} and chain {chain2} and resi {resnum2}), proptype=3 \n')


        defattr.write(f"\t#{first_structure_id}/{chain1}:{resnum1}\t{dist}\n")
        defattr.write(f"\t#{first_structure_id + 1}/{chain2}:{resnum2}\t{dist}\n")

    cxc = io.StringIO()
    #write cxc to open up models and the def attr files
    cxc.write(f"open models/{cif1_name} name {name1}-{name2} \n")
    cxc.write(f"open models/{cif2_name} name {name2}-{name1} \n")
   
    return defattr.getvalue(), cxc.getvalue(), pml.getvalue()

def chimera_movement_vectors_from_csv(df, output_path=None, fidelity=1, cmap=None, norm=None):
    """
    Export movement vectors to ChimeraX BILD format, with optional fidelity subsampling.
    """
    df = df.iloc[::fidelity]
    bild_string = generate_bild_string(df, cmap, norm)

    if output_path:
        with open(output_path, 'w') as f:
            f.write(bild_string)

    return bild_string

def generate_arrow_dicts(movement_dfs, cmap, vmin, vmax, fidelity=1):
    bild_output_dict = {}
    pml_output_dict = {}
    norm = mcolors.Normalize(vmin, vmax)

    for name, df in movement_dfs.items():
        # Group by reference chain (ChainID_Resnum1)
        if 'ChainID_Resnum1' in df.columns:
            df_grouped = df.groupby(df['ChainID_Resnum1'].apply(lambda x: str(x).split('-')[0]))
            for ref_chain, group in df_grouped:
                bild_output_dict[f"{name}_chain_{ref_chain}.bild"] = generate_bild_string(group, cmap, norm, fidelity=fidelity)
        else:
            bild_output_dict[name + ".bild"] = generate_bild_string(df, cmap, norm, fidelity=fidelity)
        pml_output_dict[name + ".pml"] = generate_pml_arrows(df.iloc[::fidelity], cmap, norm)
    return bild_output_dict, pml_output_dict


def save_bild_files_and_generate_chimerax_script(bild_output_dict, script_name="open_all_bilds.cxc"):
    """
    Generate a ChimeraX script to open all BILD files. Returns the script content as a string.
    """
    script_lines = []
    for filename in bild_output_dict.keys():
        script_lines.append(f"open {filename}")
    return '\n'.join(script_lines) + '\n'

def generate_bild_string(df, cmap=None, norm=None, fidelity=1):
    """
    Requires a matplot coloring map which is normalized to 0,1
    """
    df = df.iloc[::fidelity]
    coords1 = df['Coord1']
    coords2 = df['Coord2']
    distances = df['Distance'].apply(safe_eval)

    bild = io.StringIO()
    bild.write('.translate 0.0 0.0 0.0 \n')
    bild.write('.scale 1 \n')

    for coord1, coord2, dist in zip(coords1, coords2, distances):
        x1, y1, z1 = coord1
        x2, y2, z2 = coord2
        r, g, b = [int(255*c) for c in cmap(norm(dist))[:3]] if cmap and norm else (255, 255, 255)
        color_hex = f"{r:02x}{g:02x}{b:02x}" 
        bild.write(f'.color #{color_hex} \n')
        # Dynamically scale cone head length
        if dist is not None and dist > 0:
            head_length = min(2, max(0.1, 0.5 * math.log(dist + 1)))
        else:
            head_length = 0.1
        bild.write(f'.arrow {x1} {y1} {z1} {x2} {y2} {z2} {0.3} {head_length} \n')

    return bild.getvalue()


def generate_pml_arrows(df, cmap, norm, arrow_radius=0.2, arrow_head_ratio=0.2):
    """
    Creates pymol arrow visualization from a movement datatable.
    arrow_head_ratio specifies the percentage of the arrow taken up by the arrowhead.
    """

    coords1 = df['Coord1'].apply(safe_eval)
    coords2 = df['Coord2'].apply(safe_eval)
    distances = df['Distance'].apply(safe_eval)
    
    pml = io.StringIO()
    pml.write("from pymol.cgo import *\n")
    pml.write("from pymol import cmd\n\n")
    pml.write("all_arrows = []\n")

    for coord1, coord2, dist in zip(coords1, coords2, distances):
        r, g, b = cmap(norm(dist))[:3]
        
        vec = np.array(coord2) - np.array(coord1)
        tip = np.array(coord2) - arrow_head_ratio * vec
        
        x1, y1, z1 = coord1
        x2, y2, z2 = coord2

        cylinder = f"CYLINDER, {x1}, {y1}, {z1}, {tip[0]}, {tip[1]}, {tip[2]}, {arrow_radius}, {r}, {g}, {b}, {r}, {g}, {b}"
        cone = f"CONE, {tip[0]}, {tip[1]}, {tip[2]}, {x2}, {y2}, {z2}, {arrow_radius*4.5}, 0.0, {r}, {g}, {b}, {r}, {g}, {b}, 1.0, 0.0"
        pml.write(f"all_arrows += [{cylinder}, {cone}]\n")

    pml.write("cmd.load_cgo(all_arrows, 'arrows')\n")
    return pml.getvalue()
