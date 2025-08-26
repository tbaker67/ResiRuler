import pandas as pd
import numpy as np
import ast
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import io
from pathlib import Path
from io import StringIO
from os.path import basename

def safe_eval(val):
    """
    Extract the values from a dataframe safely and handle possibility of them being strings
    """
    if pd.isna(val):
        return np.nan
    #Already good to use
    if isinstance(val, (int, float, list, tuple, dict)):
        return val  
    
    #If a string
    if isinstance(val, str):
        try:
            return ast.literal_eval(val)
        except (ValueError, SyntaxError):
            return np.nan

    return np.nan  # fallback for unrecognized types

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

def generate_chimera_palette_string(palette, vmin, vmax):

    palette_string = f"{vmin}, {palette[0]}"

    #Evenly space the palette colors with vmin and vmax at the extremes
    #This also ensures colors match across real-valued and normalized distances as long as they are normalized along (vmin, vmax)
    for i in range(1,len(palette)):
        val = vmin + i * (vmax - vmin) / (len(palette) - 1)
        palette_string += f":{val, palette[i]}"
    
    return palette_string

def generate_pml_palette_string(palette):
    palette_string = ""
    for hex_color in palette:
        #regular hex codes have #00008B, for pymol it needs to be in the form 0x00008B
        pymol_color = hex_color.replace("#", "0x")
        palette_string += f"{pymol_color} "
    return palette_string


def generate_multiple_movement_scripts(movement_dfs, ref_name, palette, vmin, vmax):

    full_def_attr = StringIO()
    full_cxc_script = StringIO()
    full_pml_script = StringIO()

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

    chimera_palette_string = generate_chimera_palette_string(palette, vmin, vmax)
    full_cxc_script.write("open full_defattr.defattr \n")
    full_cxc_script.write(f"color #{1}-{ids - 1} grey \n")
    full_cxc_script.write(f"color byattribute r:distance #{1}-{ids - 1} target scab palette {chimera_palette_string}\n")
    full_cxc_script.write(f"key {chimera_palette_string}\n")

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

    pml = StringIO()
    pml.write(f"load models/{cif1_name}, {name1}-{name2} \n")
    pml.write(f"load models/{cif2_name}, {name2}-{name1} \n")

    defattr = StringIO()
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

    cxc = StringIO()
    #write cxc to open up models and the def attr files
    cxc.write(f"open models/{cif1_name} name {name1}-{name2} \n")
    #cxc.write(f"name {name1} #{first_structure_id} \n")
    cxc.write(f"open models/{cif2_name} name {name2}-{name1} \n")
    #cxc.write(f"name {name2} #{first_structure_id + 1} \n ")
   
    return defattr.getvalue(), cxc.getvalue(), pml.getvalue()

def chimera_movement_vectors_from_csv(df, output_path=None,):
    """
    
    """
    bild_string = generate_bild_string(df)

    if output_path:
        with open(output_path, 'w') as f:
            f.write(bild_string)

    return bild_string

def generate_arrow_dicts(movement_dfs, cmap, vmin, vmax):
    bild_output_dict = {}
    pml_output_dict = {}
    norm = mcolors.Normalize(vmin,vmax)

    for name,df in movement_dfs.items():
        bild_output_dict[name + ".bild"] = generate_bild_string(df, cmap, norm)
        pml_output_dict[name + ".pml"] = generate_pml_arrows(df, cmap, norm)
        
    return bild_output_dict, pml_output_dict

def generate_bild_string(df, cmap, norm):
    """
    Requires a matplot coloring map which is normalized to 0,1
    """
    coords1 = df['Coord1']
    coords2 = df['Coord2']
    distances = df['Distance'].apply(safe_eval)

    bild = StringIO()
    bild.write('.translate 0.0 0.0 0.0 \n')
    bild.write('.scale 1 \n')

    for coord1, coord2, dist in zip(coords1, coords2, distances):
        x1, y1, z1 = coord1
        x2, y2, z2 = coord2
        
        r, g, b = [int(255*c) for c in cmap(norm(dist))[:3]] # r,g,b on range [0,255]
        color_hex = f"{r:02x}{g:02x}{b:02x}" 
        bild.write(f'.color #{color_hex} \n')
        bild.write(f'.arrow {x1} {y1} {z1} {x2} {y2} {z2} \n')

    return bild.getvalue()


def generate_pml_arrows(df, cmap, norm, arrow_radius=0.2, arrow_head_ratio=0.2):

    """
    Creates pymol arrow visualization from a movement datatable
    arrow_head_ratio specifies the percentage of the array that is taken up by the 
    """

    coords1 = df['Coord1']
    coords2 = df['Coord2']
    distances = df['Distance'].apply(safe_eval)
    diff_vecs = df['Diff_Vec']
    
    pml = StringIO()
    pml.write("from pymol.cgo import *\n")
    pml.write("from pymol import cmd\n\n")
    pml.write("all_arrows = []\n")

    for coord1, coord2, dist, diff_vec in zip(coords1, coords2, distances, diff_vecs):
        r, g, b = cmap(norm(dist))[:3]
        vec = np.array(diff_vec)
        tip = np.array(coord2) - arrow_head_ratio*vec
        x1, y1, z1 = coord1
        x2, y2, z2 = coord2
        pml.write("all_arrows += [")
        pml.write(f"CYLINDER, {x1}, {y1}, {z1}, {tip[0]}, {tip[1]}, {tip[2]}, {arrow_radius}, {r}, {g}, {b}, {r}, {g}, {b}, ")
        pml.write(f"CONE, {tip[0]}, {tip[1]}, {tip[2]}, {x2}, {y2}, {z2}, {arrow_radius*1.5}, 0.0, {r}, {g}, {b}, {r}, {g}, {b}, 1.0, 0.0")
        pml.write("]\n\n")

    pml.write("cmd.load_cgo(all_arrows, 'arrows')\n")
    return pml.getvalue()
