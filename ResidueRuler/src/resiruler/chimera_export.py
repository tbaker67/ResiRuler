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

def generate_multiple_movement_scripts(movement_dfs, ref_name):

    vmin, vmax = float('inf'), float('-inf')
    for df in movement_dfs.values():
        vmin = min(vmin, df['Distance'].min())
        vmax = max(vmax, df['Distance'].max())

    #either going to need a counter or way to identify by name
    full_def_attr = StringIO()
    full_cxc_script = StringIO()

    ids = 1
    for tgt_structure_name, movement_df in movement_dfs.items():
        #Assumes that the names are simply the cif files without the extension
        ref_cif_name = ref_name + ".cif"
        tgt_cif_name = tgt_structure_name + ".cif"
        def_attr, cxc_script = generate_cxc_scripts(movement_df, ref_cif_name, ref_name, tgt_cif_name, tgt_structure_name, ids)
        full_def_attr.write(def_attr + '\n')
        full_cxc_script.write(cxc_script + '\n')

        ids += 2

    full_cxc_script.write("open full_defattr.defattr \n")

    full_cxc_script.write(f"color #{1}-{ids - 1} grey")
    
    full_cxc_script.write(f"color byattribute r:distance #{1}-{ids - 1} target scab palette 0,#00008B:{vmax / 5:.2f},#20073a:{(2 * vmax  / 5):.2f},#6d1950:{(3 * vmax/ 5):.2f},#bd4545:{ (4 * vmax / 5):.2f},#d48849:{vmax:.2f},#f0d171\n")

    #color bar/legend
    full_cxc_script.write(f"key #00008B:0 #20073a:{vmax / 5:.2f} #6d1950:{(2 * vmax  / 5):.2f} #bd4545:{(3 * vmax/ 5):.2f} #d48849:{ (4 * vmax / 5):.2f} #f0d171:{vmax:.2f}\n")

    return full_def_attr.getvalue(), full_cxc_script.getvalue()

def generate_cxc_scripts(df, cif1_name, structure_name1, cif2_name, structure_name2, first_structure_id):
    """
    Generate defattr files, a bild file, and a cxc chimera script to color models corresponding to distance between corresponding residues in the reference and target structures
    """
    distances = df['Distance'].apply(safe_eval)
    ids_ref = df['ChainID_Resnum1']
    ids_tgt = df['ChainID_Resnum2']

    #set up color scale
    vmin, vmax = distances.min(), distances.max()
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.cm.get_cmap('plasma')


    name1 = structure_name1
    name2 = structure_name2

    defattr = StringIO()
    
    if first_structure_id == 1:
        defattr.write("attribute: distance\nrecipient: residues\n")
    

    #write out the defattr files, basically just assign the disance as an attribute to each residue
    for id_ref, id_tgt, dist in zip(ids_ref, ids_tgt, distances):
        chain1, resnum1 = id_ref.split("-")
        chain2, resnum2 = id_tgt.split("-")

        defattr.write(f"\t#{first_structure_id}/{chain1}:{resnum1}\t{dist}\n")
        defattr.write(f"\t#{first_structure_id + 1}/{chain2}:{resnum2}\t{dist}\n")

    cxc = StringIO()
    #write cxc to open up models and the def attr files
    cxc.write(f"open models/{cif1_name} name {name1}->{name2} \n")
    #cxc.write(f"name {name1} #{first_structure_id} \n")
    cxc.write(f"open models/{cif2_name} name {name2}->{name1} \n")
    #cxc.write(f"name {name2} #{first_structure_id + 1} \n ")

    
    #actually color the residues
   
    return defattr.getvalue(), cxc.getvalue()

def chimera_movement_vectors_from_csv(df, output_path=None,):
    """
    
    """
    bild_string = generate_bild_string(df)

    if output_path:
        with open(output_path, 'w') as f:
            f.write(bild_string)

    return bild_string

def generate_multiple_bilds(movement_dfs):
    output_dict = {}
    vmin, vmax = float('inf'), float('-inf')
    for df in movement_dfs.values():
        vmin = min(vmin, df['Distance'].min())
        vmax = max(vmax, df['Distance'].max())
    for name,df in movement_dfs.items():
        coords1 = df['Coord1']
        coords2 = df['Coord2']
        distances = df['Distance'].apply(safe_eval)

        vmin, vmax = distances.min(), distances.max()
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
        cmap = plt.cm.get_cmap('plasma')

        output = io.StringIO()
        output.write('.translate 0.0 0.0 0.0 \n')
        output.write('.scale 1 \n')

        for coord1, coord2, dist in zip(coords1, coords2, distances):
            x1, y1, z1 = coord1
            x2, y2, z2 = coord2
            r, g, b = [int(255 * c) for c in cmap(norm(dist))[:3]]
            color_hex = f"{r:02x}{g:02x}{b:02x}"
            output.write(f'.color #{color_hex} \n')
            output.write(f'.arrow {x1} {y1} {z1} {x2} {y2} {z2} \n')

        output_dict[name + ".bild"] = output.getvalue()

    return output_dict

def generate_bild_string(df):
    coords1 = df['Coord1']
    coords2 = df['Coord2']
    distances = df['Distance'].apply(safe_eval)

    vmin, vmax = distances.min(), distances.max()
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.cm.get_cmap('plasma')

    output = io.StringIO()
    output.write('.translate 0.0 0.0 0.0 \n')
    output.write('.scale 1 \n')

    for coord1, coord2, dist in zip(coords1, coords2, distances):
        x1, y1, z1 = coord1
        x2, y2, z2 = coord2
        r, g, b = [int(255 * c) for c in cmap(norm(dist))[:3]]
        color_hex = f"{r:02x}{g:02x}{b:02x}"
        output.write(f'.color #{color_hex} \n')
        output.write(f'.arrow {x1} {y1} {z1} {x2} {y2} {z2} \n')

    return output.getvalue()