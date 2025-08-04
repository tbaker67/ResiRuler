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

def get_color(distance, thresholds):
    if distance <= thresholds["green"]:
        return "green"
    elif distance <= thresholds["yellow"]:
        return "yellow"
    else:
        return "red"

def draw_links(df, output_script=None, chains=None, thresholds=None):
    script = generate_chimera_link_script(df, chains=chains, thresholds=thresholds)

    if output_script:
        with open(output_script, 'w') as f:
            f.write(script)

    return script


def generate_chimera_link_script(df, chains=None, thresholds=None):
    """
    Generate a chimera cxc script to draw in links based on a desired set of thresholds
    """
    df["Coord1"] = df["Coord1"]
    df["Coord2"] = df["Coord2"]
    df["Distance"] = df["Distance"].apply(safe_eval)

    output = io.StringIO()

    for _, row in df.iterrows():
        coord1, coord2, dist = row['Coord1'], row['Coord2'], row['Distance']
        if any(pd.isna([coord1, coord2, dist])):
            continue
        color = get_color(dist, thresholds)
        output.write(f"shape cylinder radius 1 fromPoint {coord1[0]},{coord1[1]},{coord1[2]} "
                     f"toPoint {coord2[0]},{coord2[1]},{coord2[2]} color {color}\n")

    if chains:
        output.write("hide\n")
        for chain in chains:
            output.write(f"cartoon /{chain}\n")

    return output.getvalue()

def generate_cxc_scripts(df, cif1_name, cif2_name, structure_name1, structure_name2, chain_mapping=None):
    """
    Generate defattr files, a bild file, and a cxc chimera script to color models corresponding to distance between corresponding residues in the reference and target structures
    """
    distances = df['Distance'].apply(safe_eval)
    ids = df['ChainID_Resnum1']

    #set up color scale
    vmin, vmax = distances.min(), distances.max()
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.cm.get_cmap('plasma')


    name1 = structure_name1
    name2 = structure_name2

    defattr1 = StringIO()
    defattr2 = StringIO()

    defattr1.write("attribute: distance\nrecipient: residues\n")
    defattr2.write("attribute: distance\nrecipient: residues\n")

    #write out the defattr files, basically just assign the disance as an attribute to each residue
    for id_, dist in zip(ids, distances):
        chain, resnum = id_.split("_")

        defattr1.write(f"\t#1/{chain}:{resnum}\t{dist}\n")
        if chain_mapping and chain in chain_mapping:
            mapped_chain = chain_mapping[chain]
            defattr2.write(f"\t#2/{mapped_chain}:{resnum}\t{dist}\n")

    cxc = StringIO()
    #write cxc to open up models and the def attr files
    cxc.write(f"open {cif1_name} name {name1}\n")
    cxc.write(f"open {cif2_name} name {name2}\n")
    cxc.write(f"open {name1}_colors.defattr\n")
    cxc.write(f"open {name2}_colors.defattr\n")
    
    #actually color the residues
    cxc.write(f"color byattribute r:distance #1-2 target scab palette 0,#00008B:{vmax / 5:.2f},#20073a:{(2 * vmax  / 5):.2f},#6d1950:{(3 * vmax/ 5):.2f},#bd4545:{ (4 * vmax / 5):.2f},#d48849:{vmax:.2f},#f0d171\n")

    #color bar/legend
    cxc.write(f"key #00008B:0 #20073a:{vmax / 5:.2f} #6d1950:{(2 * vmax  / 5):.2f} #bd4545:{(3 * vmax/ 5):.2f} #d48849:{ (4 * vmax / 5):.2f} #f0d171:{vmax:.2f}\n")

    return defattr1.getvalue(), defattr2.getvalue(), cxc.getvalue()

def chimera_movement_vectors_from_csv(df, output_path=None,):
    """
    
    """
    bild_string = generate_bild_string(df)

    if output_path:
        with open(output_path, 'w') as f:
            f.write(bild_string)

    return bild_string

def generate_bild_string(df):
    import matplotlib.colors as mcolors
    import matplotlib.pyplot as plt

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