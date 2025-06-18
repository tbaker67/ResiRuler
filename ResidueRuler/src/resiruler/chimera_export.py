import pandas as pd
import numpy as np
import ast
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import io

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
    df["Coord1"] = df["Coord1"].apply(safe_eval)
    df["Coord2"] = df["Coord2"].apply(safe_eval)
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

def chimera_color_shift_from_csv(df, output_script=None, chain_mapping=None):
    script = generate_chimera_color_script(df, chain_mapping=chain_mapping)

    if output_script:
        with open(output_script, 'w') as f:
            f.write(script)

    return script

def generate_chimera_color_script(df, chain_mapping=None):
    coords = df['Coord1'].apply(safe_eval)
    distances = df['Distance'].apply(safe_eval)
    ids = df['ChainID_Resnum1']

    vmin, vmax = distances.min(), distances.max()
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.cm.get_cmap('plasma')

    output = io.StringIO()
    output.write('hide all \n')

    for coord, dist, id in zip(coords, distances, ids):
        x, y, z = coord
        r, g, b = [int(255 * c) for c in cmap(norm(dist))[:3]]
        color_hex = f"{r:02x}{g:02x}{b:02x}"

        chain, resnum = id.split("_")

        output.write(f"color /{chain}:{resnum} #{color_hex} \n")
        output.write(f"show /{chain}:{resnum} \n")

        if chain_mapping and chain in chain_mapping:
            mapped_chain = chain_mapping[chain]
            output.write(f"color /{mapped_chain}:{resnum} #{color_hex} \n")
            output.write(f"show /{mapped_chain}:{resnum} \n")

    return output.getvalue()

def chimera_movement_vectors_from_csv(df, output_path=None,):
    bild_string = generate_bild_string(df)

    if output_path:
        with open(output_path, 'w') as f:
            f.write(bild_string)

    return bild_string

def generate_bild_string(df):
    import matplotlib.colors as mcolors
    import matplotlib.pyplot as plt

    coords1 = df['Coord1'].apply(safe_eval)
    coords2 = df['Coord2'].apply(safe_eval)
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