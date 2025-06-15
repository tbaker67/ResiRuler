import pandas as pd
import numpy as np
import ast
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

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

def draw_links(input_csv, output_script, chains=None, thresholds=None):
    df = pd.read_csv(input_csv)
    df["Coord1"] = df["Coord1"].apply(safe_eval)
    df["Coord2"] = df["Coord2"].apply(safe_eval)
    df["Distance"] = df["Distance"].apply(safe_eval)

    with open(output_script, 'w') as out:
        for _, row in df.iterrows():
            coord1, coord2, dist = row['Coord1'], row['Coord2'], row['Distance']
            if any(pd.isna([coord1, coord2, dist])):
                continue
            color = get_color(dist, thresholds)
            out.write(f"shape cylinder radius 1 fromPoint {coord1[0]},{coord1[1]},{coord1[2]} "
                      f"toPoint {coord2[0]},{coord2[1]},{coord2[2]} color {color}\n")

        if chains:
            out.write("hide\n")
            for chain in chains:
                out.write(f"cartoon /{chain}\n")


def chimera_color_shift_from_csv(csv_path, output_script,chain_mapping=None):
    df = pd.read_csv(csv_path)

    # Parse coordinates
    coords = df['Coord1'].apply(safe_eval)
    distances = df['Distance'].apply(safe_eval)
    ids = df['ChainID_Resnum'].apply(safe_eval)

    # Normalize distances
    vmin, vmax = distances.min(), distances.max()
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.cm.get_cmap('plasma')  

    # Prepare ChimeraX script
    with open(output_script, 'w') as f:

        f.write('hide all \n')

        for (coord, dist, id) in zip(coords, distances, ids):
            x, y, z = coord
            color_rgba = cmap(norm(dist))  # Returns (r,g,b,a)
            r, g, b = [int(255 * c) for c in color_rgba[:3]]
            split = id.split("_")
            chain = split[0]
            resnum = split[1]
            color_hex = f"{r:02x}{g:02x}{b:02x}"
            f.write(f"color /{chain}:{resnum} #{color_hex} \n")
            f.write(f"show /{chain}:{resnum} \n")
            if chain_mapping and chain in chain_mapping:
                f.write(f"color /{chain_mapping[chain]}:{resnum} #{color_hex} \n")
                f.write(f"show /{chain_mapping[chain]}:{resnum} \n")

def chimera_movement_vectors_from_csv(csv_path, output_bild):
    df = pd.read_csv(csv_path)
    coords1 = df['Coord1'].apply(safe_eval)
    coords2 = df['Coord2'].apply(safe_eval)
    distances = df['Distance'].apply(safe_eval)

    #Normalize Distances
    vmin, vmax = distances.min(), distances.max()
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.cm.get_cmap('plasma')  

    with open(output_bild, 'w') as f:

        f.write('.translate 0.0 0.0 0.0 \n')
        f.write('.scale 1 \n')

        for (coord1, coord2, dist) in zip(coords1, coords2, distances):
            print(coord1)
            x1,y1,z1 = coord1[0],coord1[1],coord1[2]
            x2,y2,z2 = coord2[0],coord2[1],coord2[2]
            color_rgba = cmap(norm(dist))
            r, g, b = [int(255 * c) for c in color_rgba[:3]]
            color_hex = f"{r:02x}{g:02x}{b:02x}"

            f.write(f'.color #{color_hex} \n')
            f.write(f'.arrow {x1} {y1} {z1} {x2} {y2} {z2} \n')




