import pandas as pd
import numpy as np
import ast

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

def generate_chimerax_script(input_csv, output_script, chains=None, thresholds=None):
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
                out.write(f"show /{chain}\n")