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

def generate_chimera_coloring_palette_string(palette, positions):

    palette_string = ""
    for color, position in zip(palette, positions):
        palette_string += f"{position:.2f},{color}:"
    
    return palette_string[:-1]

def generate_chimera_key_string(palette, positions):
    return "key " + " ".join(f"{color}:{int(round(position))}" for color, position in zip(palette, positions))

def generate_pml_palette_string(palette):
    
    pymol_colors = [hex_color.replace("#", "0x") for hex_color in palette]
    return "_".join(pymol_colors)

def generate_shift_scripts(df, ref_loaded_name, cif2_name, structure_name2, first_structure_id, write_ref_load=False, ref_cif_name=None):
    """
    Generate defattr files, a cxc chimera script and a pml script to color models
    corresponding to distance between corresponding residues in the reference and
    target structures.

    ref_loaded_name: the object name the reference was already loaded under (loaded once, elsewhere)
    write_ref_load: if True, also emit the `load` line for the reference (used only for the
                    first call, so the reference file is only read from disk once)
    """
    distances = df['Distance'].apply(safe_eval)
    ids_ref = df['ChainID_Resnum1']
    ids_tgt = df['ChainID_Resnum2']

    name1 = ref_loaded_name
    name2 = structure_name2

    disp_name1 = f"{name1}_disp_{name2}"
    disp_name2 = f"{name2}_disp_{name1}"

    pml = io.StringIO()
    if write_ref_load:
        pml.write(f"load models/{ref_cif_name}, {name1} \n")
    pml.write(f"load models/{cif2_name}, {name2} \n")

    pml.write(f"create {disp_name1}, {name1} \n")
    pml.write(f"create {disp_name2}, {name2} \n")

    defattr = io.StringIO()
    if first_structure_id == 1:
        defattr.write("attribute: distance\nrecipient: residues\n")

    dist_map1 = {}
    dist_map2 = {}

    for id_ref, id_tgt, dist in zip(ids_ref, ids_tgt, distances):
        if pd.isna(dist):
            continue

        chain1, resnum1 = (s.strip() for s in id_ref.split("-"))
        chain2, resnum2 = (s.strip() for s in id_tgt.split("-"))

        dist_map1[(chain1, resnum1)] = dist
        dist_map2[(chain2, resnum2)] = dist

        defattr.write(f"\t#{disp_name1}/{chain1}:{resnum1}\t{dist}\n")
        defattr.write(f"\t#{disp_name2}/{chain2}:{resnum2}\t{dist}\n")

    # single iterate pass per object instead of one alter per residue
    pml.write("python\n")
    pml.write(f"_dist_map1 = {dist_map1!r}\n")
    pml.write(f"cmd.alter('{disp_name1}', 'b = _dist_map1.get((chain, resi), 0.0)', space={{'_dist_map1': _dist_map1}})\n")
    pml.write(f"_dist_map2 = {dist_map2!r}\n")
    pml.write(f"cmd.alter('{disp_name2}', 'b = _dist_map2.get((chain, resi), 0.0)', space={{'_dist_map2': _dist_map2}})\n")
    pml.write("python end\n")

    cxc = io.StringIO()
    if write_ref_load:
        cxc.write(f"open models/{ref_cif_name} name {name1} \n")
    cxc.write(f"open models/{cif2_name} name {name2} \n")

    return defattr.getvalue(), cxc.getvalue(), pml.getvalue()


def generate_multiple_displacement_scripts(displacement_dfs, ref_name, palette, positions):

    full_def_attr = io.StringIO()
    full_cxc_script = io.StringIO()
    full_pml_script = io.StringIO()

    ref_cif_name = ref_name + ".cif"

    ids = 1
    disp_names = []
    is_first = True
    for tgt_structure_name, displacement_df in displacement_dfs.items():
        tgt_cif_name = tgt_structure_name + ".cif"
        def_attr, cxc_script, pml_script = generate_shift_scripts(
            displacement_df,
            ref_name,
            tgt_cif_name,
            tgt_structure_name,
            ids,
            write_ref_load=is_first,   # only load reference from disk on the first pass
            ref_cif_name=ref_cif_name,
        )

        full_def_attr.write(def_attr + '\n')
        full_cxc_script.write(cxc_script + '\n')
        full_pml_script.write(pml_script + '\n')

        disp_names.append(f"{ref_name}_disp_{tgt_structure_name}")
        disp_names.append(f"{tgt_structure_name}_disp_{ref_name}")

        ids += 2
        is_first = False


    chimera_coloring_palette_string = generate_chimera_coloring_palette_string(palette, positions)
    full_cxc_script.write("open full_defattr.defattr \n")
    full_cxc_script.write(f"color #{1}-{ids - 1} grey \n")
    full_cxc_script.write(f"color byattribute r:distance #{1}-{ids - 1} target scab palette {chimera_coloring_palette_string}\n")
    chimera_key_string = generate_chimera_key_string(palette, positions)
    full_cxc_script.write(f"{chimera_key_string}\n")

    pml_palette_string = generate_pml_palette_string(palette)
    pml_selection = " or ".join(disp_names)
    min_val = min(positions)
    max_val = max(positions)
    full_pml_script.write(f'spectrum b, {pml_palette_string}, {pml_selection}, minimum={min_val}, maximum={max_val}\n')

    return full_def_attr.getvalue(), full_cxc_script.getvalue(), full_pml_script.getvalue()

def chimera_displacement_vectors_from_csv(df, output_path=None, fidelity=1, cmap=None, norm=None):
    """
    Export displacement vectors to ChimeraX BILD format, with optional fidelity subsampling.
    """
    df = df.iloc[::fidelity]
    bild_string = generate_bild_string(df, cmap, norm)

    if output_path:
        with open(output_path, 'w') as f:
            f.write(bild_string)

    return bild_string


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


def generate_pml_arrows(df, cmap, norm, arrow_radius=0.2, arrow_head_ratio=0.2, object_name="arrows", group_name=None):
    """
    Creates pymol arrow visualization from a displacement datatable.
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

    pml.write(f"cmd.load_cgo(all_arrows, '{object_name}')\n")
    if group_name:
        pml.write(f"group {group_name}, {object_name}\n")
    return pml.getvalue()


def generate_arrow_dicts(displacement_dfs, cmap, vmin, vmax, fidelity=1):
    bild_output_dict = {}
    pml_output_dict = {}
    norm = mcolors.Normalize(vmin, vmax)

    for name, df in displacement_dfs.items():
        group_name = f"arrows_{name}_group"
        if 'ChainID_Resnum1' in df.columns:
            df_grouped = df.groupby(df['ChainID_Resnum1'].apply(lambda x: str(x).split('-')[0].strip()))
            for ref_chain, group in df_grouped:
                object_name = f"arrows_{name}_chain_{ref_chain}"
                bild_output_dict[f"{name}_chain_{ref_chain}.bild"] = generate_bild_string(group, cmap, norm, fidelity=fidelity)
                pml_output_dict[f"{name}_chain_{ref_chain}.pml"] = generate_pml_arrows(
                    group.iloc[::fidelity], cmap, norm, object_name=object_name, group_name=group_name
                )
        else:
            object_name = f"arrows_{name}"
            bild_output_dict[name + ".bild"] = generate_bild_string(df, cmap, norm, fidelity=fidelity)
            pml_output_dict[name + ".pml"] = generate_pml_arrows(
                df.iloc[::fidelity], cmap, norm, object_name=object_name, group_name=group_name
            )

    return bild_output_dict, pml_output_dict


def save_pml_files_and_generate_pymol_script(pml_output_dict, script_name="run_all_arrows.pml"):
    """
    Generate a PyMOL script to run all per-chain arrow .pml files. Returns the script content as a string.
    """
    script_lines = []
    for filename in pml_output_dict.keys():
        script_lines.append(f"run {filename}")
    return '\n'.join(script_lines) + '\n'
