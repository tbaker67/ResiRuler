import py3Dmol
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from collections import defaultdict


from src.resiruler.chimera_export import get_color, safe_eval

def start_pymol_viewer(cif_file):
    if cif_file is not None:
        cif_str = cif_file.read().decode("utf-8")  # decode BytesIO to str
        view = py3Dmol.view(width=800, height=800)
        view.addModel(cif_str, 'cif')
        view.setStyle({'cartoon': {'color': 'spectrum'}})
        view.zoomTo()
        return view

def draw_links_pymol(df, view, thresholds=None):
    df = df.dropna()
    chains1 = df['Chain1_Residue1'].str.split('_').str[0]
    chains2 = df['Chain2_Residue2'].str.split('_').str[0]
    all_chains = pd.concat([chains1, chains2]).unique()

    view.setStyle({}, {}) 
    for chain_id in all_chains:
        view.setStyle({'chain': chain_id}, {'cartoon': {'colorscheme': 'chainid'}})

    starts = df['Coord1']
    ends = df['Coord2']
    distances = float(df['Distance'])
    for  start,end,dist in zip(starts,ends,distances):
        view.addCylinder({
            'start': {'x': float(start[0]), 'y': float(start[1]), 'z': float(start[2])},
            'end':   {'x': float(end[0]), 'y': float(end[1]), 'z': float(end[2])},
            'radius': 0.5,
            'fromCap': 1,
            'toCap': 1,
            'color': get_color(dist, thresholds),
            'opacity': 1.0
        })

    view.zoomTo()
    return view

def draw_movement_shift_pymol(df, view):
    df = df.dropna(subset=['ChainID_Resnum1', 'Distance'])

    vmin, vmax = df['Distance'].min(), df['Distance'].max()
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.cm.get_cmap('plasma')

    view.setStyle({}, {})  # reset styles

    for _, row in df.iterrows():
        # Parse chain and residue number from 'C_58'
        chain, resi = row['ChainID_Resnum1'].split('_')

        dist = row['Distance']
        r, g, b = [int(255 * c) for c in cmap(norm(dist))[:3]]
        color_str = f'rgb({r},{g},{b})'

        # Apply cartoon color to the residue
        view.setStyle({'chain': chain, 'resi': str(resi)}, {
            'cartoon': {'color': color_str}
        })

    view.zoomTo()
    return view


def draw_movement_vectors_py3dmol(df, view, radius=0.3, head_radius=0.5):
    """
    Draws colored vectors and heads from Coord1 to Coord2, colored by Distance.
    """
    df = df.dropna(subset=['Coord1', 'Coord2', 'Distance'])

    coords1 = df['Coord1']
    coords2 = df['Coord2']
    distances = df['Distance'].astype(float)

    view.setStyle({}, {})  
    # Normalize distances to colormap
    vmin, vmax = distances.min(), distances.max()
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.cm.get_cmap('plasma')

    for start, end, dist in zip(coords1, coords2, distances):
        r, g, b = [int(255 * c) for c in cmap(norm(dist))[:3]]
        color_str = f'rgb({r},{g},{b})'

        # Draw line
        view.addLine({
            'start': {'x': start[0], 'y': start[1], 'z': start[2]},
            'end':   {'x': end[0],   'y': end[1],   'z': end[2]},
            'color': color_str,
            'radius': radius,
            'dashed': False
        })

        # Draw head sphere at end
        view.addSphere({
            'center': {'x': end[0], 'y': end[1], 'z': end[2]},
            'color': color_str,
            'radius': head_radius
        })

    view.zoomTo()
    return view
        
        
    

