import py3Dmol
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from src.resiruler.chimera_export import get_color_discrete, get_color_gradient




def start_pymol_viewer(cif_file):
    if cif_file is not None:
        cif_str = cif_file.read().decode("utf-8")  # decode BytesIO to str
        view = py3Dmol.view(width=800, height=800)
        view.addModel(cif_str, 'cif')
        view.setStyle({'cartoon': {'color': 'spectrum'}})
        view.zoomTo()
        return view
    
def get_link_color_pymol(distance, color_mode, thresholds=None, cmap=None, min_val=None, max_val=None):
    if color_mode == "discrete":
        return get_color_discrete(distance, thresholds)
    else:
        return get_color_gradient(distance, cmap,min_val, max_val)

def draw_links_pymol(df, view, color_mode, **kwargs):
    """
    Draws a pymol viewer in the UI to display the a preview for the link visualization
    **kwargs: Extra arguments passed to get_color to work with the appropriate type i.e,
                  color_mode='gradient', cmap=<colormap>, min_val=0, max_val=10
    """
    df = df.dropna()
    chains1 = df['ChainID_Resnum1'].str.split('-').str[0]
    chains2 = df['ChainID_Resnum2'].str.split('-').str[0]
    all_chains = pd.concat([chains1, chains2]).unique()

    
    view.setStyle({}, {}) 

    
    for chain_id in all_chains:
        view.setStyle({'chain': chain_id}, {'cartoon': {'colorscheme': 'chainid'}})

    
    for _, row in df.iterrows():
        coord1, coord2, dist = row['Coord1'], row['Coord2'], row['Distance']
        if any(pd.isna([coord1, coord2, dist])):
            continue

        color = get_link_color_pymol(dist, color_mode, **kwargs)
        view.addCylinder({
            'start': {'x': float(coord1[0]), 'y': float(coord1[1]), 'z': float(coord1[2])},
            'end':   {'x': float(coord2[0]), 'y': float(coord2[1]), 'z': float(coord2[2])},
            'radius': 0.5,
            'fromCap': 1,
            'toCap': 1,
            'color': color,
            'opacity': 1.0
        })

    view.zoomTo()
    return view

def draw_movement_shift_pymol(df, view, cmap=plt.cm.get_cmap('plasma')):
    df = df.dropna(subset=['ChainID_Resnum1', 'Distance'])

    vmin, vmax = df['Distance'].min(), df['Distance'].max()
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    

    #view.setStyle({}, {})  # reset styles

    for _, row in df.iterrows():
        # Parse chain and residue number from 'C-58'
        chain, resi = row['ChainID_Resnum1'].split('-')

        dist = row['Distance']
        r, g, b = [int(255 * c) for c in cmap(norm(dist))[:3]]
        color_str = f'rgb({r},{g},{b})'

        # Apply cartoon color to the residue
        view.setStyle({'chain': chain, 'resi': str(resi[1])}, {
            'cartoon': {'color': color_str}
        })

    view.zoomTo()
    return view


def draw_movement_vectors_py3dmol(df, view, radius=0.3, head_radius=0.5, cmap = plt.cm.get_cmap('plasma')):
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
    

    for start, end, dist in zip(coords1, coords2, distances):
        r, g, b = [int(255 * c) for c in cmap(norm(dist))[:3]]
        color_str = f'rgb({r},{g},{b})'

        # Draw line
        view.addLine({
            'start': {'x': float(start[0]), 'y': float(start[1]), 'z': float(start[2])},
            'end':   {'x': float(end[0]),   'y': float(end[1]),   'z': float(end[2])},
            'color': color_str,
            #'radius': radius,
            'dashed': False
        })

        # Draw head sphere at end
        view.addSphere({
            'center': {'x': float(end[0]), 'y': float(end[1]), 'z': float(end[2])},
            'color': color_str,
            'radius': head_radius
        })

    view.zoomTo()
    return view
        
        
    

