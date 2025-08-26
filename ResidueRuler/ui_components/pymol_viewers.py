import py3Dmol
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from src.resiruler.chimera_export import get_color_discrete, get_color_gradient
import numpy as np
import plotly.graph_objects as go




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

def plot_vectors_plotly(df, cmap, vmin,vmax, show_heads=True, line_width=2, head_size=3, width=1200, height=800):
    df = df.dropna(subset=['Coord1','Coord2','Distance'])
    if df.empty:
        return go.Figure()

    coords1 = np.array(df['Coord1'].tolist(), dtype=float)
    coords2 = np.array(df['Coord2'].tolist(), dtype=float)
    distances = df['Distance'].astype(float).values

    # Normalize for colormap
    normed = mcolors.Normalize(vmin, vmax)(distances)
    colors_rgb = (np.array([cmap(n)[:3] for n in normed]) * 255).astype(int)
    color_strings = [f'rgb({r},{g},{b})' for r,g,b in colors_rgb]

    # Prepare data with NaN breaks
    x = np.empty((len(coords1)*3,))
    y = np.empty_like(x)
    z = np.empty_like(x)
    line_colors = []

    for i, (start, end, color) in enumerate(zip(coords1, coords2, color_strings)):
        idx = i*3
        x[idx:idx+3] = [start[0], end[0], np.nan]
        y[idx:idx+3] = [start[1], end[1], np.nan]
        z[idx:idx+3] = [start[2], end[2], np.nan]
        line_colors.extend([color]*3)

    fig = go.Figure()

    # Add vectors
    fig.add_trace(go.Scatter3d(
        x=x, y=y, z=z,
        mode='lines',
        line=dict(color=line_colors, width=line_width),
        hoverinfo='skip'
    ))

    # Add heads
    if show_heads:
        fig.add_trace(go.Scatter3d(
            x=coords2[:,0], y=coords2[:,1], z=coords2[:,2],
            mode='markers',
            marker=dict(size=head_size, color=color_strings, symbol='circle'),
            hoverinfo='skip'
        ))

    # Layout: remove all gridlines, axes, and background
    fig.update_layout(scene=dict(
        xaxis=dict(showbackground=False, showticklabels=False, showgrid=False, zeroline=False, visible=False),
        yaxis=dict(showbackground=False, showticklabels=False, showgrid=False, zeroline=False, visible=False),
        zaxis=dict(showbackground=False, showticklabels=False, showgrid=False, zeroline=False, visible=False),
        aspectmode='data'
    ), width=width, height=height, margin=dict(l=0,r=0,t=0,b=0), paper_bgcolor="white")

    return fig

