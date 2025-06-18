import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import ast
import os
import plotly.graph_objects as go
import plotly.express as px
from plotly.colors import sample_colorscale
import matplotlib.colors as mcolors

def parse_coord_column(series):
    return np.array([np.array(ast.literal_eval(s)) for s in series])


def remap_residue(residue, chain_map):
    if '_' not in residue:
        return residue

    try:
        chain, resnum = residue.split('_', 1)
        mapped_chain = chain_map[chain]  # this will raise KeyError if missing
        return f"{mapped_chain}_{resnum}"
    except KeyError:
        print(f"[WARNING] Chain '{chain}' not found in chain_map. Using original ID '{residue}'.")
        return residue
    except Exception as e:
        print(f"[ERROR] Failed to remap residue ID '{residue}': {e}")
        return residue


def plot_distance_difference(csv1_path, csv2_path, output_path="distance_diff.png", chain_map=None, sort=True, plotly=False):
    """
    Function to compare distance from two files and plot the differences.
    """
    
    df1 = pd.read_csv(csv1_path)
    df2 = pd.read_csv(csv2_path)

    if chain_map:
        df2['Chain1_Residue1'].apply(lambda r: remap_residue(r, chain_map))
        df2['Chain2_Residue2'].apply(lambda r: remap_residue(r, chain_map))

    #Put the two tables together and calculate the distance differences
    key_cols = ['Chain1_Residue1', 'Chain2_Residue2']
    merged = pd.merge(df1, df2, on=key_cols, suffixes=('_1', '_2'))
    merged = merged.dropna(subset=['Distance_1', 'Distance_2'])
    if merged.empty:
        print("[WARNING] No valid distance pairs to compare. Plot not generated.")
        return
    merged['Distance_Diff'] = merged['Distance_1'] - merged['Distance_2']

    #Sort by distance difference if specified, otherwise sort in alphabetcial order by pair label
    if sort:
        merged = merged.sort_values(by='Distance_Diff', key=lambda x: x, ascending=False)
    else:
        merged = merged.sort_values(by='Pair_label')

    
    

    
    merged['Pair_label'] = merged['Chain1_Residue1'] + '\n' + merged['Chain2_Residue2']

    if plotly:
        return plot_distance_difference_plotly(merged)
    else:
        plot_distance_difference_matplot(merged, output_path, csv1_path, csv2_path)
    

    columns_to_export = ['Chain1_Residue1', 'Chain2_Residue2', 'Distance_1','Distance_2','Distance_Diff']
    merged.to_csv(output_path[:-3] + 'csv', index=False, columns=columns_to_export)

def plot_distance_difference_matplot(merged,output_path, csv1_path, csv2_path):
    colors = ['orange' if diff > 0 else 'blue' for diff in merged['Distance_Diff']]

    x = np.arange(len(merged))
    bar_width = 0.4

    # Plot
    plt.figure(figsize=(max(7, len(merged) * 0.3), 6))
    bars = plt.bar(x, merged['Distance_Diff'], width=bar_width, color=colors)

    plt.axhline(0, color='gray', linewidth=0.8)

    plt.xticks(x, merged['Pair_label'], rotation=0)
    plt.ylabel('Distance Difference (Å)')
    plt.title(f"Distance Differences: {os.path.basename(csv1_path)} - {os.path.basename(csv2_path)}")

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

def plot_distance_difference_plotly(merged):
    fig = px.bar(
        merged,
        x='Pair_label',
        y='Distance_Diff',
        color=merged['Distance_Diff'] > 0,
        color_discrete_map={True: 'orange', False: 'blue'},
        labels={'Pair_label': 'Residue Pair', 'Distance_Diff': 'Δ Distance (Å)'},
        title="Distance Differences"
    )
    fig.update_layout(xaxis_tickangle=-45)
    return fig
    
def plot_movement_shift(df, output_path='movement_plot.png', plotly=False):
    coords = parse_coord_column(df['Coord1'])

    x, y, z = coords[:,0], coords[:,1], coords[:,2]
    distances = df['Distance']
    
    if plotly:
        return plotly_movement_shift(x,y,z,distances)
    else:
        matplot_movement_shift(x,y,z, distances,output_path)
    
def matplot_movement_shift(x,y,z, distances,output_path):
    norm = plt.Normalize(distances.min(), distances.max())
    cmap = plt.cm.plasma  #Maybe add option to change color scheme
    colors = cmap(norm(distances))

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')

    p = ax.scatter(x, y, z, c=distances, cmap=cmap, norm=norm, s=40)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('3D Residue Positions Colored by Distance Shift')

    # Add colorbar
    cbar = fig.colorbar(p, ax=ax, shrink=0.6)
    cbar.set_label('Distance')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.show()
    plt.close()


def plotly_movement_shift(x,y,z,distances):
    
    fig = go.Figure(data=[go.Scatter3d(
        x=x, y=y, z=z,
        mode='markers',
        marker=dict(
            size=5,
            color=distances,
            colorscale='Plasma',
            colorbar=dict(title='Distance Shift')
        )
    )])

    fig.update_layout(
        title='3D Residue Positions Colored by Distance Shift',
        scene=dict(xaxis_title='X', yaxis_title='Y', zaxis_title='Z'),
        width=800,
        height=600
    )
    return fig

def set_equal_3d_axes(ax, coords, zoom_out_factor=1.2):
    """
    Set equal aspect ratio for 3D plot and apply zoom out factor.
    """
    x, y, z = coords[:, 0], coords[:, 1], coords[:, 2]

    # Find midpoints and max range
    mid_x, mid_y, mid_z = np.mean(x), np.mean(y), np.mean(z)
    max_range = np.ptp(coords, axis=0).max() * zoom_out_factor / 2

    # Set limits
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)



def plot_movement_vectors(df, output_path='movement_vectors.png', min_distance=-0.1, plotly=False):
    # Filter to show only vectors above a threshold
    df_filtered = df[df['Distance'] >= min_distance]

    coords = parse_coord_column(df['Coord1'])
    vectors = parse_coord_column(df['Diff_Vec'])
    distances = df_filtered['Distance'].values

    
    x, y, z = coords[:, 0], coords[:, 1], coords[:, 2]
    u, v, w = vectors[:, 0], vectors[:, 1], vectors[:, 2]

    if plotly:
        return plotly_movement_vectors(x,y,z,u,v,w,distances)
    else:
        matplot_movement_vectors(x,y,z,u,v,w,distances,coords,output_path)


def matplot_movement_vectors(x,y,z,u,v,w,distances,coords,output_path):
    # Normalize distances for color mapping
    norm = plt.Normalize(vmin=distances.min(), vmax=distances.max())
    cmap = plt.cm.get_cmap('plasma')
    colors = cmap(norm(distances))

    # Plot
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    for i in range(len(x)):
        ax.quiver(x[i], y[i], z[i], u[i], v[i], w[i],
                  length=1.0,
                  normalize=False,
                  color=colors[i],
                  linewidth=0.7,
                  arrow_length_ratio=0.2)

    set_equal_3d_axes(ax, coords, zoom_out_factor=1.6)
                      
    # Add color bar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.6)
    cbar.set_label('Distance')

    ax.set_title(f'3D Residue Movement Vectors')
    ax._axis3don = False
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.show()
    plt.close()


def plotly_movement_vectors(x, y, z, u, v, w, distances, min_distance=-0.1):
    


    
    norm = plt.Normalize(vmin=distances.min(), vmax=distances.max())
    cmap = plt.cm.get_cmap('plasma')

    # Convert colormap to hex colors
    colors = [mcolors.to_hex(cmap(norm(dist))) for dist in distances]
    x_end = x + u
    y_end = y + v
    z_end = z + w

    fig = go.Figure()

    # Lines for vectors
    for i in range(len(x)):
        fig.add_trace(go.Scatter3d(
            x=[x[i], x_end[i]],
            y=[y[i], y_end[i]],
            z=[z[i], z_end[i]],
            mode='lines',
            line=dict(color=colors[i] if colors is not None else 'blue', width=4),
            showlegend=False
        ))

    # Optionally, add markers at vector heads
    fig.add_trace(go.Scatter3d(
        x=x_end, y=y_end, z=z_end,
        mode='markers',
        marker=dict(size=3, color=colors if colors is not None else 'blue'),
        name='Vector heads'
    ))

    fig.update_layout(scene=dict(
        xaxis_title='X', yaxis_title='Y', zaxis_title='Z',
    ))

    return fig