import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

def plot_distance_difference(csv1_path, csv2_path, output_path="distance_diff.png", sort=True):
    """
    Function to compare distance from two files and plot the differences.
    """
    
    df1 = pd.read_csv(csv1_path)
    df2 = pd.read_csv(csv2_path)

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

    
    colors = ['orange' if diff > 0 else 'blue' for diff in merged['Distance_Diff']]

    
    merged['Pair_label'] = merged['Chain1_Residue1'] + '\n' + merged['Chain2_Residue2']

    
    
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

    columns_to_export = ['Chain1_Residue1', 'Chain2_Residue2', 'Distance_1','Distance_2','Distance_Diff']
    merged.to_csv(output_path[:-3] + 'csv', index=False, columns=columns_to_export)
    
    
def plot_movement_shift(df, output_path='movement_plot.png'):
    coords = df['Coord1'].tolist()
    coords = np.array(coords)

    x, y, z = coords[:,0], coords[:,1], coords[:,2]
    distances = df['Distance']
    # Normalize distances for colormap
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



def plot_movement_vectors(df, output_path='movement_vectors.png',cmap_name='plasma', arrow_color='black', min_distance=-0.1):
    # Filter to show only vectors above a threshold
    df_filtered = df[df['Distance'] >= min_distance]

    coords = np.array(df_filtered['Coord1'].tolist())
    vectors = np.array(df_filtered['Diff_Vec'].tolist())
    distances = df_filtered['Distance'].values

    
    x, y, z = coords[:, 0], coords[:, 1], coords[:, 2]
    u, v, w = vectors[:, 0], vectors[:, 1], vectors[:, 2]

    # Normalize distances for color mapping
    norm = plt.Normalize(vmin=distances.min(), vmax=distances.max())
    cmap = plt.cm.get_cmap(cmap_name)
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

    ax.set_title(f'3D Residue Movement Vectors (N={len(df_filtered)})')
    ax._axis3don = False
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.show()
    plt.close()