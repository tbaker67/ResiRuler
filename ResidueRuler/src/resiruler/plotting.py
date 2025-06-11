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
    cmap = plt.cm.winter  #Maybe add option to change color scheme
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
