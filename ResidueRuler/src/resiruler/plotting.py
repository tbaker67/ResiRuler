import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import ast
import os
import plotly.graph_objects as go
from src.resiruler.distance_calc import CompareDistanceMatrix


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

def plot_colorbar(min_val, max_val, cmap_name="plasma", label="Distance (Å)"):
    fig, ax = plt.subplots(figsize=(6, 1))
    fig.subplots_adjust(bottom=0.5)

    cmap = plt.get_cmap(cmap_name)
    norm = plt.Normalize(vmin=min_val, vmax=max_val)
    cb1 = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),
                       cax=ax, orientation='horizontal')
    cb1.set_label(label)
    return fig

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
        fig =  plot_distance_difference_plotly(merged)
        return fig, merged
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
    merged = merged.copy()
    merged["Direction"] = merged["Distance_Diff"].apply(lambda x: "Increase" if x > 0 else "Decrease")

    #Make line breaks work in web interface
    merged['Stacked_Pair_label'] = merged['Pair_label'].str.replace(r'\s+', '<br>', regex=True)

    
    merged['x'] = list(range(len(merged)))
    color_map = {'Increase': 'orange', 'Decrease': 'blue'}
    bar_colors = merged['Direction'].map(color_map)

    fig = go.Figure(data=[
        go.Bar(
            x=merged['x'],
            y=merged['Distance_Diff'],
            marker_color=bar_colors,
            hovertext=merged['Pair_label'],
            hoverinfo='text+y'
        )
    ])

    fig.update_layout(
        title="Distance Differences",
        xaxis=dict(
            tickmode='array',
            tickvals=merged['x'],
            ticktext=merged['Stacked_Pair_label'],
            tickangle=0,
            tickfont=dict(size=10),
        ),
        yaxis_title="Δ Distance (Å)"
    )

    return fig

def plot_interactive_contact_map(matrix, lower_threshold=None, upper_threshold=None, title=None, min=None, max=None):
    mat, index_map = matrix.mat, matrix.index_map
    if lower_threshold is not None and upper_threshold is not None:
        mat = np.where(((mat < upper_threshold) & (mat > lower_threshold)), mat, np.nan)
    if isinstance(matrix, CompareDistanceMatrix):
        hovertemplate = "Residue 1: %{x}<br>Residue 2: %{y}<br>ΔDistance: %{z:.2f} Å<extra></extra>"
    else:
        hovertemplate = "Residue 1: %{x}<br>Residue 2: %{y}<br>Distance: %{z:.2f} Å<extra></extra>"
        
    labels = [f"{chain}-{resid[1]}{resid[2]}" for (chain, resid) in index_map.keys()]
    vmax = np.nanmax(np.abs(mat))
    fig = go.Figure(data=go.Heatmap(
        z=mat,
        x=labels,
        y=labels,
        colorscale="RdBu_r" if isinstance(matrix, CompareDistanceMatrix) else "viridis",
        zmid=0.0 if isinstance(matrix, CompareDistanceMatrix) else None,
        zmin=min if min is not None else -vmax,
        zmax=max if max is not None else vmax,
        colorbar=dict(title="ΔDistance (Å)" if isinstance(matrix, CompareDistanceMatrix) else "Distance (Å)"),
        
        hovertemplate=hovertemplate
    ))

    fig.update_layout(
        title=title or "Residue Contact Map",
        xaxis_title="Residue",
        yaxis_title="Residue",
        autosize=False,
        width=800,
        height=800
    )

    return fig

def plot_all_matrices_ensemble(ref_dm, tgt_dms_dict, compare_dms_dict, lower_threshold=None, upper_threshold=None):
    """
    This takes outputs produce by the EnsembleMapper calc_matrices functions and produces dictionaries the contain plotly objects for interactive contact maps
    THe plots are all set to be on the same value scale, and the dicitonaries can be indexed 
    """
    
    #get all the distance/compare matrices
    dist_mats = [dm.mat for dm in tgt_dms_dict.values()]
    dist_mats.append(ref_dm.mat)
    diff_mats = [cm.mat for cm in compare_dms_dict.values()]

    #get the global min and max values across all structures
    dist_min_vals = [np.nanmin(m) for m in dist_mats if not np.isnan(m).all()]
    dist_max_vals = [np.nanmax(m) for m in dist_mats if not np.isnan(m).all()]
    dist_min = min(dist_min_vals) if dist_min_vals else 0 #TODO: Adjust default values, or just throw error/print message
    dist_max = max(dist_max_vals) if dist_max_vals else 100

    
    diff_min_vals = [np.nanmin(m) for m in diff_mats if not np.isnan(m).all()]
    diff_max_vals = [np.nanmax(m) for m in diff_mats if not np.isnan(m).all()]
    if diff_min_vals and diff_max_vals:
        diff_abs_max = max(abs(min(diff_min_vals)), abs(max(diff_max_vals)))
    else:
        diff_abs_max = 1  # default fallback

    tgt_figs_dict = {}
    compare_figs_dict = {}

    # plot reference
    ref_fig = plot_interactive_contact_map(
        ref_dm,
        lower_threshold=lower_threshold,
        upper_threshold=upper_threshold,
        title="Reference Distance Matrix",
        min=dist_min,
        max=dist_max,
    )

    # plot targets
    for name, tgt_dm in tgt_dms_dict.items():
        tgt_figs_dict[name] = plot_interactive_contact_map(
            tgt_dm,
            lower_threshold=lower_threshold,
            upper_threshold=upper_threshold,
            title=f"Target Distance Matrix: {name}",
            min=dist_min,
            max=dist_max,
        )

    # Plot comparison matrices (differences) on symmetric scale
    for name, compare_dm in compare_dms_dict.items():
        compare_figs_dict[name] = plot_interactive_contact_map(
            compare_dm,
            lower_threshold=lower_threshold,
            upper_threshold=upper_threshold,
            title=f"Comparison Distance Matrix: {name}",
            min=-diff_abs_max,
            max=diff_abs_max,
        )

    return ref_fig, tgt_figs_dict, compare_figs_dict