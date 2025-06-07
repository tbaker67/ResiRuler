import pandas as pd
import numpy as np
import argparse
import yaml
from src.resiruler.distance_calc import get_header_indices, read_data
from src.resiruler.structure_parsing import load_structure, extract_CA_coords
from src.resiruler.chimera_export import generate_chimerax_script
from src.resiruler.distance_difference_plot import plot_distance_difference

import os

def parse_args():
    parser = argparse.ArgumentParser(description="ResiRuler: Calculate and visualize distances between residues.")
    parser.add_argument("-i", "--input", help="Path to Excel file.")
    parser.add_argument("-s", "--structure", help="Structure file (.cif).")
    parser.add_argument("-o", "--output", help="Output CSV file.")
    parser.add_argument("-c", "--config", help="YAML config with chain mapping and thresholds.")
    parser.add_argument("-v", "--visualize", action="store_true", help="Generate ChimeraX visualization script.")
    parser.add_argument("--compare", nargs=2, metavar=('CSV1', 'CSV2'), help="Compare distances between two result CSVs and plot differences.")
    parser.add_argument("--plot-out", default="distance_diff_plot.png", help="Output image path for comparison plot (used with --compare).")
    return parser.parse_args()

def main():
    args = parse_args()

    if args.compare:
        csv1, csv2 = args.compare
        if not all(map(os.path.exists, [csv1, csv2])):
            raise FileNotFoundError("One or both comparison CSV files do not exist.")
        plot_distance_difference(csv1, csv2, output_path=args.plot_out)
        print(f"[INFO] Distance difference plot saved to {args.plot_out}")
        return

    if not all([args.input, args.structure, args.output, args.config]):
        raise ValueError("Input (-i), structure (-s), output (-o), and config (-c) arguments are required unless using --compare mode.")

    with open(args.config) as f:
        config = yaml.safe_load(f)

    chain_mapping = config["chain_mapping"]
    #Use default thresholds if not provided in config
    thresholds = config.get("color_thresholds", {"green": 20, "yellow": 33})

    structure = load_structure(args.structure)
    if args.structure.endswith('.cif'):
        index_map, coords = extract_CA_coords(structure)

    df = pd.read_excel(args.input, sheet_name=0).replace('', pd.NA).dropna(how='all').reset_index(drop=True)
    header_indices = get_header_indices(df)

    final_df = read_data(df, header_indices, chain_mapping, index_map, coords)
    final_df.to_csv(args.output, index=False)
    print(f"CSV written to {args.output}")

    if args.visualize:
        chx_output = args.output.replace(".csv", "_chimerax.cxc")
        chains_for_viz = sorted(set(chain for sublist in chain_mapping.values() for chain in sublist))
        generate_chimerax_script(args.output, chx_output, chains=chains_for_viz, thresholds=thresholds)
        print(f"ChimeraX script written to {chx_output}")

if __name__ == "__main__":
    main()