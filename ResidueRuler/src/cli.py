import argparse
import os
import yaml
import pandas as pd
from src.resiruler.structure_parsing import load_structure, extract_CA_coords
from src.resiruler.distance_calc import get_header_indices, read_data
from src.resiruler.chimera_export import generate_chimerax_script
from src.resiruler.distance_difference_plot import plot_distance_difference

def default_command(args):
    '''
    Run the standard program which creates a CSV of distances between specified Residues
    '''
    if not all([args.input, args.structure, args.output, args.config]):
        raise ValueError("All of --input, --structure, --output, and --config are required for 'run' mode.")

    with open(args.config) as f:
        config = yaml.safe_load(f)

    chain_mapping = config["chain_mapping"]
    thresholds = config.get("color_thresholds", {"green": 20, "yellow": 33})

    structure = load_structure(args.structure)
    index_map, coords = extract_CA_coords(structure)

    df = pd.read_excel(args.input, sheet_name=0).replace('', pd.NA).dropna(how='all').reset_index(drop=True)
    header_indices = get_header_indices(df)
    final_df = read_data(df, header_indices, chain_mapping, index_map, coords)
    final_df.to_csv(args.output, index=False)

    print(f"[INFO] CSV written to {args.output}")

    if args.visualize:
        chx_output = args.output.replace(".csv", "_chimerax.cxc")
        chains_for_viz = sorted(set(chain for sublist in chain_mapping.values() for chain in sublist))
        generate_chimerax_script(args.output, chx_output, chains=chains_for_viz, thresholds=thresholds)
        print(f"[INFO] ChimeraX script written to {chx_output}")


def compare_command(args):
    '''
    Run the compare program which will take in two generated CSV files and produce a CSV file of the difference in distances between residues in both files
    '''
    if not all(map(os.path.exists, [args.csv1, args.csv2])):
        raise FileNotFoundError("One or both comparison CSV files do not exist.")
    plot_distance_difference(args.csv1, args.csv2, output_path=args.plot_out)
    print(f"[INFO] Distance difference plot saved to {args.plot_out}")


def main():
    parser = argparse.ArgumentParser(description="ResiRuler: Analyze residue distances from structure and annotation.")
    subparsers = parser.add_subparsers(dest='command', required=True)

    # Default Subcommand
    run_parser = subparsers.add_parser('run', help='Run ResiRuler with input data and structure')
    run_parser.add_argument('-i', '--input', required=True, help='Path to Excel file')
    run_parser.add_argument('-s', '--structure', required=True, help='Structure file (.cif)')
    run_parser.add_argument('-o', '--output', required=True, help='Output CSV file')
    run_parser.add_argument('-c', '--config', required=True, help='YAML config with chain mapping and thresholds')
    run_parser.add_argument('-v', '--visualize', action='store_true', help='Generate ChimeraX visualization script')
    run_parser.set_defaults(func=default_command)

    # Compare Subcommand
    compare_parser = subparsers.add_parser('compare', help='Compare two result CSVs and plot differences')
    compare_parser.add_argument('csv1', help='First CSV file to compare')
    compare_parser.add_argument('csv2', help='Second CSV file to compare')
    compare_parser.add_argument('--plot-out', default='distance_diff_plot.png', help='Output image path')
    compare_parser.set_defaults(func=compare_command)

    # Parse and Run Desired Mode
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()