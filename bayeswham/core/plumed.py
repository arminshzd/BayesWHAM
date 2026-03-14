#!/usr/bin/env python3
"""
Convert PLUMED COLVAR files from umbrella sampling simulations into
input files for bayeswham analysis.

Usage:
    bayeswham-plumed <input_file> [options]

The input file lists one umbrella window per line:
    /path/to/COLVAR  <umb_center>  <umb_k>

Lines starting with '#' are treated as comments.
"""

import os
import argparse
import numpy as np
import yaml


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Convert PLUMED COLVAR files to bayeswham input format.\n\n"
            "INPUT FILE FORMAT (one umbrella window per line):\n"
            "  /path/to/COLVAR  <umb_center>  <umb_k>\n\n"
            "Lines beginning with '#' are ignored."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "input_file",
        help="Path to input file listing COLVAR paths and umbrella parameters",
    )
    parser.add_argument(
        "--variable",
        type=str,
        default="cv1",
        help="Column name to extract from COLVAR header (default: cv1)",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="./bayeswham_input",
        help="Output directory for bayeswham files (default: ./bayeswham_input)",
    )
    parser.add_argument(
        "--nbins",
        type=int,
        default=50,
        help="Number of histogram bins (default: 50)",
    )
    parser.add_argument(
        "--skip",
        type=int,
        default=0,
        help="Number of initial frames to skip (default: 0)",
    )
    parser.add_argument(
        "--discard-fraction",
        type=float,
        default=0.1,
        help="Fraction of initial trajectory to discard for equilibration (default: 0.1)",
    )
    parser.add_argument(
        "--temperature",
        type=float,
        default=300.0,
        help="Temperature in Kelvin (default: 300.0)",
    )
    parser.add_argument(
        "--bin-range",
        type=float,
        nargs=2,
        default=None,
        metavar=("MIN", "MAX"),
        help="Manual bin range (default: auto from data)",
    )
    return parser.parse_args()


def read_input_file(input_path):
    """
    Parse the input file listing COLVAR paths and umbrella parameters.

    Each non-comment line must have the format:
        /path/to/COLVAR  <umb_center>  <umb_k>

    Parameters
    ----------
    input_path : str
        Path to the input file

    Returns
    -------
    list of dict
        Each dict contains: 'colvar_path', 'umb_loc', 'umb_k'
    """
    entries = []
    with open(input_path, 'r') as f:
        for lineno, line in enumerate(f, start=1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) != 3:
                raise ValueError(
                    f"Line {lineno} of {input_path}: expected "
                    f"'colvar_path umb_center umb_k', got: {line!r}"
                )
            colvar_path, umb_loc, umb_k = parts
            entries.append({
                'colvar_path': colvar_path,
                'umb_loc': float(umb_loc),
                'umb_k': float(umb_k),
            })
    if not entries:
        raise ValueError(f"No valid entries found in {input_path}")
    return entries


def read_colvar(colvar_path, variable_name, skip=0, discard_fraction=0.1):
    """
    Read a PLUMED COLVAR file and extract the specified variable column.

    Parameters
    ----------
    colvar_path : str
        Path to the COLVAR file
    variable_name : str
        Column name to extract (matched against the '#! FIELDS' header line)
    skip : int
        Number of initial frames to skip
    discard_fraction : float
        Fraction of the (post-skip) trajectory to discard for equilibration

    Returns
    -------
    np.ndarray
        1D array of extracted values
    """
    with open(colvar_path, 'r') as f:
        header = f.readline().strip()

    if not header.startswith("#! FIELDS"):
        raise ValueError(
            f"Expected PLUMED COLVAR header '#! FIELDS ...' in {colvar_path}, "
            f"got: {header!r}"
        )

    fields = header.replace("#! FIELDS", "").split()

    if variable_name not in fields:
        raise ValueError(
            f"Variable '{variable_name}' not found in {colvar_path}. "
            f"Available fields: {fields}"
        )

    col_idx = fields.index(variable_name)
    data = np.loadtxt(colvar_path, comments="#", usecols=col_idx)

    if skip > 0:
        data = data[skip:]

    if discard_fraction > 0:
        n_discard = int(len(data) * discard_fraction)
        data = data[n_discard:]

    return data


def compute_histogram(data, bin_edges):
    """
    Compute integer histogram counts for data over the given bin edges.

    Parameters
    ----------
    data : np.ndarray
        Data values
    bin_edges : np.ndarray
        Bin edges (length = nbins + 1)

    Returns
    -------
    np.ndarray of int
        Histogram counts
    """
    counts, _ = np.histogram(data, bins=bin_edges)
    return counts.astype(int)


def write_bayeswham_inputs(simulations, output_dir, nbins, variable,
                           bin_range=None, temperature=300.0):
    """
    Write all input files required by bayeswham and a ready-to-use config.

    Parameters
    ----------
    simulations : list of dict
        Each dict must contain: 'colvar_path', 'umb_loc', 'umb_k', 'data'
    output_dir : str
        Root output directory
    nbins : int
        Number of histogram bins
    variable : str
        CV name (used only in the summary file)
    bin_range : tuple or None
        (min, max) for histogram edges; auto-detected if None
    temperature : float
        Temperature in Kelvin for the generated bayeswham config

    Returns
    -------
    str
        Path to the generated bayeswham config YAML
    """
    output_dir = os.path.abspath(output_dir)
    hist_dir = os.path.join(output_dir, "hist")
    bias_dir = os.path.join(output_dir, "bias")
    os.makedirs(hist_dir, exist_ok=True)
    os.makedirs(bias_dir, exist_ok=True)

    # Determine bin edges
    all_data = np.concatenate([s['data'] for s in simulations])
    if bin_range is not None:
        data_min, data_max = bin_range
    else:
        data_min = all_data.min()
        data_max = all_data.max()
        padding = 0.05 * (data_max - data_min)
        data_min -= padding
        data_max += padding

    bin_edges = np.linspace(data_min, data_max, nbins + 1)

    print(f"\nHistogram range: [{data_min:.6f}, {data_max:.6f}]")
    print(f"Number of bins:  {nbins}")
    print(f"Bin width:       {(data_max - data_min) / nbins:.8f}")

    # hist_binEdges.txt (single row for 1D)
    bin_edges_file = os.path.join(output_dir, "hist_binEdges.txt")
    with open(bin_edges_file, 'w') as f:
        f.write(" ".join(f"{e:.10e}" for e in bin_edges) + "\n")
    print(f"Wrote: {bin_edges_file}")

    # harmonic_biases.txt
    biases_file = os.path.join(bias_dir, "harmonic_biases.txt")
    with open(biases_file, 'w') as f:
        for i, sim in enumerate(simulations, start=1):
            f.write(f"{i}  {sim['umb_loc']:.10e}  {sim['umb_k']:.10e}\n")
    print(f"Wrote: {biases_file}")

    # Individual histogram files
    total_counts = 0
    for i, sim in enumerate(simulations, start=1):
        counts = compute_histogram(sim['data'], bin_edges)
        total_counts += counts.sum()
        hist_file = os.path.join(hist_dir, f"hist_{i}.txt")
        with open(hist_file, 'w') as f:
            f.write(" ".join(str(c) for c in counts) + "\n")

    print(f"Wrote: {len(simulations)} histogram files to {hist_dir}/")
    print(f"Total histogram counts: {total_counts}")

    # bayeswham config YAML
    bayeswham_output_dir = os.path.join(output_dir, "bayeswham_output")
    config = {
        'umbrella_sampling': {
            'dim': 1,
            'periodicity': [0],
            'temperature': temperature,
        },
        'input_files': {
            'harmonic_biases': biases_file,
            'hist_bin_edges': bin_edges_file,
            'hist_dir': hist_dir,
        },
        'wham': {
            'tolerance': 1.0e-12,
            'max_iterations': 1000000,
        },
        'metropolis_hastings': {
            'steps': 500000,
            'save_modulus': 500,
            'print_modulus': 10000,
            'max_step_size': 1.0e-4,
            'random_seed': 42,
        },
        'prior': {
            'type': 'Dirichlet',
            'alpha': 2.0,
        },
        'output_dir': bayeswham_output_dir,
    }

    config_file = os.path.join(output_dir, "bayeswham_config.yaml")
    with open(config_file, 'w') as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)
    print(f"Wrote: {config_file}")

    # Human-readable summary
    summary_file = os.path.join(output_dir, "simulation_summary.txt")
    with open(summary_file, 'w') as f:
        f.write("# Umbrella Sampling Simulation Summary\n")
        f.write(f"# Variable: {variable}\n")
        f.write(f"# Temperature: {temperature} K\n")
        f.write(f"# Number of simulations: {len(simulations)}\n")
        f.write(f"# Histogram bins: {nbins}\n")
        f.write(f"# Bin range: [{data_min:.6f}, {data_max:.6f}]\n")
        f.write("#\n")
        f.write("# Index  UMB_LOC        UMB_K          N_samples      COLVAR\n")
        for i, sim in enumerate(simulations, start=1):
            f.write(
                f"{i:7d}  {sim['umb_loc']:12.6f}  {sim['umb_k']:12.2f}  "
                f"{len(sim['data']):12d}  {sim['colvar_path']}\n"
            )
    print(f"Wrote: {summary_file}")

    return config_file


def main():
    args = parse_args()

    print("=" * 60)
    print("PLUMED COLVAR to BayesWHAM Converter")
    print("=" * 60)
    print(f"Input file:       {args.input_file}")
    print(f"Variable:         {args.variable}")
    print(f"Skip frames:      {args.skip}")
    print(f"Discard fraction: {args.discard_fraction:.1%}")
    print(f"Temperature:      {args.temperature} K")
    print(f"Output:           {args.output_dir}")
    print("=" * 60 + "\n")

    # Parse input file
    entries = read_input_file(args.input_file)
    print(f"Found {len(entries)} entries in {args.input_file}\n")

    # Load each COLVAR
    simulations = []
    for entry in entries:
        colvar_path = entry['colvar_path']
        if not os.path.exists(colvar_path):
            print(f"Warning: Skipping {colvar_path} — file not found")
            continue
        try:
            data = read_colvar(
                colvar_path,
                args.variable,
                skip=args.skip,
                discard_fraction=args.discard_fraction,
            )
            simulations.append({**entry, 'data': data})
            print(
                f"Loaded: {colvar_path}  "
                f"(center={entry['umb_loc']:.4f}, k={entry['umb_k']:.4f}, "
                f"N={len(data)})"
            )
        except Exception as e:
            print(f"Warning: Skipping {colvar_path} — {e}")

    if not simulations:
        raise RuntimeError("No valid simulations loaded. Check your input file and COLVAR paths.")

    # Sort by umbrella center
    simulations.sort(key=lambda x: x['umb_loc'])

    print(f"\nLoaded {len(simulations)} simulations")

    # Write bayeswham inputs
    config_file = write_bayeswham_inputs(
        simulations,
        args.output_dir,
        args.nbins,
        args.variable,
        bin_range=args.bin_range,
        temperature=args.temperature,
    )

    print("\n" + "=" * 60)
    print("Done! To run bayeswham:")
    print(f"  bayeswham --config {config_file}")
    print("=" * 60)


if __name__ == "__main__":
    main()
