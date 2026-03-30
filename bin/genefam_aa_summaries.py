#!/usr/bin/env python3
import os
import math
import statistics
import pandas as pd
import argparse
from Bio import SeqIO
from Bio.SeqUtils import ProtParam
from itertools import product

# Set up command line argument parsing
parser = argparse.ArgumentParser(
    description="Analyze multiple sequence alignment in FASTA format"
)
parser.add_argument(
    "input",
    type=str,
    help="Path to input MSA file in FASTA format",
)
args = parser.parse_args()

# Define the base directories
across_fam_basedir = "aa-summary-stats/across-family-summaries/"
per_fam_basedir = "aa-summary-stats/per-family-summaries/"

# Define the order of amino acids and columns
aa_order = list("ACDEFGHIKLMNPQRSTVWY")

aa_counts_columns = ["id"] + [f"aa_composition_{aa}" for aa in aa_order]
aa_perc_columns = ["id"] + [f"aa_composition_percent_{aa}" for aa in aa_order] + ["composition_entropy"]
summary_columns = [
    "id",
    "molecular_weight",
    "aromaticity",
    "instability",
    "flexibility",
    "gravy_kd",
    "gravy_bm",
    "gravy_ro",
    "isoelectric_point",
    "charge_at_pH_3",
    "charge_at_pH_5",
    "charge_at_pH_7",
    "charge_at_pH_9",
    "charge_at_pH_11",
    "helix_fract",
    "turn_fract",
    "sheet_fract",
    "molar_ext_coef_cysteines",
    "molar_ext_coef_cystines",
]

# Define data types and statistics
data_types = [
    {
        "name": "counts",
        "columns": aa_counts_columns,
        "filename_prefix": "aa_{}_counts_per_fam.csv",
    },
    {
        "name": "percents",
        "columns": aa_perc_columns,
        "filename_prefix": "aa_{}_percents_per_fam.csv",
    },
    {
        "name": "summaries",
        "columns": summary_columns,
        "filename_prefix": "aa_{}_summaries_per_fam.csv",
    },
]

stats = ["mean", "median", "stdev"]

# Residue sets for windowed composition-based properties
HELIX_RESIDUES = set("EMALK")
TURN_RESIDUES = set("NPGSD")
SHEET_RESIDUES = set("VIYFWLT")
AROMATIC_RESIDUES = set("FWY")

# Standard pKa values for per-residue charge calculation (Lehninger)
CHARGED_PKA = {
    "D": 3.65, "E": 4.25, "C": 8.18, "Y": 10.07,
    "H": 6.00, "K": 10.53, "R": 12.48,
}

# Column definitions for windowed SD and autocorrelation tables
sd_columns = [
    "id", "flexibility_sd", "gravy_kd_sd", "gravy_bm_sd", "gravy_ro_sd",
    "aromaticity_sd",
    "charge_at_pH_3_sd", "charge_at_pH_5_sd", "charge_at_pH_7_sd",
    "charge_at_pH_9_sd", "charge_at_pH_11_sd",
    "helix_fract_sd", "turn_fract_sd", "sheet_fract_sd",
    "composition_entropy_sd",
]
autocorr_columns = [c.replace("_sd", "_autocorr") for c in sd_columns]

# Create dataframe_specs using a list comprehension
dataframe_specs = [
    {
        "key": f"{stat}_{data_type['name']}_df",
        "columns": data_type["columns"],
        "filename": os.path.join(
            across_fam_basedir, data_type["filename_prefix"].format(stat)
        ),
    }
    for data_type in data_types
    for stat in stats
]

# Initialize the dataframes dictionary using dataframe_specs
dataframes = {spec["key"]: pd.DataFrame() for spec in dataframe_specs}


def shannon_entropy(seq, aa_set="ACDEFGHIKLMNPQRSTVWY"):
    """Shannon entropy (log2) of amino acid frequencies in a sequence."""
    counts = {aa: 0 for aa in aa_set}
    total = 0
    for c in seq:
        if c in counts:
            counts[c] += 1
            total += 1
    if total == 0:
        return 0.0
    H = 0.0
    for count in counts.values():
        if count > 0:
            p = count / total
            H -= p * math.log2(p)
    return H


def residue_charge(aa, pH):
    """Net charge contribution of a single residue at a given pH."""
    if aa not in CHARGED_PKA:
        return 0.0
    pKa = CHARGED_PKA[aa]
    if aa in ("D", "E", "C", "Y"):  # acidic
        return -1.0 / (1.0 + 10.0 ** (pKa - pH))
    else:  # basic (H, K, R)
        return 1.0 / (1.0 + 10.0 ** (pH - pKa))


def compute_windowed_profiles(analyzer, seq, window=9):
    """Return dict of property_name -> list of per-window values."""
    n = len(seq)
    profiles = {}

    if n < window:
        return profiles  # empty -> caller handles short proteins

    # Scale-based (protein_scale returns per-window values directly)
    for name, scale in [
        ("flexibility", ProtParam.ProtParamData.Flex),
        ("gravy_kd", ProtParam.ProtParamData.kd),
        ("gravy_bm", ProtParam.ProtParamData.bm),
        ("gravy_ro", ProtParam.ProtParamData.ro),
    ]:
        profiles[name] = analyzer.protein_scale(scale, window)

    # Composition-based sliding windows
    aromatic, helix, turn, sheet = [], [], [], []
    for i in range(n - window + 1):
        win = seq[i : i + window]
        w = len(win)
        aromatic.append(sum(1 for c in win if c in AROMATIC_RESIDUES) / w)
        helix.append(sum(1 for c in win if c in HELIX_RESIDUES) / w)
        turn.append(sum(1 for c in win if c in TURN_RESIDUES) / w)
        sheet.append(sum(1 for c in win if c in SHEET_RESIDUES) / w)
    profiles["aromaticity"] = aromatic
    profiles["helix_fract"] = helix
    profiles["turn_fract"] = turn
    profiles["sheet_fract"] = sheet

    # Charge at different pH values (per-residue -> rolling window mean)
    for pH in [3.0, 5.0, 7.0, 9.0, 11.0]:
        per_res = [residue_charge(aa, pH) for aa in seq]
        profiles[f"charge_at_pH_{int(pH)}"] = [
            sum(per_res[i : i + window]) / window
            for i in range(n - window + 1)
        ]

    # Composition entropy per window
    entropy_vals = []
    for i in range(n - window + 1):
        win = seq[i : i + window]
        entropy_vals.append(shannon_entropy(win))
    profiles["composition_entropy"] = entropy_vals

    return profiles


def summarize_profiles(profiles):
    """Compute SD and lag-1 autocorrelation for each windowed profile."""
    sds = {}
    autocorrs = {}
    for name, vals in profiles.items():
        if len(vals) > 1:
            sds[f"{name}_sd"] = statistics.pstdev(vals)
            # Lag-1 Pearson autocorrelation
            n = len(vals)
            mean = sum(vals) / n
            var = sum((v - mean) ** 2 for v in vals) / n
            if var > 0:
                cov = sum(
                    (vals[i] - mean) * (vals[i + 1] - mean) for i in range(n - 1)
                ) / (n - 1)
                autocorrs[f"{name}_autocorr"] = cov / var
            else:
                autocorrs[f"{name}_autocorr"] = 0.0
        else:
            sds[f"{name}_sd"] = 0.0
            autocorrs[f"{name}_autocorr"] = 0.0
    return sds, autocorrs


# A function to calculate the grand average for different protein scales
def grand_average(analyzer, scale):
    scale_res = analyzer.protein_scale(scale, 9)
    return sum(scale_res) / len(scale_res) if scale_res else 0


# Function to generate all possible alternative sequences
def generate_alternative_sequences(seq):
    ambiguous_aa_mapping = {"B": ["D", "N"], "J": ["I", "L"], "Z": ["E", "Q"]}
    positions_and_replacements = [
        (i, ambiguous_aa_mapping[aa])
        for i, aa in enumerate(seq)
        if aa in ambiguous_aa_mapping
    ]
    if not positions_and_replacements:
        return [seq]
    positions, replacements_list = zip(*positions_and_replacements)
    alternative_sequences = []
    for replacements in product(*replacements_list):
        alternative_seq_list = list(seq)
        for position, replacement in zip(positions, replacements):
            alternative_seq_list[position] = replacement
        alternative_sequences.append("".join(alternative_seq_list))
    return alternative_sequences


# Get the input MSA file
msa_file = args.input

# Validate input file exists
if not os.path.isfile(msa_file):
    raise FileNotFoundError(f"Input file not found: {msa_file}")

# Create the output directories
os.makedirs(across_fam_basedir, exist_ok=True)
os.makedirs(os.path.join(per_fam_basedir, "aa-counts"), exist_ok=True)
os.makedirs(os.path.join(per_fam_basedir, "aa-proportions"), exist_ok=True)
os.makedirs(os.path.join(per_fam_basedir, "aa-physical-properties"), exist_ok=True)
os.makedirs(os.path.join(per_fam_basedir, "aa-physical-property-sds"), exist_ok=True)
os.makedirs(os.path.join(per_fam_basedir, "aa-physical-property-autocorr"), exist_ok=True)


# Function to flatten nested dictionaries
def flatten_stats(stats):
    flat_stats = {}
    for key, value in stats.items():
        if isinstance(value, dict):
            for subkey, subvalue in value.items():
                flat_stats[f"{key}_{subkey}"] = subvalue
        else:
            flat_stats[key] = value
    return flat_stats


# Function to calculate statistics
def calculate_stats(df, gene_family_name):
    stats_functions = {
        "mean": df.mean,
        "median": df.median,
        "stdev": df.std,
    }
    stats = {}
    for stat_name, func in stats_functions.items():
        stat_dict = func().to_dict()
        stat_dict["id"] = gene_family_name
        stats[stat_name] = stat_dict
    return stats


# Function to process each MSA file
def process_msa(msa_file):
    # Read in the MSA file in FASTA format
    alignment = list(SeqIO.parse(msa_file, "fasta"))

    if len(alignment) == 0:
        import sys
        print(f"WARNING: No sequences found in {msa_file}. Producing empty output files.", file=sys.stderr)
        gene_family_name = os.path.splitext(os.path.basename(msa_file))[0].split("_")[0]
        pd.DataFrame(columns=aa_counts_columns).to_csv(
            os.path.join(per_fam_basedir, "aa-counts", f"{gene_family_name}_aa_composition_counts.csv"), index=False)
        pd.DataFrame(columns=aa_perc_columns).to_csv(
            os.path.join(per_fam_basedir, "aa-proportions", f"{gene_family_name}_aa_composition_percentages.csv"), index=False)
        pd.DataFrame(columns=summary_columns).to_csv(
            os.path.join(per_fam_basedir, "aa-physical-properties", f"{gene_family_name}_summary_statistics.csv"), index=False)
        pd.DataFrame(columns=sd_columns).to_csv(
            os.path.join(per_fam_basedir, "aa-physical-property-sds", f"{gene_family_name}_summary_statistics_sd.csv"), index=False)
        pd.DataFrame(columns=autocorr_columns).to_csv(
            os.path.join(per_fam_basedir, "aa-physical-property-autocorr", f"{gene_family_name}_summary_statistics_autocorr.csv"), index=False)
        return

    sequence_stats = []
    sequence_sds = []
    sequence_autocorrs = []

    # Define the properties to calculate
    properties = [
        ("molecular_weight", lambda pa, seq: pa.molecular_weight()),
        ("aromaticity", lambda pa, seq: pa.aromaticity()),
        ("instability", lambda pa, seq: pa.instability_index()),
        (
            "flexibility",
            lambda pa, seq: grand_average(pa, ProtParam.ProtParamData.Flex),
        ),
        (
            "gravy_kd",
            lambda pa, seq: grand_average(pa, ProtParam.ProtParamData.kd),
        ),
        (
            "gravy_bm",
            lambda pa, seq: grand_average(pa, ProtParam.ProtParamData.bm),
        ),
        (
            "gravy_ro",
            lambda pa, seq: grand_average(pa, ProtParam.ProtParamData.ro),
        ),
        ("isoelectric_point", lambda pa, seq: pa.isoelectric_point()),
        ("helix_fract", lambda pa, seq: pa.secondary_structure_fraction()[0]),
        ("turn_fract", lambda pa, seq: pa.secondary_structure_fraction()[1]),
        ("sheet_fract", lambda pa, seq: pa.secondary_structure_fraction()[2]),
        (
            "molar_ext_coef_cysteines",
            lambda pa, seq: pa.molar_extinction_coefficient()[0],
        ),
        (
            "molar_ext_coef_cystines",
            lambda pa, seq: pa.molar_extinction_coefficient()[1],
        ),
    ]

    # Define pH values for charge calculations
    pH_values = [3.0, 5.0, 7.0, 9.0, 11.0]

    for seq_record in alignment:
        # Remove gaps and revert lowercase chars to upper
        original_seq = str(seq_record.seq).replace("-", "").upper()

        # Generate alternative sequences for each sequence in the MSA
        alternative_sequences = generate_alternative_sequences(original_seq)

        # Initialize cumulative stats
        cumulative_stats = {prop: 0 for prop, _ in properties}
        cumulative_stats.update({f"charge_at_pH_{int(pH)}": 0 for pH in pH_values})
        cumulative_stats["aa_composition"] = {aa: 0 for aa in aa_order}
        cumulative_stats["aa_composition_percent"] = {aa: 0 for aa in aa_order}

        for seq in alternative_sequences:
            seq = seq.replace("X", "")  # Exclude 'X' from analysis
            protein_analyzer = ProtParam.ProteinAnalysis(seq)

            # Update cumulative stats
            for prop_name, func in properties:
                cumulative_stats[prop_name] += func(protein_analyzer, seq)

            # Charges at different pH levels
            for pH in pH_values:
                cumulative_stats[f"charge_at_pH_{int(pH)}"] += (
                    protein_analyzer.charge_at_pH(pH)
                )

            # Amino acid composition
            aa_composition = protein_analyzer.count_amino_acids()
            aa_composition_perc = protein_analyzer.get_amino_acids_percent()

            for aa in aa_order:
                cumulative_stats["aa_composition"][aa] += aa_composition.get(aa, 0)
                cumulative_stats["aa_composition_percent"][aa] += (
                    aa_composition_perc.get(aa, 0)
                )

        # Calculate the mean of results for all alternative sequences
        num_alt_seqs = len(alternative_sequences)
        mean_stats = {
            key: (value / num_alt_seqs)
            for key, value in cumulative_stats.items()
            if not isinstance(value, dict)
        }
        mean_stats["aa_composition"] = {
            aa: (count / num_alt_seqs)
            for aa, count in cumulative_stats["aa_composition"].items()
        }
        mean_stats["aa_composition_percent"] = {
            aa: (percent / num_alt_seqs)
            for aa, percent in cumulative_stats["aa_composition_percent"].items()
        }

        # Whole-protein composition entropy (non-windowed)
        mean_stats["composition_entropy"] = shannon_entropy(original_seq)

        # Flatten mean_stats and append to sequence_stats
        flat_mean_stats = flatten_stats(mean_stats)
        flat_mean_stats["id"] = seq_record.id
        sequence_stats.append(flat_mean_stats)

        # Compute windowed SD and autocorrelation
        # For ambiguous AA proteins, average across alternatives (same as mean_stats)
        cumulative_sd = None
        cumulative_autocorr = None
        for seq in alternative_sequences:
            seq = seq.replace("X", "")
            protein_analyzer = ProtParam.ProteinAnalysis(seq)
            profiles = compute_windowed_profiles(protein_analyzer, seq, 9)
            if profiles:
                sd_vals, autocorr_vals = summarize_profiles(profiles)
                if cumulative_sd is None:
                    cumulative_sd = {k: 0.0 for k in sd_vals}
                    cumulative_autocorr = {k: 0.0 for k in autocorr_vals}
                for k in sd_vals:
                    cumulative_sd[k] += sd_vals[k]
                for k in autocorr_vals:
                    cumulative_autocorr[k] += autocorr_vals[k]

        if cumulative_sd:
            num_alt = len(alternative_sequences)
            sd_row = {k: v / num_alt for k, v in cumulative_sd.items()}
            autocorr_row = {k: v / num_alt for k, v in cumulative_autocorr.items()}
        else:
            # Short protein or empty
            sd_row = {c: 0.0 for c in sd_columns if c != "id"}
            autocorr_row = {c: 0.0 for c in autocorr_columns if c != "id"}

        sd_row["id"] = seq_record.id
        autocorr_row["id"] = seq_record.id
        sequence_sds.append(sd_row)
        sequence_autocorrs.append(autocorr_row)

    # Create a dataframe of the sequence statistics
    sequence_stats_df = pd.DataFrame(sequence_stats)

    # Create DataFrames with specified column order
    df_counts = sequence_stats_df[aa_counts_columns]
    df_perc = sequence_stats_df[aa_perc_columns]
    df_summary = sequence_stats_df[summary_columns]

    # Save dataframes to CSV files
    gene_family_name = os.path.splitext(os.path.basename(msa_file))[0].split("_")[0]

    df_counts.to_csv(
        os.path.join(
            per_fam_basedir,
            "aa-counts",
            f"{gene_family_name}_aa_composition_counts.csv",
        ),
        index=False,
    )
    df_perc.to_csv(
        os.path.join(
            per_fam_basedir,
            "aa-proportions",
            f"{gene_family_name}_aa_composition_percentages.csv",
        ),
        index=False,
    )
    df_summary.to_csv(
        os.path.join(
            per_fam_basedir,
            "aa-physical-properties",
            f"{gene_family_name}_summary_statistics.csv",
        ),
        index=False,
    )

    # Write windowed SD and autocorrelation tables
    df_sd = pd.DataFrame(sequence_sds)[sd_columns]
    df_autocorr = pd.DataFrame(sequence_autocorrs)[autocorr_columns]

    df_sd.to_csv(
        os.path.join(
            per_fam_basedir,
            "aa-physical-property-sds",
            f"{gene_family_name}_summary_statistics_sd.csv",
        ),
        index=False,
    )
    df_autocorr.to_csv(
        os.path.join(
            per_fam_basedir,
            "aa-physical-property-autocorr",
            f"{gene_family_name}_summary_statistics_autocorr.csv",
        ),
        index=False,
    )

    # Note: The calculate_stats function and aggregation logic were designed
    # for batch processing multiple gene families. Since we now process one
    # file at a time, per-family CSVs are the final output.


if __name__ == "__main__":
    # Process the single MSA file
    process_msa(msa_file)
