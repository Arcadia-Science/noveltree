import os
import shutil

import click
import pandas as pd
import pathlib
import tqdm

# The number of orthogroups in the OS dataset on Zenodo.
NUM_ORTHOGROUPS = 9260

# Organism names that should be excluded from the OS dataset shown in the portal.
# (These correspond to proteomes that were mislabeled.)
ORGANISM_NAMES_TO_EXCLUDE = [
    "Hypsibius-dujardini",
    "Porphyra-yezoensis",
]


@click.group()
def cli():
    pass


@cli.command()
@click.option(
    "--dataset-dirpath",
    required=True,
    type=click.Path(
        file_okay=False, dir_okay=True, exists=True, path_type=pathlib.Path
    ),
    help="Path to the directory containing the OS dataset downloaded from Zenodo.",
)
@click.option(
    "--output-filepath",
    required=True,
    type=click.Path(file_okay=True, dir_okay=False, path_type=pathlib.Path),
)
def concat(dataset_dirpath, output_filepath):
    """
    Concatenate the OS dataset files into a single file.
    There is one TSV file per orthogroup ("OG"), and there are 9,260 orthogroups in total.
    """

    if output_filepath.exists():
        click.echo(f"Removing existing output file '{output_filepath}'")
        os.remove(output_filepath)

    tsv_dirpath = (
        dataset_dirpath
        / "gf-aa-multivar-distances"
        / "final_protein_pair_summary_tables"
    )

    tsv_filepaths = list(tsv_dirpath.glob("*.tsv"))
    if len(tsv_filepaths) != NUM_ORTHOGROUPS:
        raise ValueError(
            f"Expected {NUM_ORTHOGROUPS} TSV files in the directory, but found {len(tsv_filepaths)}."
        )

    df = pd.concat(
        [
            pd.read_csv(filepath, sep="\t")
            for filepath in tqdm.tqdm(tsv_filepaths, desc="Reading TSV files")
        ]
    )
    df.to_csv(output_filepath, sep="\t", index=False)


def _merge_hgnc_gene_symbols(os_dataset, hgnc_dataset_filepath):
    """
    Merge the organism-selection dataset with the HGNC dataset to append gene symbols.

    os_dataset : pd.DataFrame
        The concatenated organism-selection dataset.
    hgnc_dataset_filepath : pathlib.Path
        Path to the processed HGNC dataset produced by the `hgnc_dataset.py` script.

    Returns the `os_dataset` dataframe with the HGNC gene symbols in the `hgnc_gene_symbol` column.
    When there is more than one gene symbol for a given row, the symbols are aggregated
    into a comma-separated list.
    """

    hgnc_dataset = pd.read_csv(hgnc_dataset_filepath, sep="\t")

    # We only want to merge the primary gene symbols, not the aliases.
    hgnc_dataset = hgnc_dataset.loc[
        hgnc_dataset.symbol_is_alias.notna() & (~hgnc_dataset.symbol_is_alias)
    ].copy()

    # Check that we still have data after filtering.
    if len(hgnc_dataset) == 0:
        raise ValueError(
            "No non-alias gene symbols found in HGNC dataset after filtering"
        )

    hgnc_dataset.rename(columns={"symbol": "hgnc_gene_symbol"}, inplace=True)

    # Some UniProt IDs map to more than one gene symbol,
    # so we aggregate them into a comma-separated list.
    hgnc_dataset = (
        hgnc_dataset.groupby("uniprot_id")["hgnc_gene_symbol"]
        .agg(lambda s: ",".join(sorted(s.dropna().unique())))
        .reset_index()
    )

    # Check that uniprot_id is unique, since we're about to merge on it.
    if hgnc_dataset.uniprot_id.nunique() != hgnc_dataset.shape[0]:
        raise ValueError("Uniprot IDs are not unique in the HGNC dataset.")

    os_dataset = pd.merge(
        os_dataset,
        hgnc_dataset,
        left_on="ref_protein",
        right_on="uniprot_id",
        how="left",
    )

    os_dataset.drop(axis=1, labels=["uniprot_id"], inplace=True)

    num_ref_proteins_missing_gene_symbols = os_dataset.loc[
        os_dataset.hgnc_gene_symbol.isnull()
    ].ref_protein.nunique()

    if num_ref_proteins_missing_gene_symbols > 0:
        click.echo(
            f"Warning: {num_ref_proteins_missing_gene_symbols} of {os_dataset.ref_protein.nunique()} "
            "reference proteins are missing gene symbols."
        )

    return os_dataset


def _append_organism_percentile_and_rank(os_dataset):
    """
    Calculate the rank of the reference proteins within each non-reference species
    and append the results to the dataframe.

    The `trait_dist` column is the distance between the reference and non-reference proteins.
    It is not comparable across different gene families, so we calculate the percentile
    of the `trait_dist` values within each gene family, then calculate the rank from the percentile.
    """
    os_dataset = os_dataset.copy()
    percentile_column = "within_gene_family_percentile"

    os_dataset[percentile_column] = os_dataset.groupby(["gene_family"])[
        "trait_dist"
    ].rank(pct=True)

    os_dataset["organism_rank"] = os_dataset.groupby(["nonref_species"])[
        percentile_column
    ].rank(method="min", ascending=True)

    return os_dataset


def _append_portfolio_rank(os_dataset):
    """
    Calculate the rank of the non-reference proteins for each reference protein
    (that is, within gene families) and append the results to the dataframe.

    This rank is calculated within each gene family, so we can use the `trait_dist` column directly.
    """
    os_dataset = os_dataset.copy()
    os_dataset["portfolio_rank"] = os_dataset.groupby(["ref_protein"])[
        "trait_dist"
    ].rank(method="min", ascending=True)
    return os_dataset


def _append_organism_homolog_rank_and_count(os_dataset):
    """
    Calculate the count and rank of the non-reference proteins for each reference protein
    and non-reference species, and append the results to the dataframe.

    This rank is calculated within each gene family, so we can use the `trait_dist` column directly.
    """
    os_dataset = os_dataset.copy()

    os_dataset["orthogroup_homolog_rank"] = os_dataset.groupby(
        ["ref_protein", "nonref_species"]
    )["trait_dist"].rank(method="min", ascending=True)

    os_dataset["orthogroup_homolog_count"] = os_dataset.groupby(
        ["ref_protein", "nonref_species"]
    )["trait_dist"].transform("count")

    return os_dataset


@cli.command()
@click.option(
    "--dataset-filepath",
    required=True,
    type=click.Path(
        file_okay=True, dir_okay=False, exists=True, path_type=pathlib.Path
    ),
    help="Path to the concatenated OS dataset TSV file.",
)
@click.option(
    "--hgnc-dataset-filepath",
    required=True,
    type=click.Path(
        file_okay=True, dir_okay=False, exists=True, path_type=pathlib.Path
    ),
    help="Path to the processed HGNC dataset produced by the `hgnc_dataset.py` script.",
)
@click.option(
    "--output-dirpath",
    required=True,
    type=click.Path(file_okay=False, dir_okay=True, path_type=pathlib.Path),
)
def process(dataset_filepath, hgnc_dataset_filepath, output_dirpath):
    """
    Process and clean up the concatenated organism-selection dataset,
    merge it with a UniProt ID mapping dataset to append gene symbols,
    then split it into per-organism and per-reference-protein files.
    """

    nonref_species_dirpath = output_dirpath / "per-nonref-species"
    ref_protein_dirpath = output_dirpath / "per-ref-protein"
    orthogroup_dirpath = output_dirpath / "per-orthogroup"

    nonref_species_dirpath.mkdir(parents=True, exist_ok=True)
    ref_protein_dirpath.mkdir(parents=True, exist_ok=True)
    orthogroup_dirpath.mkdir(parents=True, exist_ok=True)

    os_dataset = pd.read_csv(dataset_filepath, sep="\t")

    # Check for the expected number of orthogroups (or "gene families").
    if os_dataset.gene_family.nunique() != NUM_ORTHOGROUPS:
        raise ValueError(
            f"Expected {NUM_ORTHOGROUPS} unique gene families, "
            f"but found {os_dataset.gene_family.nunique()} in the concatenated dataset."
        )

    # The `trait_dist` column must be numeric.
    os_dataset["trait_dist"] = pd.to_numeric(os_dataset["trait_dist"], errors="raise")

    # Drop rows with missing values (there shouldn't be any).
    num_rows_before = os_dataset.shape[0]
    os_dataset.dropna(axis=0, how="any", inplace=True)
    num_rows_after = os_dataset.shape[0]
    click.echo(
        f"{num_rows_before - num_rows_after} rows contained missing values and were dropped."
    )

    # Drop rows corresponding to the excluded organisms.
    os_dataset = os_dataset.loc[
        ~os_dataset.nonref_species.str.lower().isin(
            [name.lower() for name in ORGANISM_NAMES_TO_EXCLUDE]
        )
    ].copy()

    # Drop unneeded columns: the `phylo_dist` column doesn't need to be shown to the user,
    # and we recalculate the ranks below, so we don't need the `rank_trait_dist` column.
    os_dataset.drop(axis=1, labels=["phylo_dist", "rank_trait_dist"], inplace=True)

    # Append the ranks to the dataset.
    os_dataset = _append_organism_percentile_and_rank(os_dataset)
    os_dataset = _append_portfolio_rank(os_dataset)
    os_dataset = _append_organism_homolog_rank_and_count(os_dataset)

    # Merge the organism-selection dataset with the HGNC dataset to append gene symbols.
    num_rows_before = os_dataset.shape[0]
    os_dataset = _merge_hgnc_gene_symbols(os_dataset, hgnc_dataset_filepath)
    num_rows_after = os_dataset.shape[0]
    click.echo(
        f"{num_rows_before - num_rows_after} rows were dropped when merging with HGNC dataset."
    )

    click.echo(f"Writing per-organism files to '{output_dirpath}'")
    for species, group in tqdm.tqdm(os_dataset.groupby("nonref_species")):
        filepath = nonref_species_dirpath / f"{species}.tsv"
        group.sort_values(by="organism_rank").to_csv(filepath, sep="\t", index=False)

    click.echo(f"Writing per-reference-protein files to '{output_dirpath}'")
    for ref_protein, group in tqdm.tqdm(os_dataset.groupby("ref_protein")):
        filepath = ref_protein_dirpath / f"{ref_protein}.tsv"
        group.sort_values(by="portfolio_rank").to_csv(filepath, sep="\t", index=False)

    click.echo(f"Writing per-orthogroup files to '{output_dirpath}'")
    for _, group in tqdm.tqdm(os_dataset.groupby("gene_family")):
        # We write the same per-orthogroup file for each reference protein in the orthogroup,
        # rather than a single file for the entire orthogroup.
        # Although this results in some files being duplicated, it simplifies the frontend code
        # by eliminating the need to map reference proteins to orthogroups.
        unique_ref_proteins = list(group.ref_protein.unique())
        filepath = orthogroup_dirpath / f"{unique_ref_proteins[0]}.tsv"
        group.sort_values(by="portfolio_rank").to_csv(filepath, sep="\t", index=False)
        for ref_protein in unique_ref_proteins[1:]:
            shutil.copy(filepath, orthogroup_dirpath / f"{ref_protein}.tsv")


@cli.command()
@click.option(
    "--dataset-filepath",
    required=True,
    type=click.Path(
        file_okay=True, dir_okay=False, exists=True, path_type=pathlib.Path
    ),
    help="Path to the concatenated OS dataset TSV file.",
)
@click.option(
    "--output-dirpath",
    required=True,
    type=click.Path(file_okay=False, dir_okay=True, path_type=pathlib.Path),
)
def sanity_check(dataset_filepath, output_dirpath):
    """
    Check that the final output files are correct, relative to the original concatenated OS dataset.
    """
    os_dataset = pd.read_csv(dataset_filepath, sep="\t")

    # Drop the excluded organisms.
    os_dataset = os_dataset.loc[
        ~os_dataset.nonref_species.str.lower().isin(
            [name.lower() for name in ORGANISM_NAMES_TO_EXCLUDE]
        )
    ].copy()

    num_os_dataset_rows = os_dataset.shape[0]

    nonref_species_filepaths = list(
        (output_dirpath / "per-nonref-species").glob("*.tsv")
    )
    ref_protein_filepaths = list((output_dirpath / "per-ref-protein").glob("*.tsv"))
    orthogroup_filepaths = list((output_dirpath / "per-orthogroup").glob("*.tsv"))

    if len(nonref_species_filepaths) != len(os_dataset.nonref_species.unique()):
        click.echo(
            f"Expected {len(os_dataset.nonref_species.unique())} non-reference species files, "
            f"but found {len(nonref_species_filepaths)} files."
        )
    else:
        click.echo("The number of non-reference-species files is correct")

    if len(ref_protein_filepaths) != len(os_dataset.ref_protein.unique()):
        click.echo(
            f"Expected {len(os_dataset.ref_protein.unique())} reference protein files, "
            f"but found {len(ref_protein_filepaths)} files."
        )
    else:
        click.echo("The number of reference-protein files is correct")

    # There should be one orthogroup file per reference protein
    # (this involves some duplication, which is by design).
    if len(orthogroup_filepaths) != len(os_dataset.ref_protein.unique()):
        click.echo(
            f"Expected {len(os_dataset.ref_protein.unique())} orthogroup files, "
            f"but found {len(orthogroup_filepaths)} files."
        )
    else:
        click.echo("The number of orthogroup files is correct")

    click.echo(
        "Checking that the total number of rows in the non-reference-species files is correct"
    )
    num_rows_in_nonref_species_files = 0
    for filepath in tqdm.tqdm(nonref_species_filepaths):
        df = pd.read_csv(filepath, sep="\t")
        num_rows_in_nonref_species_files += df.shape[0]
    if num_rows_in_nonref_species_files != num_os_dataset_rows:
        click.echo(
            f"Expected {num_os_dataset_rows} rows in the non-reference species files, "
            f"but found {num_rows_in_nonref_species_files} rows."
        )
    else:
        click.echo("The number of rows in the non-reference-species files is correct")

    click.echo(
        "Checking that the total number of rows in the reference-protein files is correct"
    )
    num_rows_in_ref_protein_files = 0
    for filepath in tqdm.tqdm(ref_protein_filepaths):
        df = pd.read_csv(filepath, sep="\t")
        num_rows_in_ref_protein_files += df.shape[0]
    if num_rows_in_ref_protein_files != num_os_dataset_rows:
        click.echo(
            f"Expected {num_os_dataset_rows} rows in the reference protein files, "
            f"but found {num_rows_in_ref_protein_files} rows."
        )
    else:
        click.echo("The number of rows in the reference-protein files is correct")


if __name__ == "__main__":
    cli()
