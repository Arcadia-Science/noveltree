import click
import pandas as pd
import pathlib


@click.group()
def cli():
    pass


@cli.command()
@click.option(
    "--output-filepath",
    required=True,
    type=click.Path(
        exists=False, file_okay=True, dir_okay=False, path_type=pathlib.Path
    ),
)
def download(output_filepath: pathlib.Path) -> None:
    """
    Download the HGNC dataset from the HGNC website.
    We use the `hgnc_complete_set.txt` file because some UniProt IDs in the OS dataset
    are not present in the `gene_with_protein_product.txt` file.
    """
    url = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"

    hgnc_dataset = pd.read_csv(url, sep="\t")

    hgnc_dataset.dropna(
        subset=["hgnc_id", "symbol", "uniprot_ids"], how="any", inplace=True
    )

    # Check that the primary symbols are unique.
    if hgnc_dataset["symbol"].nunique() != hgnc_dataset.shape[0]:
        raise ValueError("The primary symbols are not unique.")

    hgnc_dataset.to_csv(output_filepath, index=False, sep="\t")


@cli.command()
@click.option(
    "--input-filepath",
    required=True,
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False, path_type=pathlib.Path
    ),
)
@click.option(
    "--output-filepath",
    required=True,
    type=click.Path(
        exists=False, file_okay=True, dir_okay=False, path_type=pathlib.Path
    ),
)
def process(input_filepath: pathlib.Path, output_filepath: pathlib.Path) -> None:
    """
    Process the HGNC dataset into the form we need for the OS portal.

    Outputs a TSV file with one row for each unique combination of UniProt ID, primary gene symbol,
    and alias gene symbol in the HGNC dataset.

    The output file has the following columns:
    - hgnc_id: The HGNC ID.
    - uniprot_id: The UniProt ID associated with the HGNC entry.
    - symbol: The gene symbol (either primary or alias).
    - primary_symbol: The primary gene symbol (same as the `symbol` column for primary symbols).
    - symbol_is_alias: A boolean indicating whether the value in `symbol` is an alias.
    - name: The gene name.
    """

    # The separator used in the HGNC dataset for list-like fields.
    sep = "|"

    hgnc_dataset = pd.read_csv(input_filepath, sep="\t")

    # Filter out locus types that we know are not relevant to the OS dataset.
    excluded_locus_types = ["RNA, long non-coding", "RNA, ribosomal"]
    hgnc_dataset["locus_type"] = hgnc_dataset["locus_type"].str.strip().str.lower()
    hgnc_dataset = hgnc_dataset[
        ~hgnc_dataset.locus_type.isin([s.strip().lower() for s in excluded_locus_types])
    ].copy()

    # Drop all columns except the ones we need.
    hgnc_dataset = hgnc_dataset[
        ["hgnc_id", "symbol", "alias_symbol", "uniprot_ids", "name"]
    ].copy()

    # A few gene symbols have multiple UniProt IDs, so we expand them into separate rows.
    hgnc_dataset["uniprot_id"] = hgnc_dataset["uniprot_ids"].str.split(sep)
    hgnc_dataset = hgnc_dataset.explode("uniprot_id").drop(columns=["uniprot_ids"])

    # Create a dataframe of alias symbols, with one row per alias symbol.
    hgnc_alias_symbols = hgnc_dataset.loc[hgnc_dataset.alias_symbol.notnull()].copy()
    hgnc_alias_symbols["alias_symbol"] = hgnc_alias_symbols["alias_symbol"].str.split(
        sep
    )
    hgnc_alias_symbols = hgnc_alias_symbols.explode("alias_symbol").dropna(
        axis=0, how="any"
    )

    # Rename columns to allow us to concatenate the primary and alias symbols.
    hgnc_alias_symbols = hgnc_alias_symbols.rename(
        columns={"symbol": "primary_symbol"}
    ).rename(columns={"alias_symbol": "symbol"})

    # Create a dataframe of primary symbols with the same column names as the alias symbol dataframe.
    hgnc_primary_symbols = hgnc_dataset.copy()
    hgnc_primary_symbols["primary_symbol"] = hgnc_primary_symbols["symbol"].copy()

    # Concatenate the primary and alias symbol datasets.
    final_columns = [
        "hgnc_id",
        "uniprot_id",
        "symbol",
        "primary_symbol",
        "name",
    ]
    final_dataset = pd.concat(
        (hgnc_primary_symbols[final_columns], hgnc_alias_symbols[final_columns])
    )

    # For extra clarity, add a column to indicate whether the symbol is an alias.
    final_dataset["symbol_is_alias"] = (
        final_dataset["symbol"] != final_dataset["primary_symbol"]
    )

    final_dataset.sort_values(by="symbol").to_csv(
        output_filepath, index=False, sep="\t"
    )


if __name__ == "__main__":
    cli()
