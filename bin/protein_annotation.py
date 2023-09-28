#!/usr/bin/env python3
import argparse
from bioservices import UniProt
from bioservices import EUtils
from concurrent.futures import ThreadPoolExecutor
from threading import Lock
import concurrent.futures
import pandas as pd
import os

# Columns (query fields) that we will use when accessing UniProt protein annotations
COLUMNS = ['organism_name', 'organism_id',  'accession', 'gene_names', 
    'gene_primary', 'gene_synonym', 'protein_name', 'length', 'mass', 
    'cc_mass_spectrometry', 'virus_hosts', 'organelle', 'cc_rna_editing',
    'go_p', 'go_c', 'go_f', 'go_id','absorption', 'ft_act_site', 'cc_activity_regulation',
    'ft_binding', 'cc_catalytic_activity', 'cc_cofactor','ft_dna_bind', 'ec', 
    'cc_function', 'kinetics', 'cc_pathway','ph_dependence', 'redox_potential', 
    'rhea', 'ft_site','temp_dependence','cc_interaction', 'cc_subunit',
    'xref_biogrid', 'xref_corum', 'xref_complexportal', 'xref_dip','xref_elm', 
    'xref_intact', 'xref_mint', 'xref_string','cc_allergen','cc_biotechnology', 
    'cc_disruption_phenotype','cc_disease','ft_mutagen', 'cc_pharmaceutical',
    'cc_toxic_dose','ft_intramem', 'cc_subcellular_location', 'ft_topo_dom', 
    'ft_transmem','ft_chain', 'ft_crosslnk', 'ft_disulfid', 'ft_carbohyd',
    'ft_init_met', 'ft_lipid', 'ft_mod_res', 'ft_peptide','cc_ptm', 'ft_propep', 
    'ft_signal', 'ft_transit','xref_carbonyldb', 'xref_depod', 'xref_glyconnect', 
    'xref_glygen','xref_metosite', 'xref_phosphositeplus', 'xref_swisspalm', 
    'xref_unicarbkb','xref_iptmnet','ft_coiled', 'ft_compbias', 'cc_domain', 
    'ft_domain','ft_motif', 'protein_families', 'ft_region', 'ft_repeat',
    'ft_zn_fing','xref_cdd', 'xref_disprot', 'xref_gene3d', 'xref_hamap',
    'xref_ideal', 'xref_interpro', 'xref_panther', 'xref_pirsf','xref_prints', 
    'xref_prosite', 'xref_pfam', 'xref_prodom','xref_sfld', 'xref_smart', 
    'xref_supfam', 'xref_ccds', 'xref_embl', 'xref_pir', 
    'xref_refseq','xref_alphafolddb', 'xref_bmrb', 'xref_pcddb', 'xref_pdb',
    'xref_pdbsum', 'xref_sasbdb', 'xref_smr', 'xref_brenda', 'xref_biocyc', 
    'xref_pathwaycommons', 'xref_plantreactome','xref_reactome', 'xref_sabio-rk', 
    'xref_signor', 'xref_signalink','xref_unipathway','xref_genetree', 'xref_hogenom', 
    'xref_inparanoid', 'xref_ko','xref_oma', 'xref_orthodb', 'xref_phylomedb', 
    'xref_treefam','xref_eggnog']

    # Now define the different annotation sets - we will pull out tables for each.
ANNOTATION_SETS = {
    'cogeqc': ['organism_name', 'organism_id', 'accession', 'xref_interpro', 'xref_oma'],
    'seq_meta': ['organism_name', 'organism_id',  'accession', 'gene_names', 'gene_primary', 'gene_synonym', 'protein_name', 'length', 'mass', 'cc_mass_spectrometry', 'virus_hosts', 'organelle', 'cc_rna_editing'],
    'seq_gos': ['organism_name', 'organism_id', 'accession', 'go_p', 'go_c', 'go_f', 'go_id'],
    'seq_funcs': ['organism_name', 'organism_id', 'accession', 'absorption', 'ft_act_site', 'cc_activity_regulation','ft_binding', 'cc_catalytic_activity', 'cc_cofactor','ft_dna_bind', 'ec', 'cc_function', 'kinetics', 'cc_pathway','ph_dependence', 'redox_potential', 'rhea', 'ft_site','temp_dependence'],
    'seq_inter': ['organism_name', 'organism_id', 'accession', 'cc_interaction', 'cc_subunit'],
    'seq_interdbs': ['organism_name', 'organism_id', 'accession', 'xref_biogrid', 'xref_corum', 'xref_complexportal', 'xref_dip','xref_elm', 'xref_intact', 'xref_mint', 'xref_string'],
    'seq_biotech': ['organism_name', 'organism_id', 'accession', 'cc_allergen','cc_biotechnology', 'cc_disruption_phenotype','cc_disease','ft_mutagen', 'cc_pharmaceutical','cc_toxic_dose'],
    'seq_localize': ['organism_name', 'organism_id', 'accession', 'ft_intramem', 'cc_subcellular_location', 'ft_topo_dom', 'ft_transmem'],
    'seq_ptm': ['organism_name', 'organism_id', 'accession', 'ft_chain', 'ft_crosslnk', 'ft_disulfid', 'ft_carbohyd','ft_init_met', 'ft_lipid', 'ft_mod_res', 'ft_peptide','cc_ptm', 'ft_propep', 'ft_signal', 'ft_transit'],
    'seq_xptmdb': ['organism_name', 'organism_id', 'accession', 'xref_carbonyldb', 'xref_depod', 'xref_glyconnect', 'xref_glygen','xref_metosite', 'xref_phosphositeplus', 'xref_swisspalm', 'xref_unicarbkb','xref_iptmnet'],
    'seq_domains': ['organism_name', 'organism_id', 'accession', 'ft_coiled', 'ft_compbias', 'cc_domain', 'ft_domain','ft_motif', 'protein_families', 'ft_region', 'ft_repeat','ft_zn_fing'],
    'seq_xfamdom': ['organism_name', 'organism_id', 'accession', 'xref_cdd', 'xref_disprot', 'xref_gene3d', 'xref_hamap','xref_ideal', 'xref_interpro', 'xref_panther', 'xref_pirsf','xref_prints', 'xref_prosite', 'xref_pfam', 'xref_prodom','xref_sfld', 'xref_smart', 'xref_supfam'],
    'seq_xseqdbs': ['organism_name', 'organism_id', 'accession', 'xref_ccds', 'xref_embl', 'xref_pir', 'xref_refseq'],
    'seq_x3ddbs': ['organism_name', 'organism_id', 'accession', 'xref_alphafolddb', 'xref_bmrb', 'xref_pcddb', 'xref_pdb','xref_pdbsum', 'xref_sasbdb', 'xref_smr'],
    'seq_xenzpath': ['organism_name', 'organism_id', 'accession', 'organism_name', 'organism_id', 'accession', 'xref_brenda', 'xref_biocyc', 'xref_pathwaycommons', 'xref_plantreactome','xref_reactome', 'xref_sabio-rk', 'xref_signor', 'xref_signalink','xref_unipathway'],
    'seq_xortho': ['organism_name', 'organism_id', 'accession', 'xref_genetree', 'xref_hogenom', 'xref_inparanoid', 'xref_ko','xref_oma', 'xref_orthodb', 'xref_phylomedb', 'xref_treefam','xref_eggnog']
}

# Suffixes for each annotation filename
ANNOT_NAMES = ['_cogeqc_annotations.tsv', '_prot_metadat.tsv', '_prot_gene_ontologies.tsv', 
    '_prot_functions.tsv', '_prot_interactions.tsv', '_prot_interactions_xref.tsv',
    '_prot_biotech_annots.tsv', '_prot_localization.tsv', '_prot_post_trans_mods.tsv', 
    '_prot_post_trans_mods_xref.tsv', '_prot_fams_domains.tsv', '_prot_fams_domains_xref.tsv',
    '_prot_seq_dbs_xref.tsv', '_prot_3d_dbs_xref.tsv', '_prot_enzyme_paths_xref.tsv', 
    '_prot_orthology_dbs.tsv']
        
def fetch_batch(batch_accessions, organism_name, columns):
    uniprot = UniProt()
    query = " OR ".join([f"accession:{acc}" for acc in batch_accessions])
    result = uniprot.search(query, frmt="tsv", columns=",".join(columns), limit=None)

    annotations = []
    if result and isinstance(result, str):
        lines = result.strip().split('\n')[1:]  # Ignore the header row
        for line in lines:
            fields = line.split('\t')
            annotations.append(fields)
    return annotations

def get_annotations(organism_name, input_file, columns, num_workers=None):
    uniprot = UniProt()
    annotations = []

    with open(input_file, 'r') as file:
        accessions = file.read().splitlines()

    # Split the accessions into batches of 100
    batch_size = 100
    batches = [accessions[i:i + batch_size] for i in range(0, len(accessions), batch_size)]
    
    if num_workers is None:
        num_workers = min(os.cpu_count(), len(batches))
    
    # Print out some info:
    print(f"Pulling down annotations for {organism_name}")
    
    # Create a lock for thread-safe list updates
    annotations_lock = Lock()

    # Run the loop in parallel using ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        future_results = [executor.submit(fetch_batch, batch_accessions, organism_name, columns) for batch_accessions in batches]
        for i, future in enumerate(concurrent.futures.as_completed(future_results)):
            batch_annotations = future.result()
            with annotations_lock:
                annotations.extend(batch_annotations)
            print(f"Completed batch {i + 1} of {len(batches)}")

    return annotations

def main():
    parser = argparse.ArgumentParser(description='Retrieve UniProt annotations for a given species and list of protein accessions.')
    parser.add_argument('spp', help='Species to be queried (in Genus_species format)')
    parser.add_argument('ids', help='File containing the UniProt protein accessions')

    args = parser.parse_args()

    organism_name = args.spp
    input_file = args.ids

    # Begin by pulling down all annotations that we may want
    annotations = get_annotations(organism_name, input_file, COLUMNS)
    annots_df = pd.DataFrame(annotations, columns=COLUMNS)
    
    # Create a dictionary of data frames for each annotation set. 
    annot_df_dict = {key: annots_df[cols] for key, cols in ANNOTATION_SETS.items()}
    
    # Save each DataFrame subset to a TSV file with the appropriate prefix and suffix
    for idx, (key, df) in enumerate(annot_df_dict.items()):
        output_filename = f"{organism_name}{ANNOT_NAMES[idx]}"
        df.to_csv(output_filename, sep='\t', index=False)

if __name__ == "__main__":
    main()
    
