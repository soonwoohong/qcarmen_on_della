# General imports
# this requests crawling keeps failing on Della. In my opinion, it accesses the source too fast & too many.
# replace this with Bio.Entrez
#import requests
import xml.etree.ElementTree as ET
from io import StringIO
import os
import time

# Biopython imports
from Bio import SeqIO, Entrez

# Gets genbank files associated with a given gene name
def get_gbs(
    gene_name,
    organism,
    ncbi_key: str,
    user_email: str
):
    print("Searching for", gene_name)
    Entrez.api_key = ncbi_key
    Entrez.email = user_email
    # The base query
    # Build your query string exactly as before
    query = (
        f'({gene_name}[Gene Name] NOT PREDICTED[Title] AND srcdb_refseq[PROP]) '
        f'AND "{organism}"[porgn] AND srcdb_refseq[PROP] AND biomol_mrna[PROP]'
    )
    # 1) ESearch: get up to `retmax` IDs
    with Entrez.esearch(db="nuccore", term=query) as es:
        search_results = Entrez.read(es)

    id_list = search_results.get("IdList", [])
    if not id_list:
        print("  → no hits.")
        return []

    # 2) EFetch: retrieve GenBank (gbwithparts) for those IDs
    ids = ",".join(id_list)
    with Entrez.efetch(
            db="nuccore",
            id=ids,
            rettype="gbwithparts",
            retmode="text"
    ) as fetch_handle:
        records = list(SeqIO.parse(fetch_handle, "genbank"))

    # Filter out any predicted entries just in case
    filtered = [r for r in records if "PREDICTED" not in r.description.upper()]
    print(f"  → retrieved {len(filtered)} records (of {len(records)})")
    return filtered

# Save a set of Genbank files for a given gene provided a root directory for gbs
def save_gbs(gene_name, gb_list, gbs_root):
    # Root dir
    gene_dir = gbs_root + gene_name + "/"

    # Folder name
    if not os.path.exists(gene_dir):
        os.makedirs(gene_dir)

    # Loop through genbank files in each sublist
    for gb_file in gb_list:
        output_file = open(gene_dir + gb_file.id + ".gb", 'w')
        SeqIO.write(gb_file, output_file, 'genbank')
        output_file.close()

def read_gbs(gb_path):
    """
    Reads Genbank files from a given path and returns a list of Genbnk records.
    """
    gbs = []

    mod_path = gb_path + "/" if not gb_path.endswith('/') else gb_path
    for file in os.listdir(gb_path):
        gb_indiv = SeqIO.read(mod_path + file, "genbank")
        gbs.append(gb_indiv)

    return gbs