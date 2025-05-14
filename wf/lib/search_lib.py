# General imports
import requests
import xml.etree.ElementTree as ET
from io import StringIO
import os
import time

# Biopython imports
from Bio import SeqIO

# Gets genbank files associated with a given gene name
def get_gbs(
    gene_name,
    organism,
    ncbi_key: str,
):
    print("Searching for", gene_name)
    # The base URL to call NIH nuccore API
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

    esearch_url = (base_url 
        + "esearch.fcgi" 
        + "?db=nuccore" 
        + "&api_key=" + ncbi_key)

    # We add onto the base URL using the gene of interest
    esearch_url += '&term=(' + gene_name + '[Gene Name] NOT PREDICTED [Title] AND srcdb_refseq[PROP]) AND "' + organism + '"[porgn] AND srcdb_refseq[PROP] AND biomol_mrna[PROP]'

    # Then, we call the API to get IDs for the transcripts
    response = requests.get(esearch_url, timeout=5)
    time.sleep(1)
    # Convert the response (in XML form) into an XML tree
    tree = ET.fromstring(response.text)

    # Check number of results returned
    if tree[0].text == "0": return []

    # Collect transcript IDs
    transcript_ids = [el.text for el in tree[3]]

    epost_url = (base_url 
        + "epost.fcgi" 
        + "?db=nuccore" 
        + "&api_key=" + ncbi_key 
        + "&id=" 
        + ",".join(transcript_ids))

    epost_res = requests.get(epost_url, timeout=5)
    time.sleep(1)
    epost_tree = ET.fromstring(epost_res.text)

    query_key = epost_tree[0].text
    web_env = epost_tree[1].text

    efetch_url = (base_url 
        + "efetch.fcgi" 
        + "?db=nuccore" 
        + "&api_key=" + ncbi_key 
        + "&query_key=" + query_key
        + "&WebEnv=" + web_env
        + "&rettype=gbwithparts" 
        + "&retmode=text"
        + "&retmax=100")

    efetch_res = requests.get(efetch_url, timeout=5)
    time.sleep(1)
    if efetch_res.status_code != 200:
        efetch_res.raise_for_status()
        raise RuntimeError(f"Request returned status code {efetch_res.status_code}")

    # Convert returned fastas into a String
    gbs = efetch_res.text

    # print("GBs for gene", gene_name, gbs)

    retrieved_records = list(SeqIO.parse(StringIO(gbs), "genbank"))

    print("Number of retrieved records for", gene_name, len(retrieved_records))

    # Return all records
    return [record for record in retrieved_records if not ("PREDICTED".lower() in record.description.lower())]

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