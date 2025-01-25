#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import requests
import argparse
import os


def search_genes(gene_name, organism="Homo sapiens"):
    """Search for gene IDs using NCBI's eSearch."""
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    parameters = {
        "db": "gene",
        "term": f'"{gene_name}"[gene] AND "{organism}"[organism]',  # Corrected query format
        "retmode": "json"
    }

    try:
        response = requests.get(base_url, params=parameters)
        response.raise_for_status()  # Check if the request was successful
        data = response.json()

        # Return the list of gene IDs from the response
        return data.get("esearchresult", {}).get("idlist", [])
    except requests.exceptions.RequestException as e:
        print(f"Error searching for {gene_name}: {e}")
        return []


def fetch_genes(gene_id):
    """Fetch gene data (FASTA sequence) using NCBI's eFetch."""
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    parameters = {
        "db": "nuccore",
        "id": gene_id,
        "rettype": "fasta",
        "retmode": "text"
    }

    try:
        response = requests.get(base_url, params=parameters)
        response.raise_for_status()  # Check if the request was successful
        return response.text
    except requests.exceptions.RequestException as e:
        print(f"Error fetching data for Gene ID {gene_id}: {e}")
        return None


def fetch_genes_from_file(input_file, organism, output_dir):
    """Fetch gene sequences for genes listed in the input file and save to output directory."""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with open(input_file, "r") as file:
        genes = file.read().splitlines()

    for gene in genes:
        print(f"Searching for {gene}...")
        gene_ids = search_genes(gene, organism)

        if not gene_ids:
            print(f"No gene found for {gene}.")
            continue

        for gene_id in gene_ids:
            print(f"Fetching sequence for Gene ID {gene_id}...")
            sequence = fetch_genes(gene_id)

            if sequence:
                output_file = os.path.join(output_dir, f"{gene}_{gene_id}.fasta")
                with open(output_file, "w") as f:
                    f.write(sequence)
                print(f"Saved to {output_file}")
            else:
                print(f"Failed to fetch sequence for Gene ID {gene_id}.")


def main():
    """Main function to handle command-line arguments and execute fetching."""
    parser = argparse.ArgumentParser(description="Fetch gene sequences from NCBI")
    parser.add_argument("input", help="Input file containing gene names")
    parser.add_argument("--organism", default="Homo sapiens", help="Organism name")
    parser.add_argument("--output_dir", default="gene_sequences", help="Directory to save gene sequences.")

    args = parser.parse_args()

    # Fetch genes based on input file and save to output directory
    fetch_genes_from_file(args.input, args.organism, args.output_dir)


if __name__ == "__main__":
    main()
