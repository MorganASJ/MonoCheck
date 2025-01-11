import argparse
from Bio import Entrez
import pandas as pd
import os
from ete3 import Tree

Entrez.email = "Morgan.Jones@bristol.ac.uk"

def request_taxonomy_data(taxa_list):
    taxonomy_data = []

    # For each taxon make a request to Entrez NCBI Taxonomy
    for taxon in taxa_list:
        try:
            # Search for the taxon in NCBI Taxonomy
            search = Entrez.esearch(db="taxonomy", term=taxon)
            search_result = Entrez.read(search)
            search.close()

            if not search_result["IdList"]:
                print(f"Erorr: Taxon not found: {taxon}")
                continue

            # Fetch taxonomy details
            tax_id = search_result["IdList"][0]
            fetch = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
            fetch_result = Entrez.read(fetch)
            fetch.close()

            #Retreive taxonomy details
            tax_data = fetch_result[0]
            lineage_dict = {rank["Rank"]: rank["ScientificName"] for rank in tax_data["LineageEx"]}
            lineage_dict["Species"] = tax_data["ScientificName"]
            lineage_dict["Taxon"] = taxon
            taxonomy_data.append(lineage_dict)

        except Exception as e:
            print(f"Error fetching data for {taxon}: {e}")

    return taxonomy_data

def fetch_taxonomy(taxa, taxonomy_file="taxonomy.tax"):

    taxonomy_data = {}
    # parts = len(taxa[0].split("_"))

    taxa_list = ["_".join((taxon[0],taxon[1])) for taxon in [taxon.split("_") for taxon in taxa]]
    print(taxa_list)

    # First check if a file exists
    if os.path.exists(taxonomy_file):
        print(f"Using existing taxonomy file: {taxonomy_file}")
        try:
            records = pd.read_csv(taxonomy_file).to_dict(orient='records')       
        except:
            raise Exception("Failure to read Taxonomy File")
    
    else:

        # return 0
        records = request_taxonomy_data(taxa_list)
        # Create a DataFrame and save to CSV
        df = pd.DataFrame(records).fillna("N/A")
        df.to_csv(taxonomy_file, index=False)
        print(f"Taxonomy data saved to {taxonomy_file}")

    for rec in records:
        taxonomy_data[rec['Taxon']] = rec

    return taxonomy_data

# Root tree based on taxa listed in file
def root_tree(tree, outgroup_file):
    with open(outgroup_file) as f:
        outgroup_taxa = f.read().strip().split("\n")

    print(outgroup_taxa) # Works

    # if set of taxa, find commmon ancestor node, else root on leaf node
    if len(outgroup_taxa) > 1:
        tree.set_outgroup(tree.get_common_ancestor(outgroup_taxa))
    else:
        tree.set_outgroup(outgroup_taxa[0])

    # return transformed tree object
    return tree

# TODO
def is_monophyletic(tree, group, taxonomy):
    # ete3 is is_monophyletic?
    return True

# Entry point
def main():

    # Parse arguments python monocheck [tree file] -group [taxonomic group] -outgroup [outgroup file OR outgroup taxon]
    parser = argparse.ArgumentParser(description="Check monophyly of a group in a phylogenetic tree.")
    parser.add_argument("tree", help="Newick file containing the tree")
    parser.add_argument("-group", required=True, help="Group to check for monophyly")
    # parser.add_argument("-family", required=False, help="Group to check for monophyly")
    # parser.add_argument("-order", required=False, help="Group to check for monophyly")
    # parser.add_argument("-genus", required=False, help="Group to check for monophyly")
    # parser.add_argument("-species", required=False, help="Group to check for monophyly")
    # parser.add_argument("-phylum", required=False, help="Group to check for monophyly")
    parser.add_argument("-outgroup", required=True, help="File containing outgroup taxa (one per line)")
    args = parser.parse_args()

    # Load tree
    tree = Tree(args.tree, format=1)

    # print(tree) # Works

    # Root the tree
    tree = root_tree(tree, args.outgroup)

    # print(tree) # Works

    # Load or generate taxonomy
    taxonomy_file = "taxonomy.tax"
    taxonomy = fetch_taxonomy([leaf.name for leaf in tree.get_leaves()], taxonomy_file)

    # print(taxonomy)

        # Check monophyly
    result = is_monophyletic(tree, args.group, taxonomy)
    print(f"Is {args.group} monophyletic? {'Yes' if result else 'No'}")

if __name__ == "__main__":
    main()

    