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

def fetch_taxonomy(taxa_list, taxonomy_file="taxonomy.tax"):

    taxonomy_data = {}
    # parts = len(taxa[0].split("_"))

    
    print(taxa_list)

    # First check if a file exists
    if os.path.exists(taxonomy_file):
        print(f"Using existing taxonomy file: {taxonomy_file}")
        try:
            # records = pd.read_csv(taxonomy_file).to_dict(orient='records')
            taxonomy_data = pd.read_csv(taxonomy_file)   
   
        except:
            raise Exception("Failure to read Taxonomy File")
    
    else:

        # return 0
        records = request_taxonomy_data(taxa_list)
        # Create a DataFrame and save to CSV
        taxonomy_data = pd.DataFrame(records).fillna("N/A")
        taxonomy_data.to_csv(taxonomy_file, index=False)
        print(f"Taxonomy data saved to {taxonomy_file}")

        

    # for rec in records:
        # taxonomy_data[rec['Taxon']] = rec

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
def is_monophyletic(tree, group, taxonomy, leaf_taxon):

    ranks = taxonomy.columns[taxonomy.isin([group]).any()]

    if ranks.empty:
        raise Exception("Taxonomic Group not Present in Tree")
    elif len(ranks) > 1:
        raise Exception("Taxonomic group appears under multiple ranks")
    
    print(ranks)

    group_records = taxonomy[["Taxon",ranks[0]]]

    print(group_records)

    mono_taxa = group_records[taxonomy[ranks[0]]==group]

    mono_leaves = [leaf for leaf in leaf_taxon if leaf.rsplit('_', 1)[0] in set(mono_taxa['Taxon'])]

    print(mono_leaves)

    isMono = tree.check_monophyly(values=mono_leaves, target_attr="name")

    print(isMono)

    return isMono

    return 1
    # ete3 is is_monophyletic
    # for taxon, rec in taxonomy:
        # if group in rec.values:     ?
    
    # leaf_names = 
    # tree.check_monophyly(values=["a", "e", "i", "o", "u"], target_attr="name")

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

    print(tree) # Works
    
    # Load or generate taxonomy
    taxonomy_file = "taxonomy.tax"
    leaf_names = [leaf.name for leaf in tree.get_leaves()]
    taxa_list = ["_".join((taxon[0],taxon[1])) for taxon in [taxon.split("_") for taxon in leaf_names]]
    taxonomy = fetch_taxonomy(taxa_list, taxonomy_file)

    print(taxonomy)

    print(taxonomy['Taxon'])

        # Check monophyly
    result = is_monophyletic(tree, args.group, taxonomy, leaf_names)
    return 1
    print(f"Is {args.group} monophyletic? {'Yes' if result else 'No'}")

if __name__ == "__main__":
    main()

    