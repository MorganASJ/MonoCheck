import argparse
from Bio import Entrez
import pandas as pd
import os
from ete3 import Tree
import sys, time, threading

Entrez.email = "Morgan.Jones@bristol.ac.uk"

#Spinner while taxonomy calculates
def animated_loading(message):
    chars = "/—\|" 
    for char in chars:
        sys.stdout.write('\r' + char + message +'\r')
        time.sleep(.1)
        sys.stdout.flush()

def request_taxonomy_data(taxa_list):
    taxonomy_data = []

    start_time = time.time()

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

    print("--- %s seconds ---" % (time.time() - start_time))

    return taxonomy_data

def request_taxonomy_data1(taxa_list):
    taxonomy_data = []

    start_time = time.time()

    sterm = " OR ".join(taxa_list)

    # Search for the taxon in NCBI Taxonomy
    search = Entrez.esearch(db="taxonomy", term=sterm)
    search_result = Entrez.read(search)
    search.close()

    if not search_result["IdList"]:
        print(f"Erorr: Taxon not found: {sterm}")
        raise Exception("Id List Empty")

    # For each taxon make a request to Entrez NCBI Taxonomy

    # Fetch taxonomy details

    for tax_id in search_result["IdList"]:
        fetch = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
        fetch_result = Entrez.read(fetch)
        fetch.close()

        #Retreive taxonomy details
        tax_data = fetch_result[0]
        lineage_dict = {rank["Rank"]: rank["ScientificName"] for rank in tax_data["LineageEx"]}
        lineage_dict["Species"] = tax_data["ScientificName"]
        lineage_dict["Taxon"] = "_".join(lineage_dict["Species"].split())
        taxonomy_data.append(lineage_dict)

    for f in taxonomy_data:
        print(f['Taxon'])

    print("--- %s seconds ---" % (time.time() - start_time))

    return taxonomy_data

def fetch_taxonomy(taxa_list, taxonomy_file, wrapper):
    taxonomy_data = {}
    # parts = len(taxa[0].split("_"))

    
    # print(taxa_list)

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

        
    taxonomy_data.columns = map(str.lower, taxonomy_data.columns)
    # for rec in records:
        # taxonomy_data[rec['Taxon']] = rec
    wrapper.append(taxonomy_data)
    return taxonomy_data

def read_taxonomy(filepath):
    print(f"Using existing taxonomy file: {filepath}")
    try:
        taxonomy_data = pd.read_csv(filepath)   

    except:
        raise Exception("Failure to read Taxonomy File")
    
    return taxonomy_data

# Root tree based on taxa listed in file
def root_tree(tree, outgroup_file):
    with open(outgroup_file) as f:
        outgroup_taxa = f.read().strip().split("\n")

    # print(outgroup_taxa) # Works

    # if set of taxa, find commmon ancestor node, else root on leaf node
    if len(outgroup_taxa) > 1:
        tree.set_outgroup(tree.get_common_ancestor(outgroup_taxa))
    else:
        tree.set_outgroup(outgroup_taxa[0])

    # return transformed tree object
    return tree

def is_monophyletic(tree, group, taxonomy, leaf_taxon):

    ranks = taxonomy.columns[taxonomy.apply(lambda col: col.str.lower().isin([group.lower()]).any())]

    # print(ranks[0])

    if ranks.empty:
        raise Exception("Taxonomic Group not Present in Tree")
    elif len(ranks) > 1:
        raise Exception("Taxonomic group appears under multiple ranks")

    group_records = taxonomy[["Taxon",ranks[0]]]

    mono_taxa = group_records[taxonomy[ranks[0]].str.lower() == group.lower()]

    mono_leaves = [leaf for leaf in leaf_taxon if leaf.rsplit('_', 1)[0] in set(mono_taxa['Taxon'])]

    isMono = tree.check_monophyly(values=mono_leaves, target_attr="name")

    return isMono

# Entry point
def main():

    # Parse arguments python monocheck [tree file] -group [taxonomic group] -outgroup [outgroup file OR outgroup taxon]
    parser = argparse.ArgumentParser(description="Check monophyly of a group in a phylogenetic tree.")
    parser.add_argument("tree", help="Newick file containing the tree")
    parser.add_argument("-clade", required=True, help="Clade to check for monophyly")
    parser.add_argument("-outgroup", required=True, help="File containing outgroup taxa (one per line)")
    parser.add_argument("-taxonomy", required=False, help="File containing taxonomy information")
    parser.add_argument("--showtree", action=argparse.BooleanOptionalAction, required=False, help="Shows the rooted tree prior to checking for Monophyly of clade")
    args = parser.parse_args()

    taxonomy=[]
    leaf_names = []
    taxa_list = []

    # create tree object from tree
    try:
        if not os.path.exists(args.tree):
            raise Exception("Tree file does not exist at path:", args.tree)
        tree = Tree(args.tree, format=1)

        leaf_names = [leaf.name for leaf in tree.get_leaves()]
        taxa_list = ["_".join((taxon[0],taxon[1])) for taxon in [taxon.split("_") for taxon in leaf_names]]
    except Exception as err:
        print("ERROR PROCESSING TREE: ", err)
        sys.exit(1)

    # Read in outgroup file
    try:
        # Root the tree
        if not os.path.exists(args.outgroup):
            raise Exception("Outgroup file does not exist at path:", args.outgroup)
        tree = root_tree(tree, args.outgroup)

        if args.showtree:
            print(tree) # Works

    except Exception as err:
        print("ERROR PROCESSING OUTGROUP FILE: ", err)
        sys.exit(1)
    
    if args.taxonomy:

        try:
            if not os.path.exists(args.taxonomy):
                raise Exception("taxonomy file does not exist at path:", args.taxonomy)
            
            taxonomy = read_taxonomy(args.taxonomy)

            if taxa_list not in taxonomy['Taxon']:
                raise Exception("Taxonomy file does not contain all the taxa in tree file\nMissing taxa:", taxa_list - taxonomy['Taxon'])

        except Exception as err:
            print("ERROR READING TAXONOMY FILE:", err)
            sys.exit(1)

    else:

        try:
            found = False
            for file in os.listdir():
                if file.endswith(".tax"):
                    if found:
                        raise Exception("Multiple taxonomy files present")
                    # print(file)
                    taxonomy = read_taxonomy(file)
                    # print(taxonomy)
                    found = True
        except Exception as err:
            print("ERROR SEARCHING FOR TAXONOMY FILE:", err)
            sys.exit(1)

        try:

            if not found:
                # Calculate new taxonomy
                wrapper = []
                ret_taxonomy = threading.Thread(name='retrieve_taxonomy', target=fetch_taxonomy, args=(taxa_list, "taxonomy.tax", wrapper))
                ret_taxonomy.start()
                while ret_taxonomy.is_alive():
                    animated_loading('Retriving taxonomy information for', len(taxa_list),'taxa...')
                if not wrapper:
                    raise Exception("No taxonomy returned by retrieval thread")
                    
                
                taxonomy = wrapper[0]

        except Exception as err:
            print("ERROR GENERATING TAXONOMY:", err)
            sys.exit(1)
            
    # # Load or generate taxonomy
    # taxonomy_file = "taxonomy.tax"
    
    # taxa_list = ["_".join((taxon[0],taxon[1])) for taxon in [taxon.split("_") for taxon in leaf_names]]

    # # taxonomy = fetch_taxonomy(taxa_list, taxonomy_file)

    # wrapper = []

    # the_process = threading.Thread(name='process', target=fetch_taxonomy, args=(taxa_list, taxonomy_file, wrapper))

    # the_process.start()

    # time.sleep(0.1)

    # while the_process.is_alive():
    #     animated_loading()

    # if not wrapper:
    #     raise Exception("Failure to retrieve Taxonomy")
    # else:
    #     taxonomy = wrapper[0]

        # Check monophyly
    result = is_monophyletic(tree, str.lower(args.clade), taxonomy, leaf_names)
    # return 1
    print(f"Is {args.clade} monophyletic? {'Yes' if result[0] else 'No'}")

    if not result[0]:
        print(f'Monophyly is prevented by the inclusion of the following leaves:')
        for leaf in result[2]:
            print(leaf)

if __name__ == "__main__":
    main()

    