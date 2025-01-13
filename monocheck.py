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
        sys.stdout.write('\r' + char + message)
        time.sleep(.1)
        sys.stdout.flush()

def request_taxonomy_data(taxa_list):
    taxonomy_data = []

    # For each taxon make a request to Entrez NCBI Taxonomy
    for taxon in taxa_list:
        try:
            # Search for the taxon in NCBI Taxonomy
            search = Entrez.esearch(db="taxonomy", term=taxon)
            search_result = Entrez.read(search)
            search.close()

            # If no results returned
            if not search_result["IdList"]:
                print(f"\rWARNING: Unable to find records on NCBI Taxonomy for: {taxon}. \nPlease check tip labels or manually insert information into taxonomy file.")
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

def fetch_taxonomy(taxa_list, taxonomy_file, wrapper):
    # Request taxonomy records from NCBI
    records = request_taxonomy_data(taxa_list)
    
    # Create a DataFrame and save to CSV
    taxonomy_data = pd.DataFrame(records).fillna("N/A")
    taxonomy_data.to_csv(taxonomy_file, index=False)

    # Place taxonomy in mutable object so that it can be accessed by the main thread upon completion
    wrapper.append(taxonomy_data)
    return 1

def read_taxonomy(filepath):
    print(f"Reading in existing taxonomy file: {filepath}")
    try:
        taxonomy_data = pd.read_csv(filepath)   
    except:
        raise Exception("Taxonomy file could not be read")
    
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
def main1():

    # Parse arguments python monocheck [tree file] -group [taxonomic group] -outgroup [outgroup file OR outgroup taxon]
    parser = argparse.ArgumentParser(description="Check monophyly of a group in a phylogenetic tree.")
    parser.add_argument("tree", help="Newick file containing the tree")
    parser.add_argument("-clade", required=True, help="Clade to check for monophyly")
    parser.add_argument("-outgroup", required=True, help="File containing outgroup taxa (one per line)")
    parser.add_argument("-taxonomy", required=False, help="File containing taxonomy information")
    parser.add_argument("--showtree", action=argparse.BooleanOptionalAction, required=False, help="Shows the rooted tree prior to checking for Monophyly of clade")
    parser.add_argument("--resettaxonomy", action=argparse.BooleanOptionalAction, required=False, help="Resets the taxonomy file by redownloading information from NCBI Taxonomy")
    args = parser.parse_args()

    taxonomy=[]
    leaf_names = []
    taxa_list = []
    trees = []


    # create tree object from tree
    try:
        if not os.path.exists(args.tree):
            raise Exception(f"Tree file does not exist at path: {args.tree}")
        tree = Tree(args.tree, format=1)

        leaf_names = [leaf.name for leaf in tree.get_leaves()]
        taxa_list = ["_".join((taxon[0],taxon[1])) for taxon in [taxon.split("_") for taxon in leaf_names]]
    except IndexError as err:
        print("ERROR PROCESSING TIP LABEL: Ensure tip labels are in format Genus_species*")
        sys.exit(1)
    except Exception as err:
        print("ERROR PROCESSING TREE: ", err.args[0])
        sys.exit(1)

    # Read in outgroup file
    try:
        # Root the tree
        if not os.path.exists(args.outgroup):
            raise Exception(f"Outgroup file does not exist at path: {args.outgroup}")
        tree = root_tree(tree, args.outgroup)

        if args.showtree:
            print(tree) # Works

    except Exception as err:
        print("ERROR PROCESSING OUTGROUP FILE: ", err.args[0])
        sys.exit(1)
    
    if args.taxonomy:

        try:
            if not os.path.exists(args.taxonomy):
                raise Exception(f"Taxonomy file does not exist at path: {args.taxonomy}")
            
            taxonomy = read_taxonomy(args.taxonomy)

        except Exception as err:
            print("ERROR READING TAXONOMY FILE:", err.args[0])
            sys.exit(1)

    else:
        found = False
        if not args.resettaxonomy:
            try:
                
                for file in os.listdir():
                    if file.endswith(".tax"):
                        if found:
                            raise Exception("Multiple taxonomy files present")
                        # print(file)
                        taxonomy = read_taxonomy(file)
                        # print(taxonomy)
                        found = True
            except Exception as err:
                print("ERROR SEARCHING FOR TAXONOMY FILE:", err.args[0])
                sys.exit(1)

        try:

            if not found:
                # Calculate new taxonomy
                wrapper = []
                start_time = time.time()
                ret_taxonomy = threading.Thread(name='retrieve_taxonomy', target=fetch_taxonomy, args=(taxa_list, "taxonomy.tax", wrapper))
                ret_taxonomy.start()
                while ret_taxonomy.is_alive():
                    animated_loading(f' Retriving taxonomy information for {len(taxa_list)} taxa...')
                
                sys.stdout.flush()
                print(f'\rRetrieving taxonomy information from NCBI Taxonomy database took {round(time.time() - start_time,1)} seconds.')

                if not wrapper:
                    raise Exception("No taxonomy returned by retrieval thread")
                    
                taxonomy = wrapper[0]

        except Exception as err:
            print("ERROR GENERATING TAXONOMY:", err.args[0])
            sys.exit(1)

    # Check no taxa are missing from taxonomy before proceeding
    try:
        missing_taxa = set(taxa_list) - set(taxonomy['Taxon'])
        if missing_taxa:
            raise Exception(f"The following taxa are missing: {missing_taxa}")
    except Exception as err:
        print("ERROR READING TAXONOMY: ", err.args[0])
        sys.exit(1)

    # Check monophyly
    try:
        result = is_monophyletic(tree, str.lower(args.clade), taxonomy, leaf_names)

        print(f"Is {args.clade} monophyletic? {'Yes' if result[0] else 'No'}")

        if not result[0]:
            print(f'Monophyly is prevented by the inclusion of the following leaves:')
            for leaf in result[2]:
                print(leaf.name)

    except Exception as err:
        print("ERROR CALCULATING MONOPHYLY:", err.args[0])
        sys.exit(1)    



def main():

    # Parse arguments python monocheck [tree file] -group [taxonomic group] -outgroup [outgroup file OR outgroup taxon]
    parser = argparse.ArgumentParser(description="Check monophyly of a group in a phylogenetic tree.")
    parser.add_argument("tree", help="Newick file containing the tree")
    parser.add_argument("-clade", required=True, help="Clade to check for monophyly")
    parser.add_argument("-outgroup", required=True, help="File containing outgroup taxa (one per line)")
    parser.add_argument("-taxonomy", required=False, help="File containing taxonomy information")
    parser.add_argument("--showtree", action=argparse.BooleanOptionalAction, required=False, help="Shows the rooted tree prior to checking for Monophyly of clade")
    parser.add_argument("--resettaxonomy", action=argparse.BooleanOptionalAction, required=False, help="Resets the taxonomy file by redownloading information from NCBI Taxonomy")
    args = parser.parse_args()

    taxonomy=[]
    leaf_names = []
    taxa_list = []
    trees = []
    paths = []

    # get trees from directory file path
    if os.path.isdir(args.tree):
        for file in os.listdir(args.tree):
            # try:
            print(file)
            trees.append(Tree(file, format=1))
            paths.append(file)
            # except Exception as err:
            #     print(f"WARNING: File {file} could not be parsed")
            #     continue
        print(f'Found {len(trees)} in directory')
        if len(trees) < 1:
            print("Exiting...")
            sys.exit(0)
    else:
        # Otherwise get the single tree given
        try:
            if not os.path.exists(args.tree):
                raise Exception(f"Tree file does not exist at path: {args.tree}")
            trees.append(Tree(args.tree, format=1))
            paths.append(args.tree)
        except Exception as err:
            print("ERROR PROCESSING TREE: ", err.args[0])
            sys.exit(1)

    # Read outgroup root file
    if not os.path.exists(args.outgroup):
        print(f"ERROR PROCESSING OUTGROUP FILE: Outgroup file does not exist at path: {args.outgroup}")
    
    for tree in trees:
        # Try to root
        try:
            tree = root_tree(tree, args.outgroup)
            if args.showtree:
                print(tree)
        except:
            print("ERROR PROCESSING OUTGROUP: Outgroup must be a file containing tip labels as seen on trees")
            sys.exit(1)

    # Base tree for taxonomy work
    tree = trees[0]

    # Get tips and leaves
    try:
        leaf_names = [leaf.name for leaf in tree.get_leaves()]
        taxa_list = ["_".join((taxon[0],taxon[1])) for taxon in [taxon.split("_") for taxon in leaf_names]] 
    except IndexError as err:
        print("ERROR PROCESSING TIP LABEL: Ensure tip labels are in format Genus_species*")
        sys.exit(1)

    # Begin taxonomy process
    if args.taxonomy:

        try:
            if not os.path.exists(args.taxonomy):
                raise Exception(f"Taxonomy file does not exist at path: {args.taxonomy}")
            
            taxonomy = read_taxonomy(args.taxonomy)

        except Exception as err:
            print("ERROR READING TAXONOMY FILE:", err.args[0])
            sys.exit(1)

    else:
        found = False
        if not args.resettaxonomy:
            try:
                
                for file in os.listdir():
                    if file.endswith(".tax"):
                        if found:
                            raise Exception("Multiple taxonomy files present")
                        # print(file)
                        taxonomy = read_taxonomy(file)
                        # print(taxonomy)
                        found = True
            except Exception as err:
                print("ERROR SEARCHING FOR TAXONOMY FILE:", err.args[0])
                sys.exit(1)

        try:

            if not found:
                # Calculate new taxonomy
                wrapper = []
                start_time = time.time()
                if len(trees) > 1:
                    ret_taxonomy = threading.Thread(name='retrieve_taxonomy', target=fetch_taxonomy, args=(taxa_list, args.tree + "/taxonomy.tax", wrapper))
                else:
                    ret_taxonomy = threading.Thread(name='retrieve_taxonomy', target=fetch_taxonomy, args=(taxa_list,"/taxonomy.tax", wrapper))
                ret_taxonomy.start()
                while ret_taxonomy.is_alive():
                    animated_loading(f' Retriving taxonomy information for {len(taxa_list)} taxa...')
                
                sys.stdout.flush()
                print(f'\rRetrieving taxonomy information from NCBI Taxonomy database took {round(time.time() - start_time,1)} seconds.')

                if not wrapper:
                    raise Exception("No taxonomy returned by retrieval thread")
                    
                taxonomy = wrapper[0]

        except Exception as err:
            print("ERROR GENERATING TAXONOMY:", err.args[0])
            sys.exit(1)

    # Check no taxa are missing from taxonomy before proceeding
    try:
        missing_taxa = set(taxa_list) - set(taxonomy['Taxon'])
        if missing_taxa:
            raise Exception(f"The following taxa are missing: {missing_taxa}")
    except Exception as err:
        print("ERROR READING TAXONOMY: ", err.args[0])
        sys.exit(1)

    
    # Check monophyly

    for t in range(len(trees)):
        try:
            tree_leaves = [leaf.name for leaf in trees[t].get_leaves()]
            if tree_leaves != leaf_names:
                raise Exception(f"Leaf names are not consistent. Missing expected tips: {tree_leaves - leaf_names}")
            result = is_monophyletic(trees[t], str.lower(args.clade), taxonomy, leaf_names)

            print(f"{paths[0]}>>> Is {args.clade} monophyletic? {'Yes' if result[0] else 'No'}")

            if not result[0]:
                print(f'Monophyly is prevented by the inclusion of the following leaves:')
                for leaf in result[2]:
                    print(leaf.name)

        except Exception as err:
            print(f"ERROR CALCULATING MONOPHYLY IN TREE {paths[t]}: {err.args[0]}")

    
    

if __name__ == "__main__":
    main()