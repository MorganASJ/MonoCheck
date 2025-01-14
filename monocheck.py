from Bio import Entrez
import pandas as pd
from ete3 import Tree
import argparse, os, sys, time, threading

Entrez.email = ""

#Spinner while taxonomy calculates
def animated_loading(message):
    chars = "/—\\|" 
    for char in chars:
        sys.stdout.write('\r' + char + message)
        time.sleep(.1)
        sys.stdout.flush()

# Fetch taxonomy data from NCBI
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

# Fetch taxonomy data from NCBI, convert to DataFrame and save to CSV
def fetch_taxonomy(taxa_list, taxonomy_file, wrapper):
    # Request taxonomy records from NCBI
    records = request_taxonomy_data(taxa_list)
    
    # Create a DataFrame and save to CSV
    taxonomy_data = pd.DataFrame(records).fillna("N/A")
    taxonomy_data.to_csv(taxonomy_file, index=False)

    # Place taxonomy in mutable object so that it can be accessed by the main thread upon completion
    wrapper.append(taxonomy_data)
    return 1

# Read in existing taxonomy file from path
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

    leaf_names = [leaf.name for leaf in tree.get_leaves()]

    outgroup_taxa = list(set(leaf_names).intersection(outgroup_taxa))

    # if set of taxa, find commmon ancestor node, else root on leaf node
    if len(outgroup_taxa) > 1:
        tree.set_outgroup(tree.get_common_ancestor(outgroup_taxa))
    else:
        tree.set_outgroup(outgroup_taxa[0])

    # return transformed tree object
    return tree

# Checks if a clade is mono or paraphyletic in a given tree
def is_monophyletic(tree, group, taxonomy, leaf_taxon):

    # Get the column which features the group/clade being assessed, and set to lower case
    ranks = taxonomy.columns[taxonomy.apply(lambda col: col.str.lower().isin([group.lower()]).any())]

    # Check if group is present in taxonomy
    if ranks.empty:
        raise Exception("Taxonomic Group not present in taxonomy file")
    elif len(ranks) > 1:
        raise Exception("Taxonomic group appears under multiple ranks")

    # Gets Taxon and Rank columns
    group_records = taxonomy[["Taxon",ranks[0]]]

    # Gets taxa with clade membership for that rank
    mono_taxa = group_records[taxonomy[ranks[0]].str.lower() == group.lower()]

    # Get the leaves that are in the clade
    mono_leaves = [leaf for leaf in leaf_taxon if "_".join(leaf.split("_")[:2]) in set(mono_taxa['Taxon'])]

    # Use ete3 function to check monophyly of those leaves
    return tree.check_monophyly(values=mono_leaves, target_attr="name", ignore_missing=True)

# Main function
def main():

    # Parse arguments python monocheck [tree file] -group [taxonomic group] -outgroup [outgroup file OR outgroup taxon]
    parser = argparse.ArgumentParser(description="Check monophyly of a group in a phylogenetic tree.")
    parser.add_argument("tree", help="Newick file containing the tree")
    parser.add_argument("-clade", nargs='+', required=True, help="Clade(s) to check for monophyly")
    parser.add_argument("-outgroup", required=True, help="File containing outgroup taxa (one per line)")
    parser.add_argument("-email", required=True, help="Email for")
    parser.add_argument("-taxonomy", required=False, help="File containing taxonomy information")
    parser.add_argument("-out", required=False, help="Write infomration about monophyly to a csv file")
    parser.add_argument("--showtree", action=argparse.BooleanOptionalAction, required=False, help="Shows the rooted tree prior to checking for Monophyly of clade")
    parser.add_argument("--resettaxonomy", action=argparse.BooleanOptionalAction, required=False, help="Resets the taxonomy file by redownloading information from NCBI Taxonomy")
    args = parser.parse_args()

    # Declare variables
    taxonomy=[]
    leaf_names = []
    taxa_list = []
    trees = []
    paths = []

    # Set email for Entrez
    Entrez.email = args.email
    
    # Normalise the tree path
    args.tree = os.path.normpath(args.tree)

    # get trees from directory file path
    if os.path.isdir(args.tree):
        for file in os.listdir(args.tree):
            try:
                print('Reading',file,'...')
                path = args.tree + '/' + file
                trees.append(Tree(path, format=1))
                paths.append(path)
            except Exception as err:
                print(f"WARNING: File {file} could not be parsed as a tree file")
                continue
        print(f'Found {len(trees)} tree file in directory')
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
        # Check if outgroup file is in the same directory as the tree
        if os.path.exists(args.tree + "/" + args.outgroup):
            args.outgroup = args.tree + "/" + args.outgroup
        else:
            print(f"ERROR PROCESSING OUTGROUP FILE: Outgroup file does not exist at path: {args.outgroup}")
            sys.exit(1)
    
    # Root the trees
    for i in range(len(trees)):
        # Try to root
        try:
            # Use ete3 function to root tree
            tree = root_tree(trees[i], args.outgroup)
            if args.showtree:
                print(paths[i])
                print(trees[i])
        except Exception as err:
            print("ERROR PROCESSING OUTGROUP: Outgroup must be a file containing tip labels as seen on trees")
            print(err)
            sys.exit(1)

    # Base tree for taxonomy work - always uses first tree. If users want to use multiple trees, they should ensure the same taxa are present in each tree.
    tree = trees[0]

    # Get tips and leaves, probably should replace this and calculate for each tree when needed. For now leaving as is.
    try:
        leaf_names = [leaf.name for leaf in tree.get_leaves()]
        taxa_list = ["_".join((taxon[0],taxon[1])) for taxon in [taxon.split("_") for taxon in leaf_names]] 
    except IndexError as err:
        print("ERROR PROCESSING TIP LABEL: Ensure tip labels are in format Genus_species_*")
        sys.exit(1)

    # Begin taxonomy process

    # Check if taxonomy file is provided in user input
    if args.taxonomy:
        
        # Check if the taxonomy file exists in working directory or tree directory
        try:
            if not os.path.exists(args.taxonomy):
                if os.path.exists(args.tree + '/' + args.taxonomy):
                    args.taxonomy = args.tree + '/' + args.taxonomy
                else:
                    raise Exception(f"Taxonomy file does not exist at path: {args.taxonomy}")
            
            # Read in the taxonomy file
            taxonomy = read_taxonomy(args.taxonomy)

        except Exception as err:
            print("ERROR READING TAXONOMY FILE:", err.args[0])
            sys.exit(1)

    # If no taxonomy file provided
    else:
        found = False
        # If not reseting taxonomy, check for one locally
        if not args.resettaxonomy:
            try:   
                # Check for taxonomy file in working directory or tree directory
                if os.path.isdir(args.tree):
                    path = args.tree + '/'
                else:
                    path = '.'
                for file in os.listdir(path):
                    if file.endswith(".tax"):
                        # Ensuring there is no confusion with multiple taxonomy files
                        if found:
                            raise Exception("Multiple taxonomy files present")
                        taxonomy = read_taxonomy(path + '/' + file)
                        found = True
            except Exception as err:
                print("ERROR SEARCHING FOR TAXONOMY FILE:", err.args[0])
                sys.exit(1)

        # If no taxonomy file found, or resetting taxonomy, fetch from NCBI
        try:
            if not found:

                # Calculate new taxonomy
                wrapper = []
                # Timer for user information
                start_time = time.time()
                # Start thread to fetch taxonomy
                if os.path.isdir(args.tree):
                    ret_taxonomy = threading.Thread(name='retrieve_taxonomy', target=fetch_taxonomy, args=(taxa_list, os.path.join(args.tree, 'taxonomy.tax'), wrapper))
                else:
                    ret_taxonomy = threading.Thread(name='retrieve_taxonomy', target=fetch_taxonomy, args=(taxa_list,"./taxonomy.tax", wrapper))
                ret_taxonomy.start()

                # Spinner while taxonomy is being calculated in main thread
                while ret_taxonomy.is_alive():
                    animated_loading(f' Retriving taxonomy information for {len(taxa_list)} taxa...')
                
                # Once thread is complete clear the spinner
                sys.stdout.flush()

                # Of main thread has not received taxonomy data, raise exception
                if not wrapper:
                    raise Exception("No taxonomy returned by retrieval thread")
                
                # Inform user how much time it took to retreive taxonomy
                print(f'\rRetrieving taxonomy information from NCBI Taxonomy database took {round(time.time() - start_time,1)} seconds.')
                
                # Set taxonomy to the returned data
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
    print("\n----- Checking Monophyly of clades -----")

    # Collect data rows in a list
    output_rows = []

    # For each tree, check monophyly of the clade(s)
    for t in range(len(trees)):
        # if there is more than one clade to check, we will output the results in a clearer format
        if len(args.clade) > 1:
            for clade in args.clade:
                try:
                    result = is_monophyletic(trees[t], str.lower(clade), taxonomy, leaf_names)
                    if result[0]:
                        print(f"{paths[t]}>>> {clade} | {result[1]} | *")
                    else:
                        print(f"{paths[t]}>>> {clade} | {result[1]} | {len(result[2])}")

                    # Append data row to the list
                    output_rows.append({
                        'Tree': paths[t],
                        'Clade': clade,
                        'Status': result[1],
                        'NumberDisrupting': len(result[2]),
                        'DisruptiveTaxa': [leaf.name for leaf in result[2]]
                    })

                except Exception as err:
                    print(f"{paths[t]}>>> {clade} | ERROR | Clade not in tree")
                    print(err)
        # Single clade to check we just say yes or no
        else:
            result = is_monophyletic(trees[t], str.lower(args.clade[0]), taxonomy, leaf_names)
            print(f"{paths[0]}>>> Is {args.clade[0]} monophyletic? {'Yes' if result[0] else 'No'}")

            if not result[0]:
                print(f'Monophyly is prevented by the inclusion of the following leaves:')
                for leaf in result[2]:
                    print(leaf.name)

            # Append data row to the list
            output_rows.append({
                'Tree': paths[0],
                'Clade': args.clade[0],
                'Status': result[1],
                'NumberDisrupting': len(result[2]),
                'DisruptiveTaxa': [leaf.name for leaf in result[2]]
            })

        print("----------------------------------------")

    # Create the DataFrame from the collected rows
    output = pd.DataFrame(output_rows, columns=['Tree', 'Clade', 'Status', 'NumberDisrupting', 'DisruptiveTaxa'])

    # Write output to CSV
    try:
        if args.out:
            print("Writing monophyly information to CSV file at:", args.out)
            output.to_csv(args.out, index=False)
    except Exception as e:
        print("ERROR WRITING OUTPUT CSV: Unable to write to file path", args.out)
        print(e)
        exit(1)

    # End program
    
# Run the main function
if __name__ == "__main__":
    main()